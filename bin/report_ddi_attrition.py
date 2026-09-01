#!/usr/bin/env python3
"""One report answering "where did the DDIs of source X go?".

``ppi-splitting``'s own DDI attrition waterfall starts at the CSV
``EXPORT_SPLIT_DDIS`` handed it, so it structurally cannot show the DDIs that
never reached the splitting stage -- and under ``instance_tier =
human_reviewed`` that is the large majority of 3did.  A run where 3did contributes ~20 000 pairs
and the splitting stage reports ~3 200 input DDIs looks like 17 000 DDIs
vanished with no accounting anywhere.  They are accounted for, but the counters
were spread over three ``.command.log`` files that go away with the work
directory.  This collects them into one published TSV.

One row per source **token**.  ``domain_domain_interaction.source`` is a
comma-joined list (``"PPIDM,PPIDM_Gold"``), so a DDI contributed by both PPIDM
and its Gold class is counted under each token; token rows therefore overlap and
do not sum to the ``ALL_SOURCES`` row, which counts distinct DDI rows.

Columns, left to right in the order a DDI travels:

``offered``
    Contributions this source made before any deduplication -- for 3did the
    pairs read out of its SQLite (from ``INSERT_3DID``'s counts file), for an
    external source the rows of its normalized TSV.  ``NA`` for
    ``sampled_negative``, whose contributions arrive as ppi-splitting CSVs and
    are only ever counted after dedup.
``in_db``
    DDI rows carrying this token after insert, i.e. ``pruned + surviving``.
    ``offered - in_db`` is dedup plus the pairs a protected source already owned
    plus the label conflicts in the ``conflict`` column.
``to_splitting``
    DDIs handed to ppi-splitting.  Only the split population has one -- the
    ``--splitting-source`` token, 3did -- everything else is ``NA`` by design:
    external sources are the held-out evaluation set and are never split.
``pruned``
    Deleted by ``PRUNE_UNREPRESENTED_DDIS`` because a family reached zero domain
    instances.  This is where ``instance_tier`` shows up.
``surviving``
    DDI rows in the published master.
``conflict``
    Contributions dropped or deleted by the external-source label-conflict rule
    (``source_conflicts.tsv``).  Not a subset of any other column: these pairs
    never entered the database.
"""

import argparse
import csv
import sqlite3
import sys
from collections import Counter

from ddi_db_utils import split_sources

COLUMNS = ["source", "offered", "in_db", "to_splitting", "pruned", "surviving", "conflict"]

#: The second report: what the widened protein universe actually bought.
#:
#: ``--instance_tier`` exists to stop discarding DDIs whose families have no human
#: Swiss-Prot representative, and the only way to answer "was that worth it" is to
#: count the survivors by the stratum they rest on. Same four strata as upstream's
#: ``TIERS``, same order -- reviewed outranks human.
TIER_COLUMNS = ["scope", "tier", "count"]

TIERS = ("human_reviewed", "other_reviewed", "human_unreviewed", "other_unreviewed")

HUMAN_TAXON = "9606"

#: A family with no instance at all, and a DDI touching one. Should not survive
#: PRUNE_UNREPRESENTED_DDIS; carried so the report cannot silently lose a row.
NO_INSTANCE = "no_instance"

#: `protein.reviewed` is NULL on a database written before the column existed.
UNKNOWN_TIER = "unknown_review_status"

TOTAL_ROW = "ALL_SOURCES"

NA = "NA"


def tier_of(taxon_id, reviewed):
    """Which of upstream's four strata one protein belongs to.

    Mirrors ``fetch_domains.py``'s ``tier_of()`` deliberately: if the two ever
    disagree, this report would describe a stratification the sampler did not use.
    """
    flag = (reviewed or "").strip().lower()
    if flag not in ("reviewed", "unreviewed"):
        return UNKNOWN_TIER
    is_human = (taxon_id or "").strip() == HUMAN_TAXON
    if flag == "reviewed":
        return TIERS[0] if is_human else TIERS[1]
    return TIERS[2] if is_human else TIERS[3]


def tier_rank(tier):
    """Sort key: the four strata in fill order, then the two escape hatches."""
    return TIERS.index(tier) if tier in TIERS else len(TIERS) + (tier == NO_INSTANCE)


def tier_breakdown(conn):
    """``(rows, best_by_family)``: per-tier instance/protein/family/DDI counts.

    A family's tier is the *best* one any of its instances reaches, and a DDI's is
    the worse of its two families' -- so ``max(best_a, best_b)``. That is exactly
    the "would this DDI have survived a narrower universe" question: a DDI above
    ``human_reviewed`` has at least one family with no human Swiss-Prot instance,
    which under ``--instance_tier human_reviewed`` would have reached zero
    instances and taken the DDI with it in PRUNE_UNREPRESENTED_DDIS.
    """
    instances = Counter()
    proteins = {}
    best_by_family = {}
    for pfam, uniprot, taxon, reviewed in conn.execute(
        "SELECT d.pfam_id, p.uniprot_id, p.taxon_id, p.reviewed "
        "FROM domain_protein_map AS dpm "
        "JOIN domain AS d ON d.id = dpm.domain_id "
        "JOIN protein AS p ON p.id = dpm.protein_id "
        "WHERE dpm.instance_id IS NOT NULL"
    ):
        tier = tier_of(taxon, reviewed)
        instances[tier] += 1
        proteins[uniprot] = tier
        if pfam not in best_by_family or tier_rank(tier) < tier_rank(best_by_family[pfam]):
            best_by_family[pfam] = tier

    ddis = Counter()
    for pfam_a, pfam_b in conn.execute(
        "SELECT da.pfam_id, db.pfam_id FROM domain_domain_interaction AS ddi "
        "JOIN domain AS da ON da.id = ddi.domain_id_a "
        "JOIN domain AS db ON db.id = ddi.domain_id_b"
    ):
        tier_a = best_by_family.get(pfam_a, NO_INSTANCE)
        tier_b = best_by_family.get(pfam_b, NO_INSTANCE)
        ddis[tier_a if tier_rank(tier_a) >= tier_rank(tier_b) else tier_b] += 1

    scopes = {
        "instances": instances,
        "proteins": Counter(proteins.values()),
        "families": Counter(best_by_family.values()),
        "ddis": ddis,
    }
    rows = []
    for scope in ("ddis", "families", "proteins", "instances"):
        counts = scopes[scope]
        for tier in sorted(counts, key=tier_rank):
            rows.append({"scope": scope, "tier": tier, "count": counts[tier]})
    return rows, ddis


def write_tier_report(path, rows, ddis):
    with open(path, "w", newline="") as fh:
        for line in TIER_DOC.strip().splitlines():
            fh.write(f"# {line}".rstrip() + "\n")
        writer = csv.DictWriter(fh, fieldnames=TIER_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    beyond = sum(n for tier, n in ddis.items() if tier != TIERS[0])
    total = sum(ddis.values())
    print(
        f"[attrition] tiers: {total} surviving DDIs, {beyond} of them resting on a "
        "stratum beyond human_reviewed -- those are the DDIs a human-reviewed-only "
        "universe would have pruned",
        flush=True,
    )
    for row in rows:
        print("[attrition] tier {scope}/{tier} = {count}".format(**row), flush=True)


TIER_DOC = """What the run's protein universe bought, per stratum.

`scope` is what is being counted: `instances` and `proteins` are the domain
instances and their parent proteins; `families` counts each Pfam family under the
*best* stratum any of its instances reaches; `ddis` counts each surviving DDI
under the *worse* of its two families' best strata.

That last one is the number `--instance_tier` exists for. A DDI in any row but
`human_reviewed` has at least one family with no human Swiss-Prot instance, so a
`--instance_tier human_reviewed` run would have left that family with zero
instances and PRUNE_UNREPRESENTED_DDIS would have deleted the DDI.

Strata are upstream's, in fill order -- reviewed outranks human, so a family with
no human Swiss-Prot member takes a curated non-human sequence before an
auto-annotated human one.

`no_instance` should not appear: PRUNE_UNREPRESENTED_DDIS deletes those DDIs.
`unknown_review_status` means `protein.reviewed` is NULL.
"""


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--db", required=True, help="pruned master DB; opened read-only")
    p.add_argument("--pruned-ddis", required=True, help="PRUNE_UNREPRESENTED_DDIS' pruned_ddis.tsv")
    p.add_argument("--split-ddis", required=True, help="EXPORT_SPLIT_DDIS' split_ddis.csv")
    p.add_argument(
        "--external-ddis",
        nargs="*",
        default=[],
        help="normalized external-source TSVs (bin/external_ddi_tsv.py), for the offered column",
    )
    p.add_argument(
        "--offered-counts",
        nargs="*",
        default=[],
        help="counts TSVs from a protected source's inserter (columns: source, offered, ...)",
    )
    p.add_argument("--conflicts", help="INSERT_EXTERNAL_SOURCES' source_conflicts.tsv")
    p.add_argument(
        "--splitting-source",
        default="3did",
        help="the token whose DDIs form the split population (default: 3did)",
    )
    p.add_argument("--out", required=True, help="ddi_source_attrition.tsv")
    p.add_argument("--tier-out", required=True, help="ddi_tier_breakdown.tsv")
    p.add_argument("--versions", required=True)
    p.add_argument("--process-name", required=True)
    return p.parse_args()


def read_tsv_rows(path):
    """Yield dict rows of a tab-separated file with a header. Missing file -> nothing."""
    try:
        fh = open(path, newline="")
    except FileNotFoundError:
        return
    with fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            yield row


def count_tokens(rows, column="source"):
    """``(per_token, n_rows)`` -- token counts and the number of rows seen.

    A row whose ``source`` holds several tokens increments each of them once, so
    the token counts overlap while ``n_rows`` stays a count of distinct rows.
    """
    per_token, n_rows = Counter(), 0
    for row in rows:
        value = (row.get(column) or "").strip()
        if not value:
            continue
        n_rows += 1
        for token in split_sources(value):
            per_token[token] += 1
    return per_token, n_rows


def surviving_tokens(conn):
    """Token counts and total row count over the master's DDI table."""
    rows = ({"source": src} for (src,) in conn.execute("SELECT source FROM domain_domain_interaction"))
    return count_tokens(rows)


def offered_counts(external_paths, counts_paths):
    """``{token: offered}`` from the external TSVs and any inserter counts files."""
    offered = Counter()
    for path in external_paths:
        per_token, _n = count_tokens(read_tsv_rows(path))
        offered.update(per_token)
    for path in counts_paths:
        for row in read_tsv_rows(path):
            source = (row.get("source") or "").strip()
            value = (row.get("offered") or "").strip()
            if source and value:
                offered[source] += int(value)
    return offered


def count_csv_rows(path):
    """Data rows of a CSV with one header line."""
    with open(path, newline="") as fh:
        return max(sum(1 for _ in csv.reader(fh)) - 1, 0)


def build_rows(counts):
    """Assemble the report rows, ordered by surviving count then name."""
    tokens = sorted(
        counts["tokens"],
        key=lambda t: (-counts["surviving"].get(t, 0), -counts["in_db"].get(t, 0), t),
    )
    out = []
    for token in tokens:
        pruned = counts["pruned"].get(token, 0)
        surviving = counts["surviving"].get(token, 0)
        out.append(
            {
                "source": token,
                "offered": counts["offered"].get(token, NA),
                "in_db": pruned + surviving,
                "to_splitting": counts["to_splitting"] if token == counts["splitting_source"] else NA,
                "pruned": pruned,
                "surviving": surviving,
                "conflict": counts["conflict"].get(token, 0),
            }
        )
    out.append(
        {
            "source": TOTAL_ROW,
            # A row sum, not a distinct count: one pair offered by two sources is
            # two contributions. Every other column on this row is distinct DDIs.
            "offered": counts["offered_rows"] or NA,
            "in_db": counts["pruned_rows"] + counts["surviving_rows"],
            "to_splitting": counts["to_splitting"],
            "pruned": counts["pruned_rows"],
            "surviving": counts["surviving_rows"],
            "conflict": counts["conflict_rows"],
        }
    )
    return out


def write_report(path, rows):
    with open(path, "w", newline="") as fh:
        for line in __doc__.strip().splitlines():
            fh.write(f"# {line}".rstrip() + "\n")
        writer = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()

    conn = sqlite3.connect(f"file:{args.db}?mode=ro", uri=True)
    surviving, surviving_rows = surviving_tokens(conn)
    tier_rows, tier_ddis = tier_breakdown(conn)
    conn.close()

    pruned, pruned_rows = count_tokens(read_tsv_rows(args.pruned_ddis))
    conflict, conflict_rows = (Counter(), 0)
    if args.conflicts:
        conflict, conflict_rows = count_tokens(read_tsv_rows(args.conflicts))

    offered = offered_counts(args.external_ddis, args.offered_counts)
    offered_rows = sum(offered.values())

    in_db = Counter(pruned)
    in_db.update(surviving)

    counts = {
        "tokens": set(surviving) | set(pruned) | set(conflict) | set(offered),
        "offered": offered,
        "offered_rows": offered_rows,
        "in_db": in_db,
        "to_splitting": count_csv_rows(args.split_ddis),
        "pruned": pruned,
        "pruned_rows": pruned_rows,
        "surviving": surviving,
        "surviving_rows": surviving_rows,
        "conflict": conflict,
        "conflict_rows": conflict_rows,
        "splitting_source": args.splitting_source,
    }

    rows = build_rows(counts)
    write_report(args.out, rows)
    write_tier_report(args.tier_out, tier_rows, tier_ddis)

    if args.splitting_source not in surviving:
        print(
            f"[attrition] WARNING: --splitting-source {args.splitting_source} has no surviving DDI "
            "in the master; the to_splitting column is attached to a source that is not there",
            file=sys.stderr,
        )
    for row in rows:
        print(
            "[attrition] {source}: offered={offered} in_db={in_db} to_splitting={to_splitting} "
            "pruned={pruned} surviving={surviving} conflict={conflict}".format(**row),
            flush=True,
        )

    with open(args.versions, "w") as f:
        f.write(f'"{args.process_name}":\n')
        f.write(f"    python: {sys.version.split()[0]}\n")
        f.write(f"    sqlite3: {sqlite3.sqlite_version}\n")


if __name__ == "__main__":
    main()
