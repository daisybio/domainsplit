#!/usr/bin/env python3
"""Groovy-level checks on the `instance_tier` derivation (needs `nextflow`).

`instance_tier` is the only protein-universe knob a user sets; three more values
are derived from it *in `nextflow.config`* rather than in the workflow, because
FETCH_DOMAIN_META lives in the read-only submodule and reads
`params.instance_tiers` straight out of its own script block. That derivation is
config-level Groovy, so it is invisible to every other test in this directory --
and it depends on a Nextflow behaviour worth pinning: a command-line
`--instance_tier X` IS visible to the config, because Nextflow merges CLI params
into the binding before parsing it. If that ever stops being true, every derived
value silently reverts to the default and a run asks for a universe nobody
requested.

Driven through `nextflow run -preview`, which compiles the DAG and runs the
workflow body without submitting a task, and asserted against the one-line
universe summary the guard block logs.

**These currently stop at the pinned submodule's own guard** (`--instance_tier
must be 'any' or 'human_only'`), which is the pre-bump state: the derivation and
this repo's guards have already run and logged by then. The assertions are
written against the log line and the error text, not the exit status, so they
keep working once the pin moves.

Skipped when `nextflow` is not on PATH.
"""

import os
import re
import shutil
import subprocess

import pytest

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

pytestmark = pytest.mark.skipif(
    shutil.which("nextflow") is None, reason="nextflow is not installed"
)

#: Our four-valued knob -> the stratum list ppi-splitting is given. The one place
#: the two vocabularies meet, restated here so a change to it is deliberate.
EXPECTED_TIERS = {
    "human_reviewed": "human_reviewed",
    "all_species_reviewed": "human_reviewed,other_reviewed",
    "human_any_review_status": "human_reviewed,human_unreviewed",
}

EXPECTED_TAXA = {
    "human_reviewed": "9606",
    "all_species_reviewed": "all species",
    "human_any_review_status": "9606",
}

UNIVERSE_RE = re.compile(
    r"protein universe: instance_tier=(?P<tier>\S+) -> tiers=(?P<tiers>\S+), "
    r"taxa=(?P<taxa>[^,]+), (?P<n>\d+) UniProt flat file"
)


def preview(tmp_path, *args):
    env = dict(os.environ, NXF_OFFLINE="1")
    # Any existing file: without it the preview stops at the absent local
    # gurobi.lic, and a clean exit status would mean nothing.
    return subprocess.run(
        ["nextflow", "run", REPO, "-profile", "test", "-preview",
         "--outdir", str(tmp_path), "--gurobi_license", os.path.join(REPO, "nextflow.config"), *args],
        cwd=REPO, env=env, capture_output=True, text=True, timeout=600,
    )


@pytest.mark.parametrize("tier", sorted(EXPECTED_TIERS))
def test_tier_resolves_to_the_right_universe(tier, tmp_path):
    proc = preview(tmp_path, "--instance_tier", tier)
    match = UNIVERSE_RE.search(proc.stdout + proc.stderr)
    assert match, proc.stdout + proc.stderr

    assert match["tier"] == tier
    assert match["tiers"] == EXPECTED_TIERS[tier]
    assert match["taxa"] == EXPECTED_TAXA[tier]
    # -profile test pins one fixture flat file for every tier, so the count is
    # the profile's, not the derived default's. What matters here is that the
    # override survived the derivation instead of being clobbered by it.
    assert match["n"] == "1"


def test_unimplemented_tier_fails_before_any_task(tmp_path):
    """`all_species_any_review_status` needs the 110 GB TrEMBL file, which
    FETCH_DOMAIN_META parses into an in-memory dict. It has to fail in the guard
    block, not at the end of a 4.7 GB regions stream."""
    proc = preview(tmp_path, "--instance_tier", "all_species_any_review_status")
    out = proc.stdout + proc.stderr

    assert proc.returncode != 0
    assert "all_species_any_review_status is not implemented" in out, out
    assert "daisybio/domainsplit/issues/4" in out, out
    # Nothing downstream may have run: the point is that it costs nothing.
    assert not UNIVERSE_RE.search(out), out


def test_retired_vocabulary_is_rejected_with_the_mapping(tmp_path):
    """`any` used to mean a Swiss-Prot universe over every species -- the new
    `all_species_reviewed`, NOT `all_species_any_review_status`. Aliasing it would
    silently change the universe of a working invocation, so it errors with the
    mapping instead."""
    proc = preview(tmp_path, "--instance_tier", "any")
    out = proc.stdout + proc.stderr

    assert proc.returncode != 0
    assert "any -> all_species_reviewed" in out, out
    assert "human_only -> human_reviewed" in out, out


# --external_test_sources decides which external sources are fetched at all: a
# source it leaves out must not appear in the DAG (a Negatome download failing
# on a run that never uses Negatome is the bug this pins).
SOURCE_PROCESSES = {"DOWNLOAD_NEGATOME", "PARSE_NEGATOME", "PARSE_SINGLE_DOMAIN_PPI",
                    "PARSE_PPIDM", "PARSE_INVERSE_PPIDM",
                    "EXTRACT_PPI_MITAB", "BUILD_PPIDM_SOURCES", "COMPUTE_PPIDM"}
PPIDM_CHAIN = {"EXTRACT_PPI_MITAB", "BUILD_PPIDM_SOURCES", "COMPUTE_PPIDM"}


def dag_processes(tmp_path, *args):
    dag = tmp_path / "dag.mmd"
    proc = preview(tmp_path, "-with-dag", str(dag), *args)
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert dag.exists(), "preview wrote no DAG"
    text = dag.read_text()
    return {p for p in SOURCE_PROCESSES if p in text}


# None = the pipeline default (PPIDM_Gold,INVERSE_PPIDM_Gold); "" = every source,
# through a params file and through `--external_test_sources=` (a bare
# `--external_test_sources ''` arrives as the String "true"). -profile test's
# url_negatome is a local file, so it is staged straight into PARSE_NEGATOME and
# DOWNLOAD_NEGATOME (URLs only) must not appear.
@pytest.mark.parametrize("sources, expected", [
    (None, PPIDM_CHAIN | {"PARSE_PPIDM", "PARSE_INVERSE_PPIDM"}),
    ("", SOURCE_PROCESSES - {"DOWNLOAD_NEGATOME"}),
    ("=", SOURCE_PROCESSES - {"DOWNLOAD_NEGATOME"}),
    ("PPIDM_Gold", PPIDM_CHAIN | {"PARSE_PPIDM"}),
    ("INVERSE_PPIDM_Silver", PPIDM_CHAIN | {"PARSE_INVERSE_PPIDM"}),
    ("PPIDM_Gold,negatome", PPIDM_CHAIN | {"PARSE_PPIDM", "PARSE_NEGATOME"}),
    ("negatome", {"PARSE_NEGATOME"}),
])
def test_external_test_sources_prunes_fetches(sources, expected, tmp_path):
    if sources is None:
        args = []
    elif sources == "":
        pf = tmp_path / "params.yml"
        pf.write_text('external_test_sources: ""\n')
        args = ["-params-file", str(pf)]
    elif sources == "=":
        args = ["--external_test_sources="]
    else:
        # `none`: a one-label selection is refused at startup under `uniform`
        args = ["--external_test_sources", sources, "--external_test_balance", "none"]
    assert dag_processes(tmp_path, *args) == expected


@pytest.mark.parametrize("sources", ["PPIDM_Gold", "single_domain_ppi,PPIDM", "negatome,INVERSE_PPIDM_Gold"])
def test_uniform_balance_refuses_a_one_label_selection(sources, tmp_path):
    proc = preview(tmp_path, "--external_test_sources", sources, "--external_test_balance", "uniform")
    out = proc.stdout + proc.stderr
    assert proc.returncode != 0
    assert "--external_test_balance uniform needs both labels" in out, out


def test_inverse_ppidm_refused_with_legacy_ppidm_tsv(tmp_path):
    proc = preview(tmp_path, "--ppidm_tsv", str(tmp_path / "x.tsv"))
    assert proc.returncode != 0
    assert "--ppidm_tsv is set" in proc.stdout + proc.stderr


def test_external_test_sources_rejects_unknown_token(tmp_path):
    proc = preview(tmp_path, "--external_test_sources", "PPIDM_gold")
    assert proc.returncode != 0
    assert "unknown source PPIDM_gold" in proc.stdout + proc.stderr


@pytest.mark.parametrize("name", ["structure_data_methods", "embedding_models_domainbench",
                                  "external_test_sources", "swissprot_taxon_ids"])
def test_bare_empty_off_switch_fails_with_the_working_form(name, tmp_path):
    proc = preview(tmp_path, f"--{name}", "")
    out = proc.stdout + proc.stderr
    assert proc.returncode != 0
    assert f"`--{name}=`" in out, out


def test_default_preview_exits_cleanly(tmp_path):
    """-profile test gives most `url_*` as local paths, several of them the same
    file (url_hippie is also a PPIDM input): the manifest's local-file staging
    must handle both. `.values().unique()` on the map threw
    UnsupportedOperationException here."""
    proc = preview(tmp_path)
    assert proc.returncode == 0, proc.stdout + proc.stderr


# Which classes PARSE_PPIDM / PARSE_INVERSE_PPIDM are handed (logged at startup).
# Explicit tokens narrow --ppidm_classes only when every token is a PPIDM one:
# next to negatome or single_domain_ppi an unparsed class would hide a label
# conflict from INSERT_EXTERNAL_SOURCES.
@pytest.mark.parametrize("args, expected", [
    ([], {"PPIDM": "Gold", "INVERSE_PPIDM": "Gold"}),
    (["--ppidm_classes", "BRONZE,GOLD"], {"PPIDM": "GOLD", "INVERSE_PPIDM": "GOLD"}),
    (["--external_test_sources", "PPIDM_Gold,negatome"], {"PPIDM": "Bronze,Silver,Gold"}),
    (["--external_test_sources", "PPIDM,INVERSE_PPIDM_Silver"],
     {"PPIDM": "Bronze,Silver,Gold", "INVERSE_PPIDM": "Silver"}),
])
def test_ppidm_parse_classes(args, expected, tmp_path):
    proc = preview(tmp_path, *args)
    out = proc.stdout + proc.stderr
    assert proc.returncode == 0, out
    got = dict(re.findall(r"PARSE_(PPIDM|INVERSE_PPIDM) classes: (\S+?)(?=PARSE_|--|\s|$)", out))
    assert got == expected, out


def test_ppidm_class_not_configured_fails_at_startup(tmp_path):
    """Gold is configured, Gold and Silver requested: the partial overlap used to
    pass here and fail only in BUILD_EXTERNAL_TEST at the end of the run."""
    proc = preview(tmp_path, "--ppidm_classes", "Gold",
                   "--external_test_sources", "PPIDM_Gold,PPIDM_Silver,INVERSE_PPIDM_Gold")
    out = proc.stdout + proc.stderr
    assert proc.returncode != 0
    assert "selects PPIDM_Silver, but --ppidm_classes (Gold) does not parse that class" in out, out
