#!/usr/bin/env python3
"""InterPro counters per Pfam family of ddi_split_membership (the interpro_*
columns of domain_metadata.csv), ported from master_ct's
download_interpro_metadata.sh with a release-keyed cache.

InterPro serves only its live release, and its counters (proteins, structures,
AlphaFold models, ...) change with every release. So every response is cached
under <cache>/interpro-<release>/<pfam>.json, the release is read from the API
root before fetching, and a release change while fetching is fatal: a cache
directory never mixes releases.

--release pins a release. When the live release differs, nothing is fetched:
every family must already be in <cache>/interpro-<pin>/, else the run fails.

Per family: an entry the API does not have (HTTP 204/404) gives empty
(NaN) counters, never 0 (the helper filled 0). For an existing entry a
counter key the response lacks is also empty, while a null `literature`
is 0: it counts the references cited in the entry, and null means none.

Counts per family are study-intensity proxies (literature, structures):
they can confound any label they correlate with.

Output: --out TSV (pfam_id + interpro_* columns) and --summary JSON
(embedded in run_manifest.json as interpro_metadata).
"""

import argparse
import csv
import json
import os
import re
import sqlite3
import sys
import tempfile
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from utils_struct import extract_data_from_instance_id  # noqa: E402

API = "https://www.ebi.ac.uk/interpro/api/"
COUNTERS = ["subfamilies", "domain_architectures", "interactions", "matches", "pathways", "proteins",
            "proteomes", "sets", "structures", "taxa", "literature", "alphafold"]
PFAM_RE = re.compile(r"^PF\d{5}$")
RETRY_STATUS = {408, 429, 500, 502, 503, 504}
RELEASE_RE = re.compile(r"^\d+(\.\d+)*$")  # also keeps the cache path safe
RATE_LIMIT_TRIES = 12  # 429 is the API throttling us, not a failure: wait it out (~30 min max)
_pause = {"until": 0.0}  # shared by all worker threads: one 429 pauses every request
_pause_lock = threading.Lock()


def _wait_for_pause():
    left = _pause["until"] - time.monotonic()
    if left > 0:
        time.sleep(left)


class EmptyResponse(Exception):
    """HTTP 200 without a body: transient, retried, never cached as missing."""


def get(url, tries=6, timeout=60):
    """(status, parsed JSON or None); None only for 204/404. Retries transient
    failures with backoff, honouring a numeric Retry-After."""
    attempt = 0
    while True:
        wait = min(60, 2 ** attempt)
        _wait_for_pause()
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers={"Accept": "application/json"}),
                                        timeout=timeout) as r:
                if r.status == 204:
                    return 204, None
                data = json.loads(r.read() or b"null")
                if data is None:
                    raise EmptyResponse(url)
                return r.status, data
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return 404, None
            limit = RATE_LIMIT_TRIES if e.code == 429 else tries
            if e.code not in RETRY_STATUS or attempt >= limit - 1:
                raise
            ra = e.headers.get("Retry-After", "") if e.headers else ""
            if e.code == 429:
                wait = min(300, 5 * 2 ** attempt)
            wait = max(wait, min(300, int(ra))) if ra.isdigit() else wait
            if e.code == 429:  # pause every thread, not just this one
                with _pause_lock:
                    _pause["until"] = max(_pause["until"], time.monotonic() + wait)
        except (urllib.error.URLError, TimeoutError, ConnectionError, EmptyResponse, json.JSONDecodeError):
            if attempt >= tries - 1:
                raise
        time.sleep(wait)
        attempt += 1


def live_release():
    _s, root = get(API)
    dbs = root["databases"]
    rel = str(dbs["interpro"]["version"])
    if not RELEASE_RE.fullmatch(rel):
        raise SystemExit(f"unexpected InterPro release string from the API: {rel!r}")
    return rel, str(dbs["pfam"]["version"])


def parse_entry(status, data):
    """Cache record for one API response."""
    if status in (204, 404):
        return {"status": "missing"}
    meta = data["metadata"]
    counters = dict(meta.get("counters") or {})
    models = counters.pop("structural_models", None) or {}
    if "alphafold" in models:
        counters["alphafold"] = models["alphafold"]
    counters["literature"] = len(meta.get("literature") or {})
    return {"status": "ok", "counters": {k: counters[k] for k in COUNTERS if k in counters}}


def valid(record):
    """A cache record as parse_entry writes it; anything else is not a hit."""
    if not isinstance(record, dict):
        return False
    if record.get("status") == "missing":
        return set(record) == {"status"}
    c = record.get("counters")
    return (record.get("status") == "ok" and isinstance(c, dict) and "literature" in c
            and set(c) <= set(COUNTERS)
            and all(isinstance(v, int) and not isinstance(v, bool) and v >= 0 for v in c.values()))


def read_cached(path):
    try:
        with open(path) as f:
            rec = json.load(f)
    except (OSError, ValueError):
        return None
    return rec if valid(rec) else None


def write_atomic(path, record):
    fd, tmp = tempfile.mkstemp(dir=os.path.dirname(path), prefix=".tmp-", suffix=".json")
    with os.fdopen(fd, "w") as f:
        json.dump(record, f)
    os.replace(tmp, path)


def membership_pfams(db):
    conn = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    try:
        ids = {i for row in conn.execute("SELECT instance_id_a, instance_id_b FROM ddi_split_membership")
               for i in row}
    finally:
        conn.close()
    pfams = sorted({extract_data_from_instance_id(i)[0] for i in ids})
    bad = [p for p in pfams if not PFAM_RE.match(p)]
    if bad:
        raise SystemExit(f"not a Pfam accession: {bad[:5]}")
    return pfams


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--db", required=True)
    p.add_argument("--cache_dir", default="", help="cache root; empty = ./ (work dir, not kept)")
    p.add_argument("--release", default="", help="pinned InterPro release, e.g. 111.0")
    p.add_argument("--out", required=True)
    p.add_argument("--summary", required=True)
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--versions")
    p.add_argument("--process_name", default="FETCH_INTERPRO_METADATA")
    args = p.parse_args(argv)

    if args.release and not RELEASE_RE.fullmatch(args.release):
        raise SystemExit(f"--interpro_release {args.release!r} is not a release number like 111.0")
    pfams = membership_pfams(args.db)
    live, live_pfam = live_release()
    release = args.release or live
    offline = release != live
    if offline and not args.cache_dir:
        raise SystemExit(f"--interpro_release {release} differs from the live release {live}; "
                         "its counters can only come from cache_dir, which is not set")
    cache = os.path.join(args.cache_dir or ".", f"interpro-{release}")

    records, todo = {}, []
    for pf in pfams:
        rec = read_cached(os.path.join(cache, f"{pf}.json"))
        if rec is None:
            todo.append(pf)  # absent, truncated or malformed: refetched when live
        else:
            records[pf] = rec
    if offline and todo:
        raise SystemExit(f"--interpro_release {release} is not live (live: {live}) and {len(todo)} "
                         f"famil(ies) are missing from or invalid in {cache}, e.g. {todo[:5]}")

    # Fetched records stay in memory until the release is confirmed unchanged,
    # so a failed or straddling run never writes into the cache.
    fetched = {}
    with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
        for n, (pf, rec) in enumerate(pool.map(lambda pf: (pf, parse_entry(*get(f"{API}entry/pfam/{pf}/"))),
                                               todo), 1):
            fetched[pf] = rec
            if n % 500 == 0:
                print(f"[interpro] fetched {n}/{len(todo)}", flush=True)
    if todo and live_release()[0] != release:
        raise SystemExit(f"InterPro release changed from {release} while fetching; nothing was cached, re-run")
    bad = [pf for pf, rec in fetched.items() if not valid(rec)]
    if bad:
        raise SystemExit(f"unexpected InterPro counters for {len(bad)} famil(ies), e.g. "
                         f"{[(pf, fetched[pf]) for pf in bad[:3]]}")
    if fetched:
        os.makedirs(cache, exist_ok=True)
        for pf, rec in fetched.items():
            write_atomic(os.path.join(cache, f"{pf}.json"), rec)
    records.update(fetched)

    with open(args.out, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["pfam_id"] + [f"interpro_{c}" for c in COUNTERS])
        for pf in pfams:
            c = records[pf].get("counters", {})
            w.writerow([pf] + [c.get(k, "") for k in COUNTERS])

    missing = [pf for pf in pfams if records[pf]["status"] != "ok"]
    summary = {
        "interpro_release": release,
        "live_interpro_release": live,
        "live_pfam_release": live_pfam,
        "pinned": bool(args.release),
        "families": len(pfams),
        "fetched": len(todo),
        "from_cache": len(pfams) - len(todo),
        "missing_entries": len(missing),
        "missing_examples": missing[:20],
        "missing_counter_values": {k: sum(1 for pf in pfams if records[pf]["status"] == "ok"
                                          and k not in records[pf]["counters"]) for k in COUNTERS},
    }
    with open(args.summary, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"[interpro] {json.dumps(summary)}", flush=True)

    if args.versions:
        with open(args.versions, "w") as f:
            f.write(f'"{args.process_name}":\n')
            f.write(f"    python: {sys.version.split()[0]}\n")
            f.write(f"    interpro: {release}\n")


if __name__ == "__main__":
    main()
