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
    return subprocess.run(
        ["nextflow", "run", REPO, "-profile", "test", "-preview",
         "--outdir", str(tmp_path), *args],
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
