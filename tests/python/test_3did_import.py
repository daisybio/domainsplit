#!/usr/bin/env python3
"""bin/3did_import.py keeps only DDI1 and Domain, deciding per statement.

The skipped table's INSERT spans several lines whose continuation rows look like
nothing in particular, and a quoted value holds ``;`` plus a newline -- a
line-prefix filter would mis-handle both.

Run directly (`python3 tests/python/test_3did_import.py`) or via pytest.
"""

import os
import sqlite3
import subprocess
import sys
import tempfile

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRIPT = os.path.join(REPO, "bin", "3did_import.py")

DUMP = """PRAGMA synchronous = OFF;
BEGIN TRANSACTION;
CREATE TABLE `Domain` (
  `domain_id` integer NOT NULL
,  `pfam_id` varchar(20) NOT NULL
);
CREATE TABLE `residues` (
  `id` integer NOT NULL
,  `note` text
);
CREATE TABLE `DDI1` (
  `ddi_id` integer NOT NULL
,  `domain1` integer
,  `domain2` integer
);
INSERT INTO `Domain` VALUES (1,'PF00001'),(2,'PF00002');
INSERT INTO `residues` VALUES
(1,'a;
b'),
(2,'it''s;');
INSERT INTO `DDI1` VALUES (10,1,2);
INSERT INTO `DDI1` VALUES (10,1,2);
CREATE INDEX "idx_residues_id" ON "residues" (`id`);
CREATE TRIGGER [residues_touch] AFTER UPDATE ON `residues` FOR EACH ROW BEGIN UPDATE residues SET id = id; END;
CREATE INDEX "idx_DDI1_d1" ON "DDI1" (`domain1`);
END TRANSACTION;
"""


def test_only_ddi1_and_domain_are_imported():
    with tempfile.TemporaryDirectory() as tmp:
        dump, db = os.path.join(tmp, "dump.sql"), os.path.join(tmp, "3did.sqlite3")
        with open(dump, "w") as fh:
            fh.write(DUMP)
        subprocess.run([sys.executable, SCRIPT, dump, db], check=True, capture_output=True)
        conn = sqlite3.connect(db)
        tables = {r[0] for r in conn.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        assert tables == {"Domain", "DDI1"}
        assert conn.execute("SELECT * FROM Domain ORDER BY domain_id").fetchall() == [(1, "PF00001"), (2, "PF00002")]
        assert conn.execute("SELECT * FROM DDI1").fetchall() == [(10, 1, 2), (10, 1, 2)]
        assert conn.execute("SELECT name FROM sqlite_master WHERE type='index'").fetchall() == [("idx_DDI1_d1",)]
        conn.close()


if __name__ == "__main__":
    test_only_ddi1_and_domain_are_imported()
    print("OK")
