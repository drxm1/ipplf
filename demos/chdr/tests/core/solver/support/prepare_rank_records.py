"""Create fresh fixture directories, retaining only the immediately preceding run.

Usage: prepare_rank_records.py ROOT. The script moves ROOT/serial and ROOT/mpi2
of an earlier run to ROOT/previous-ranks and creates empty ROOT/serial and
ROOT/mpi2. CTest runs it as chdr.solver.ranks.prepare and
chdr.solver.medium.ranks.prepare, the fixture setup before the record producers.
It checks nothing; the empty directories make a missing producer output fail the
comparison in compare_ranks.py.
"""
from pathlib import Path
import sys
import shutil

root = Path(sys.argv[1])
root.mkdir(parents=True, exist_ok=True)
existing = [root / name for name in ("serial", "mpi2") if (root / name).exists()]
if existing:
    archive = root / "previous-ranks"
    if archive.exists():
        shutil.rmtree(archive)
    archive.mkdir()
    for path in existing:
        path.rename(archive / path.name)
for name in ("serial", "mpi2"):
    (root / name).mkdir()
print(f"Fresh rank records: {root / 'serial'}, {root / 'mpi2'}")
