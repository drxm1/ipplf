"""Create fresh fixture directories, retaining only the immediately preceding run."""
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
