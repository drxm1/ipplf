#!/usr/bin/env bash
# Run E18 (TE) and E19 (TM) on each of the four specified meshes.
set -euo pipefail

if (( $# != 2 )); then
    echo "Usage: bash run_study.sh EXECUTABLE NEW_STUDY_DIRECTORY" >&2
    exit 2
fi

executable=$(realpath -e -- "$1")
[[ -x $executable ]] || { echo "Not executable: $executable" >&2; exit 1; }
inputs=$(cd -- "$(dirname -- "$0")/inputs" && pwd -P)
mkdir -p -- "$(dirname -- "$2")"
mkdir -- "$2"
cd -- "$2"
sha256sum -- "$executable" > executable.sha256
export OMP_NUM_THREADS=4 OMP_PLACES=threads OMP_PROC_BIND=spread

for cells in 24 48 96 192; do
    for experiment in E18 E19; do
        name="$experiment-n$cells"
        echo "Starting $name"
        mpiexec -n 8 --map-by slot:PE=4:HWTCPUS --bind-to hwthread --report-bindings \
            "$executable" --input "$inputs/$name.yaml" > "$name.log" 2>&1
        test -s "$name/run.json"
        cp -- "$inputs/$name.yaml" "$name/input.yaml"
    done
done
