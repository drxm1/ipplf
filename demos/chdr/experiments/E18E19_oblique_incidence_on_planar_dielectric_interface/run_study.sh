#!/usr/bin/env bash
# Run E18/E19, or their S cases, on each of the four specified meshes.
set -euo pipefail

if (( $# < 2 || $# > 3 )) || [[ ${3:-} != "" && ${3:-} != --s ]]; then
    echo "Usage: bash run_study.sh EXECUTABLE NEW_STUDY_DIRECTORY [--s]" >&2
    exit 2
fi

experiments=(E18 E19)
if [[ ${3:-} == --s ]]; then
    experiments=(E18-S E19-S)
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
    for experiment in "${experiments[@]}"; do
        name="$experiment-n$cells"
        echo "Starting $name"
        mpiexec -n 8 --map-by slot:PE=4:HWTCPUS --bind-to hwthread --report-bindings \
            "$executable" --input "$inputs/$name.yaml" > "$name.log" 2>&1
        test -s "$name/run.json"
        cp -- "$inputs/$name.yaml" "$name/input.yaml"
    done
done
