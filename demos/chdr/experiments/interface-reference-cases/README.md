# Planar-interface reference tools

These tools run explicit `planar-interface` inputs, read their recorded results,
and evaluate continuum and aligned-lattice references through a small C++ program.
[references.org](references.org) derives the coefficients, field conventions,
phase transfer and applicability limits. Calculating a reference or completing a
simulation does not by itself establish agreement with the physical problem.

The C++ utility uses only standard C++ and the independent functions in
[`InterfaceOracles.h`](../../tests/apps/planar-interface/InterfaceOracles.h).
It does not initialize MPI or an IPPL simulation. Python performs orchestration,
record validation and analysis; it does not reimplement the Fresnel formulas.

## Build and requirements

In an IPPL build configured with the ChDR applications, build:

```sh
cmake --build /path/to/build --target planar-interface interface-reference
```

The Python readers need Python 3.10 or newer and its standard library.
`scan.py` additionally needs PyYAML and a POSIX system. An MPI launcher is needed only when requesting
more than one rank. The focused CTest check is `chdr.interface.reference-tools`.

## Run explicit inputs

For example, an aligned TE row can be requested with:

```yaml
schema_version: 1
case: plane
units: {system: normalized}
parameters:
  polarization: TE
  incidence: vacuum
  epsilon: 2.13
  angle_degrees: 30
  cells_per_wavelength: 20
  cell_shape_y: 1
  cell_shape_z: 1
```

Numerical resolution, modes, pulse placement and gates are calculated by the C++
application. The runner changes only `output/directory`; it preserves the spelling
of every physical parameter, including exponent notation. Inputs must be explicit
YAML mappings without anchors or aliases. An existing `output` mapping must contain
a scalar `directory`; otherwise omit `output` and the runner adds it.

```sh
python3 demos/chdr/experiments/interface-reference-cases/scan.py \
  --binary /path/to/immutable/planar-interface \
  --output-root /path/to/build/interface-reference-cases \
  --ranks 1 --threads 1 --jobs 2 --max-cores 2 \
  normal.yaml oblique.yaml
```

Use a copied immutable executable if the source or build is changing during a
campaign. Its SHA-256 is checked before and after each run. The runner records the
binary's compiled identity from `run.json` separately from the analysis checkout's
identity at launch; those two revisions need not be identical.

`jobs * ranks * threads` must fit the explicit CPU budget, which cannot exceed the
current process's local CPU affinity. Each job receives `OMP_NUM_THREADS` and
`OMP_PROC_BIND=false`. Multi-rank runs use `--mpi` (default `mpiexec`) and `-n`;
extra arguments can be supplied as repeated `--mpi-arg=VALUE` options. Use
`--timeout SECONDS` for a per-attempt time limit. Each subprocess is invoked with an
argument list, without a command interpreter. Each attempt has its own process
session. On timeout the runner sends `SIGTERM` to the local process group, allows
two seconds for shutdown, then sends `SIGKILL` if needed. Cleanup waits are bounded,
and the launcher is reaped. The failure record reports whether the process group
still exists after cleanup; killed orphan zombies may remain until the system
reaper collects them. The retained stdout/stderr include output written before
timeout. Process-group signals cover local descendants that remain in that group;
remote MPI workers or workers that start a separate session require their
launcher's job cleanup.

Pressing `Ctrl-C` cancels the scan and returns exit status 130. The runner blocks
new process launches, terminates all registered local attempt groups, then joins
its workers. Running and queued attempts are recorded as `cancelled`; rows that
already completed or failed retain those statuses. The final manifest includes
the cancellation reason and group-cleanup results, alongside the retained logs
and inputs. A further interrupt is ignored while this bounded cleanup and record
write finish. Cancellation does not require a per-attempt timeout.

Each row receives a fresh directory:

```text
OUTPUT_ROOT/
  scan-SCAN_ID.json
  ROW_ID/SCAN_ID/
    input.original.yaml
    input.yaml
    attempt.json
    stdout.txt
    stderr.txt
    output/
      run.json
      spectra.csv
      ...other selected application outputs...
```

Original inputs are copied before any simulation starts. Repeating a scan creates
new attempts; existing results are not replaced. Preparation, execution, timeout
and output-reading failures remain in the manifest with their reasons. CLI exit
status is one if any attempt failed. A completed entry means that a completed run
record and a consistent spectrum were read. Its `physics_validation` remains
`not_performed`; agreement, boundary isolation and angular errors require separate
analysis.

## Read and evaluate one completed run

```sh
python3 demos/chdr/experiments/interface-reference-cases/read_run.py RUN_OUTPUT
python3 demos/chdr/experiments/interface-reference-cases/references.py RUN_OUTPUT \
  --binary /path/to/interface-reference --output NEW_REFERENCE_DIRECTORY
```

The reader supports aligned-row records with plane schema zero or one and tilted
records with plane schema two. An absent version denotes the first aligned record.
Each version requires explicit normalized effective metadata and its declared
geometry. Unsupported geometry/schema combinations are rejected rather than
interpreted from directory names. Version one validates the declared transform,
raw units and separate probe/coefficient phase origins.
SI frequency columns are converted with the recorded output factor and checked
against `omega_code`. Raw transforms remain in normalized solver units.

An invalid frequency bin is retained with its reason. Undefined coefficients are
represented as `None` by the Python reader; invalid raw non-finite values remain
available for diagnosis. Use `valid_m` when selecting data for comparisons. Do not
turn an invalid-bin mask into a zero reflected or transmitted signal.

Schema two requires `geometry: periodic_tilted_plane`, the normal/tangent basis,
material-interface anchor, periodic geometry, probe layout and the recorded qualification
quantities. `spectra.csv` stores the fitted incident, reflected and transmitted
interface amplitudes, three signed normal roots and their group velocities, fit
residuals, and the derived coefficients. These columns have different meanings
from the raw point transforms of schema one.
The requested material-interface position is a signed normal offset; the recorded anchor must
lie on that plane. Both x and z windings are retained when lifting probe points.

The required `probe-spectra.csv` retains every point's wrapped component position,
integer x/z box windings, lifted normal/tangent coordinates, and raw incident and
total transforms. The reader checks coordinate consistency, point counts and
coverage of both probe sides at every frequency. A declared-valid bin must satisfy
the recorded fit-residual limit and positive Floquet exclusion margins. These
checks verify the record's consistency; they do not establish boundary isolation
or bound the remaining evanescent contribution.

The reference directory is new and contains `frequencies.txt`, `reference.csv`,
`stderr.txt` and `reference.json`. The JSON record preserves the exact command,
input/output identities, reference executable hash and compiled revision, analysis
source identity, and reference-document hash. It also records discrepancies between
published coefficients and a reconstruction from raw spectra when that audit is
applicable. For schema two, the audit repeats the recorded spatial fits with the
three supplied directional roots, checks all interface amplitudes and normalized fit
residuals, and then checks their ratios. It does not replace the reflected root
by the negative incident root. The raw point file is included in the provenance
hashes. These discrepancies test extraction consistency, not physical accuracy.

## Direct reference-program interface

Write one normalized angular frequency per line in `frequencies.txt`, then run:

```sh
interface-reference --polarization TM \
  --epsilon-incident 1 --epsilon-transmitted 2.13 \
  --spacing 1,1,1 --dt 0.5196152422706632 --kz 0 \
  --frequencies frequencies.txt > reference.csv
```

`--frequencies -` reads stdin. `--identity` reports the reference schema and compiled
source identity. Each output row contains `omega_code`, `kz`, continuum complex roots
and coefficients, continuum `R` and `T`, the continuum incidence angle, and aligned
lattice roots/coefficients. Continuum and lattice validity/reasons are separate.
A bin can have propagating lattice roots without a propagating continuum incident
root; no physical flux normalization is emitted for the latter. Complex transmitted
roots can still have valid reference coefficients, even when a measurement method
does not support that branch. See the derivations for the exact meanings.

`--continuum-only` disables the aligned-lattice result. For a rotated continuum
reference, `--kz` denotes the signed tangential wavenumber in the interface's local
normal/tangent basis; it must not be copied blindly from a laboratory wavevector.
The aligned lattice formulas are not an exact oracle for a rotated material map.
The Python wrapper also rejects a requested material interface between cell faces for this
unshifted lattice formula; use `--continuum-only` for the requested-interface shifted
comparison. It checks the recorded Standard/staggered/O1VS configuration before
requesting the aligned lattice reference.
For a schema-two run, call the Python wrapper with `--continuum-only`. It supplies
the recorded `kappa` to the C++ local-tangent input, checks the central reference
angle against the actual recorded angle, and rejects an aligned-lattice request.
The directional roots used in the raw-fit audit are recorded propagation roots;
they are not an exact scattering oracle for the staircased interface.

For a dense continuum curve, write one incidence angle in degrees per line and use:

```sh
interface-reference --polarization TM \
  --epsilon-incident 1 --epsilon-transmitted 2.13 \
  --spacing 1,1,1 --dt 0.5196152422706632 \
  --angles angles_degrees.txt --omega 0.3 --continuum-only > angle-reference.csv
```

Angles lie in `[0, 90)` and are relative to the material normal. The C++ program
calculates the tangential wavenumber from each angle and the supplied normalized
frequency, as derived in [references.org](references.org). This mode requires
`--continuum-only` and excludes `--kz` and `--frequencies`. It samples continuum
references; it does not resolve a new simulation row or define a lattice angle.

## Python interfaces for case analyses

- `read_run.load_run(path)` returns `PlaneRun`. Its `plane`, `spacing`, `dt` and
  `frequencies` properties expose normalized effective values. `samples_m` retains
  bin order and `central_bin_m` identifies the central bin.
- Each `SpectrumSample` has `omega_m`, `reflection_m`, `transmission_m`, `valid_m`,
  `reason_m`, and optional raw `incident_m`, `total_m`, `transmitted_m` transforms.
  Schema two instead provides `face_amplitudes_m`, `normal_wavenumbers_m`,
  `normal_group_velocities_m` and `fit_residuals_m`, each ordered incident,
  reflected, transmitted. Its raw points are retained in `PlaneRun.probes_m`.
- `references.evaluate(run, binary, new_directory, continuum_only=False)` returns
  a `ReferenceTable` with `rows_m` and its provenance in `record_m`.
- `references.angle_curve(run, binary, angles_degrees, new_directory)` returns a
  continuum `ReferenceTable` at the run's recorded central frequency, polarization
  and media. It archives `angles_degrees.txt`, the C++ invocation, table and source
  identities. Python supplies the angle vector; C++ computes the wavevectors and
  reference coefficients. For example, `[index / 2 for index in range(161)]` gives
  the 0–80 degree range at half-degree intervals.
- `references.coefficient(row, "continuum" or "lattice", "r" or "t" or "k1" or
  "k2")` returns the corresponding complex reference or `None` when unavailable.
  Power columns are `continuum_R` and `continuum_T`; the angle column is
  `continuum_angle_degrees`. For measured transmission power use the explicitly
  calculated `continuum_transmission_weight` times the squared magnetic/electric
  amplitude magnitude in the recorded basis; do not infer that weight by dividing
  two potentially small reference quantities.
- `references.extraction_residuals(run, rows)` audits the probe-to-interface phase transfer
  using the C++ lattice roots for aligned data, or the recorded directional roots
  and point fits for tilted data. The selected measured bins must be valid.

Keep the run input, resolved record, spectra, reference invocation and output,
comparison data and relevant control runs together when retaining a result. Store
the selected frame/series data needed to reproduce reported field profiles or
movies, not only their rendered pictures. Artifact retention and storage policy
belong to the surrounding experiment workflow; these programs have no external
research-directory dependency.
