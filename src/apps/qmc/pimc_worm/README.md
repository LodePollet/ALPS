# pimc_worm: path integral Monte Carlo with worm-type updates

A continuous-time path integral Monte Carlo code with worm-type updates for
sign-positive Bose-Hubbard and XXZ models, described in

> N. Sadoune and L. Pollet, *Efficient and scalable Path Integral Monte Carlo
> Simulations with worm-type updates for Bose-Hubbard and XXZ models*,
> SciPost Phys. Codebases 9 (2022), [doi:10.21468/SciPostPhysCodeb.9](https://doi.org/10.21468/SciPostPhysCodeb.9).

It is unrelated to the older worm code in [`../worms`](../worms).

## Provenance

Copied from <https://github.com/LodePollet/worm> at commit
`fb738408c5b1c52c0b2a27b614b0285e45b4c554` (2026-08-16). The sources in `src/`
of that repository are the `.cpp` and `.hpp` files here; `parameter_files/`,
`tools/` and `docs/reset_statistics_on_restore.md` are copied unchanged, and
`test/` holds that repository's `test_mpi/` cases (parameter files,
exact-diagonalization references and `compare.py`) without their build scripts.

The upstream repository is licensed under the GPL-3.0. Its copyright holders,
Jonas Greitemann, Nicolas Sadoune and Lode Pollet, relicense the code here
under the MIT License of ALPS. Files added in ALPS in 2026 (`alps_lattice.hpp`,
`alps_input.hpp`, `test/lattice_test.cpp`, `test/input_test.cpp`) are by Lode
Pollet alone. Changes made for ALPS:

- MIT/SPDX file headers, and an MIT notice in the start-up banner
  (`copyright.hpp`), which prints the `pimc_worm` entry of the ALPS citation
  catalog (`CITATIONS.yaml`).
- `CMakeLists.txt` rewritten for the ALPS build: it links `ALPS::alps`, whose
  `alps::mcbase`, `alps::params` and accumulators are the ALPSCore ones the
  code was written for, instead of a separate ALPSCore installation.
- `[[maybe_unused]]` on parameters that some lattices or the uniform build
  do not read, so that the code compiles without warnings under
  `-Wall -Wextra`.
- A detailed-balance fix in `DELETEKINK` (`worm.update.cpp`): its acceptance
  ratio used the coordination number of the site the worm head leaves, where
  the reverse `INSERTKINK` picks the bond with 1/(coordination number of the
  site the head lands on). The two differ only with open boundaries or mixed
  coordination; there soft-core bosons were biased (5 % in the kinetic
  energy of an open four-site chain), periodic lattices were not affected.

## Building

`pimc_worm` is built with the ALPS applications when Eigen 3.3 or newer is
found. The lattice fixes the maximum coordination number, which sizes every
diagram element, so it is a compile-time choice: one executable is built per
lattice, `pimc_worm_<lattice>`, and `pimc_worm_<lattice>_mpi` as well when
`ALPS_ENABLE_MPI=ON`.

`pimc_worm_alps` takes its lattice from the ALPS lattice library instead
(`alps_lattice.hpp`): any lattice whose sites have at most `PIMC_WORM_ZCMAX`
neighbors, selected with `LATTICE`, `L` and, where the lattice has them, `W`
and `H` in place of `Lx`, `Ly`, `Lz` and `pbcx`, `pbcy`, `pbcz`:

```ini
LATTICE = "triangular lattice"
L = 3
W = 3
```

The lattice library is found as in other ALPS programs, or through
`LATTICE_LIBRARY`; in a build tree, set `ALPS_XML_PATH` to
`src/alps/resources`. Sites are numbered as in the ALPS lattice, with the last
dimension running fastest, so per-site output is ordered differently from
`pimc_worm_<lattice>`; the correlation functions are indexed as before. ALPS's
`honeycomb lattice` has the primitive vectors `a1` and `a2` where worm's
`honeycomb` has `a2` and `a2 - a1`, so the two winding-number components
differ, while their lattices agree on an `L x L` torus. See
[docs/alps_lattice_design.md](docs/alps_lattice_design.md).

| CMake option | Default | Effect |
|---|---|---|
| `PIMC_WORM_LATTICES` | all | Lattices to build: `chain`, `square`, `cubic`, `ladder`, `triangular`, `honeycomb`, `alps` |
| `PIMC_WORM_ZCMAX` | `6` | `pimc_worm_alps`: largest number of neighbors of a site |
| `PIMC_WORM_DSTRUC` | `LIST` | Diagram data structure: `AVL`, `LIST` or `LIST_STACK` |
| `PIMC_WORM_UNIFORM` | `ON` | Uniform parameters; `OFF` for site- and bond-dependent ones |
| `PIMC_WORM_CWINDOW` | `OFF` | Canonical window for the Green's function |
| `PIMC_WORM_MATSUBARA` | `OFF` | `G(k, omega_n)` and binned `G(k=0, tau)` measurements |
| `PIMC_WORM_DEBUG` | `OFF` | Check the configuration after every update |

## Running

```sh
pimc_worm_chain parameter_files/BoseHubbard.ini
mpiexec -n 4 pimc_worm_chain_mpi parameter_files/BoseHubbard.ini
```

A run resumes from its checkpoint with `pimc_worm_chain job.clone.h5`;
`--reset-statistics` discards the measurements while keeping the restored
configuration. Results are written to `job.out.h5`.

### From pyalps, and with ALPS parameter names

`pimc_worm` also takes its input the way ALPS's other programs do
(`alps_input.hpp`), so the same runs can be driven from pyalps:

```python
files = pyalps.writeInputH5Files('parm', [{'LATTICE': 'square lattice', 'L': 4,
    'MODEL': 'boson Hubbard', 'T': 0.1, 't': 0.05, 'U': 1.0, 'mu': 0.5, 'Nmax': 2,
    'SWEEPS': 500000, 'runtimelimit': 600}])
pyalps.runApplication('pimc_worm_square', files)
data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix='parm'), 'Stiffness')
```

- An HDF5 file whose `/parameters` group has no ALPSCore attributes, as
  `pyalps.writeInputH5Files` writes it, is an input file, not a checkpoint.
  `parm.taskN.in.h5` writes `parm.taskN.out.h5` and checkpoints to
  `parm.taskN.clone.h5`; values on the command line take precedence.
- The arguments `runApplication` adds are understood: `-T n` is
  `runtimelimit`, `--Tmin`/`--TMax` go to the MPI driver's schedule, and
  `--write-xml` is ignored. With `MPI=n`, run `pimc_worm_<lattice>_mpi`; the
  serial programs refuse `--mpi` rather than run n copies writing the same
  files. XML input (`writeInputFiles`) is not read.
- The ALPS names are accepted in input files, `.ini` files and on the command
  line alike; giving both an ALPS name and the `pimc_worm` name it stands for
  is an error.

| ALPS | `pimc_worm` |
|---|---|
| `MODEL = "boson Hubbard"` | `model = "BoseHubbard"` |
| `T` | `beta` = 1/`T` |
| `t`, or `t0`, `t1`, ... if all equal | `t_hop` |
| `V`, or `V0`, `V1`, ... if all equal | `V_nn` |
| `U`, `Nmax`, `mu` | `U_on`, `nmax`, `mu` |
| `SWEEPS`, `THERMALIZATION`, `SEED` | `sweeps`, `thermalization`, `seed` |
| `LATTICE`, `L`, `W`, `H` | `Lx`, `Ly`, `Lz` of `pimc_worm_<lattice>`, whose `LATTICE` must match (`"square lattice"` for `pimc_worm_square`); `pimc_worm_alps` reads them itself |

XXZ runs keep `pimc_worm`'s names (`model = "XXZ"`, `Jpm`, `Jzz`, `h`). A
`pimc_worm` sweep is one worm cycle, not a sweep over the lattice, and a run
also ends after `runtimelimit` seconds (default 60).

Besides its own observables, `pimc_worm` measures three under the names of
ALPS's older `worm`: `Energy` (= `Total_Energy`), `Density`
(= `Number_of_particles` / number of sites) and `Stiffness`
= Σ_d ⟨W_d²⟩ / (d β), with `W_d` the winding numbers and `d` the number of
periodic directions. `Stiffness` is the superfluid stiffness in two
dimensions; multiply it by L^(2-d) otherwise. Like `Winding_number_squared` it
is measured with `PIMC_WORM_UNIFORM=ON`. It equals `worm`'s on fully periodic
lattices; with open directions `worm` divides by the lattice dimension
instead. `pimc_worm_honeycomb` takes the winding numbers along worm's
primitive vectors, which differ from ALPS's, so its `Stiffness` differs from
`worm`'s and from `pimc_worm_alps`'s on `honeycomb lattice`. A checkpoint written before these
observables were added still resumes; they then hold only the measurements
made after the resume.

[tutorials/mc-10-pimc-worm](../../../../tutorials/mc-10-pimc-worm) computes
the stiffness of the Bose-Hubbard model both ways.

## Tests

`test/lattice_test.cpp` checks that `alps_lattice` reproduces each of the
hand-written lattices: the same bonds, the same winding increment of every
hop, neighbor tables that agree from both ends of a bond, and the same
relative numbering. `test/input_test.cpp` checks the pyalps input and the
ALPS parameter names. `test/open_chain_test.cpp` runs a soft-core
Bose-Hubbard chain with open boundaries for about 20 s and compares its
energies and particle number with exact diagonalization; sites with fewer
neighbors are where the `DELETEKINK` fix matters. All three are registered
with CTest when `ALPS_BUILD_TESTING=ON`:

```sh
ctest -R pimc_worm
```

Each directory in `test/` compares a short run with exact diagonalization:

```sh
cd test/BoseHubbard_GrandCanonical
mkdir run && cd run
pimc_worm_chain ../parameter_files/chain.ini
cd .. && python compare.py chain
```

## Not yet integrated

The models stay in `model.hpp`. `pimc_worm_alps` is the first step towards
ALPS lattices; the remaining steps (one executable for all lattices,
correlation functions on ALPS distance classes, couplings by site and bond
type) are described in [docs/alps_lattice_design.md](docs/alps_lattice_design.md).
The exact-diagonalization comparisons are not yet registered with CTest.
