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
under the MIT License of ALPS. Changes made for ALPS:

- MIT/SPDX file headers, and an MIT notice in the start-up banner
  (`copyright.hpp`), which prints the `pimc_worm` entry of the ALPS citation
  catalog (`CITATIONS.yaml`).
- `CMakeLists.txt` rewritten for the ALPS build: it links `ALPS::alps`, whose
  `alps::mcbase`, `alps::params` and accumulators are the ALPSCore ones the
  code was written for, instead of a separate ALPSCore installation.
- `[[maybe_unused]]` on parameters that some lattices or the uniform build
  do not read, so that the code compiles without warnings under
  `-Wall -Wextra`.

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

## Tests

`test/lattice_test.cpp` checks that `alps_lattice` reproduces each of the
hand-written lattices: the same bonds, the same winding increment of every
hop, neighbor tables that agree from both ends of a bond, and the same
relative numbering. It is registered with CTest when `ALPS_BUILD_TESTING=ON`:

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
