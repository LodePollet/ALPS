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

| CMake option | Default | Effect |
|---|---|---|
| `PIMC_WORM_LATTICES` | all | Lattices to build: `chain`, `square`, `cubic`, `ladder`, `triangular`, `honeycomb` |
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

Each directory in `test/` compares a short run with exact diagonalization:

```sh
cd test/BoseHubbard_GrandCanonical
mkdir run && cd run
pimc_worm_chain ../parameter_files/chain.ini
cd .. && python compare.py chain
```

## Not yet integrated

The code uses its own lattices (`lattice.hpp`) and models (`model.hpp`)
rather than the ALPS lattice and model libraries. Moving to those, and
registering the tests with CTest, are left for later.
