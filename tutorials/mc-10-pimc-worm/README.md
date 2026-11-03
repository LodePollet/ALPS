# MC-10: Bosons with pimc_worm

The superfluid stiffness of the two-dimensional Bose-Hubbard model,

$$H = -t \sum_{\langle ij\rangle} b_i^\dagger b_j + \frac{U}{2}\sum_i n_i(n_i-1) - \mu \sum_i n_i,$$

on a square lattice with $U = 1$, $\mu = 0.5$ and at most two bosons per site,
computed with `pimc_worm`, the path integral Monte Carlo code with worm-type
updates of N. Sadoune and L. Pollet
([SciPost Phys. Codebases 9 (2022)](https://doi.org/10.21468/SciPostPhysCodeb.9)).
It is the physics of [MC-05](../mc-05-bosons), which uses the older `worm`
code, run in the two ways `pimc_worm` takes its input:

- `stiffness_native.py` writes one `pimc_worm` input file (`.ini`) per point,
  runs `pimc_worm_square` on each, in parallel, and reads the results from the
  `.out.h5` files with h5py. It needs h5py, matplotlib, and `pimc_worm_square`
  in `PATH` or in the directory `PIMC_WORM_BIN`.
- `stiffness_pyalps.py` is [MC-05](../mc-05-bosons)'s `tutorial5a.py` and
  `tutorial5b.py` with `pimc_worm_square` in place of `worm`: the same ALPS
  parameters (`LATTICE`, `L`, `T`, `t`, `U`, `mu`, `Nmax`, `SWEEPS`), written
  with `pyalps.writeInputH5Files`, and the same `Stiffness` observable. It runs
  the points one after the other.

```sh
python stiffness_native.py a     # or b, or ab (default)
python stiffness_pyalps.py a
```

## a: stiffness against hopping

$\rho_s = \sum_d \langle W_d^2\rangle / (d\beta)$ from the winding numbers
$W_d$ in the $d = 2$ periodic directions, on a $4\times4$ lattice at $T = 0.1$. It rises from zero in the Mott
insulator at small $t$ to the superfluid at large $t$. The native script takes
about a minute on several cores, the pyalps one about three; with the same
seeds the two give identical numbers. `worm` and `pimc_worm` agree within
their error bars:

| $t$ | `worm` | `pimc_worm` |
|---|---|---|
| 0.02 | 0.00019(3) | 0.00018(4) |
| 0.04 | 0.00491(16) | 0.00481(16) |
| 0.06 | 0.0358(3) | 0.0348(4) |
| 0.08 | 0.0853(5) | 0.0860(7) |
| 0.10 | 0.1266(6) | 0.1256(9) |

## b: finite-size scaling

$\rho_s L$ against $t$ for $L = 4, 6, 8$ at $T = 0.05$. At the superfluid-Mott
transition, which has dynamical exponent $z = 1$, $\rho_s L^{d+z-2} = \rho_s L$
is independent of $L$, so the curves cross there, between $t = 0.055$ and
$0.0575$. The native script takes about four minutes; the pyalps one runs the
21 points one after the other.

## Parameters

`pimc_worm` has its own parameter names and accepts the ALPS ones as well:

| ALPS | `pimc_worm` |
|---|---|
| `LATTICE`, `L`, `W`, `H` | `Lx`, `Ly`, `Lz`, `pbcx`, `pbcy`, `pbcz` (or `pimc_worm_alps` for any ALPS lattice) |
| `MODEL = "boson Hubbard"` | `model = "BoseHubbard"` |
| `T` | `beta` $= 1/T$ |
| `t` (or `t0`, `t1`, ... if all equal) | `t_hop` |
| `U`, `V`, `mu`, `Nmax` | `U_on`, `V_nn`, `mu`, `nmax` |
| `SWEEPS`, `THERMALIZATION`, `SEED` | `sweeps`, `thermalization`, `seed` |

A `pimc_worm` sweep is one worm cycle, not a sweep over the lattice; 500000 of
them give about the error bars of `worm`'s 500000 sweeps here. A run also ends
after `runtimelimit` seconds (default 60), so both scripts set it.
See [the pimc_worm README](../../src/apps/qmc/pimc_worm/README.md) for the other
observables and the build options.
