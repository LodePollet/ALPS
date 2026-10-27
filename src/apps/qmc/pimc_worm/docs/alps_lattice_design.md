# pimc_worm on ALPS lattices: design notes

Status: step 1 of section 6 is implemented (`alps_lattice.hpp`,
`pimc_worm_alps`, `test/lattice_test.cpp`); the rest is a proposal.
Feedback from the ALPS developers:
the models (Bose-Hubbard, XXZ) stay hard-coded in `model.hpp`; the lattice
should come from the ALPS lattice library (`alps::graph_helper<>`).

## 1. What worm needs from a lattice today

Everything goes through `lattice.hpp` and is copied into flat tables in the
`worm` constructor (`worm.cpp`, around line 123):

| Worm quantity | Today | Used for |
|---|---|---|
| `Nsites`, `get_Nbonds()` | `lattice<DIM,N_BASIS>` | sizes, inhomogeneous parameter vectors |
| `nb[s][d]`, `-1` if absent | per-site map, direction `d < zcmax` | every update, `mAssoc` |
| `opposite_direction[s][d]` | `BondDirections` map | re-linking `mAssoc` on the neighbor |
| `bond_index[s][d]` | `BondIndices` map | non-uniform bond couplings |
| `zcmax` (compile time) | `LATTICE::zcmax` | size of `Element::mAssoc[ZC]` in every vertex |
| `winding_element[d]` | `directions()`, integer crystal vectors | `mWinding += w[d] * dn` in `INSERTKINK`/`DELETEKINK` |
| `get_Ls(i)` | `Lx, Ly, Lz` | `W^2 = mWinding^2 / L^2`, Matsubara k grid |
| `relNumbering(s1,s2)` | translation index, `N_BASIS == 1` only | density matrix, density-density, G(k, omega) |
| `pbcx/y/z`, `n_basis` | parameters | switch off correlation functions unless fully periodic Bravais |

## 2. Mapping onto `alps::graph_helper<>`

Construction: `alps::graph_helper<> g(alps::make_deprecated_parameters(parameters));`
(`<alps/ngs/make_deprecated_parameters.hpp>` turns the ALPSCore `alps::params`
worm uses into the `alps::Parameters` the lattice library reads). The lattice is
then selected with the usual `LATTICE`, `L`, `W`, `H`, `LATTICE_LIBRARY`.

### 2.1 Neighbor tables, both directions

Build all three tables in one pass over the bonds, in bond-index order:

```
for each bond b = (s, t):
    ds = deg[s]++ ; dt = deg[t]++
    nb[s][ds] = t ; nb[t][dt] = s
    opposite_direction[s][ds] = dt ; opposite_direction[t][dt] = ds
    bond_index[s][ds] = bond_index[t][dt] = b
    winding[s][ds] = +bond_vector_relative(b) ; winding[t][dt] = -bond_vector_relative(b)
```

- The direction index is just the slot in the site's neighbor list; it no longer
  means "+x". Nothing in the updates relies on a geometric meaning except the
  winding, which gets its own per-(site, slot) table below.
- Open boundaries, depleted or irregular lattices give a smaller `deg[s]`; the
  remaining slots stay `-1`, as worm's pseudo-site already handles.
- The order is deterministic (bond index), so `mAssoc` slots in a checkpoint stay
  meaningful.
- Reject what worm cannot represent: two bonds between the same pair of sites
  (ALPS produces them for periodic extents of 1 or 2, where worm already
  throws) and self-loops. `linkdir` lookups (`nb[isite][zci] == adj_site`) assume
  uniqueness.

### 2.2 `zcmax` and the per-vertex `mAssoc[ZC]`

`ZC` must stay a compile-time constant: a runtime-sized `mAssoc` in every
`Element` would cost an allocation or indirection per vertex. Options:

- **A (recommended):** template `Element`/`worm` on `ZC`, instantiate a few
  values (e.g. 2, 3, 4, 6, 8, 12) and dispatch at run time on the lattice's
  maximum degree. One executable `pimc_worm` replaces the six per-lattice ones;
  more compile time, same run-time cost as today.
- B: one CMake option `PIMC_WORM_ZCMAX`, with a run-time check that the
  lattice's maximum degree fits.

### 2.3 Winding numbers

ALPS stores for every bond `bond_vector_relative`: the bond vector in lattice
(crystal) coordinates divided by the extent in that direction
(`latticegraph.h`). Summed over a closed world line it is the winding number
itself.

Step 1 uses instead the bond's displacement in whole cells: the difference of
the cell offsets of its ends plus the extent for each periodic boundary it
crosses (`boundary_crossing_t`). That is exactly worm's own gauge (zero for a
bond inside a cell), so the winding increments match worm's bond for bond and
the division by `L^2` stays. The points below apply if `bond_vector_relative`
is adopted later:

- `mWinding` becomes a run-time vector of length `g.dimension()`
  (`std::vector<double>`); Eigen is no longer needed.
- `W_i^2 = mWinding_i^2`, without the division by `L_i^2`.
- The gauge differs from worm's for honeycomb (ALPS uses the basis-site
  positions, worm put the whole displacement on the inter-cell bonds). Both
  give the same integers for closed loops, and `measure()` is only called in
  the diagonal sector (`update()` loops until `worm_diag`), so results agree.
  Checkpoints written with the old gauge must not be resumed.
- Winding numbers and the superfluid stiffness need periodic boundaries along
  the Bravais vectors, not translation invariance: hopping and potential terms
  may break it (disorder, superlattices, traps) and the estimator stays valid.
  Gate it on the periodic directions only, never on homogeneity or `UNISYS`.
- Open directions have no wrapping bonds, so their winding is zero; measure
  only the periodic components and normalize the stiffness by their number.
- A plain `GRAPH` without a lattice has no bond vectors: omit winding and
  stiffness there, measure everything else, and say so in one start-up line
  that names the reason ("no periodic lattice directions", not "no
  translation symmetry"). Today `alps_lattice` assumes cells and would fail
  on a `GRAPH`; this must become a clean path with a test.
- The winding estimator is currently compiled only with `UNISYS`, but it does
  not need uniform couplings; it can be measured always.

### 2.4 Distances and correlation functions

`relNumbering` becomes `g.distance(s1, s2)`, with `g.num_distances()` and
`g.distance_labels()` written to the output. This lifts the present limits
(single-site unit cell, all directions periodic): ALPS distance classes include
the basis sites, and on inhomogeneous lattices fall back to all site pairs.

Two estimators change in substance:

- `DensDens_CorrFun` now takes one snapshot relative to site 0. With distance
  classes it can average over all pairs in a class, which is better statistics
  at O(N^2) cost per measurement (every `Nmeasure2` steps). Keep the
  reference-site version as an option if that cost matters.
- `G(k, omega_n)` uses an `Lx * Ly` FFT. With ALPS, take `g.momenta()` and
  coordinates; keep it to Bravais lattices first.

`compare.py` in `test/` indexes these vectors by worm's `relNumbering`; it must
map through the distance labels instead.

## 3. `UNIFORM` and the coupling back into the model

`UNISYS` currently decides three separate things at compile time:

1. **Where couplings come from:** scalars from the parameters, or per-site and
   per-bond vectors filled by `user_set_params()` from values hard-coded in
   `worm.cpp` (a placeholder).
2. **Speed:** with `UNISYS` the weights are precomputed tables
   (`bond_weight_diag[n][m]`), without it every weight is a virtual call.
3. **Which measurements exist:** the translation-invariant ones only with `UNISYS`.

Proposal:

- **Couplings indexed by ALPS site and bond types.** The model stays hard-coded,
  but each coupling has a value per type: `t_hop`, `V_nn`, `Jxy`, `Jzz` per bond
  type, `U_on`, `mu`, `h` per site type, read as `t_hop0`, `t_hop1`, ... with the
  untyped name as default (the ALPS `t#`/`J#` convention). `g.site_type(s)` and
  `g.bond_type(b)` are the coupling from lattice to model.
- **Precomputed tables per type:** `site_weight_diag[stype][n]`,
  `bond_weight_diag[btype][n][m]`, and the off-diagonal ones. Same cost as
  `UNISYS` now (one extra index), so the compile-time flag can go. A uniform
  system is the one-type case; ladders (ALPS's `ladder` already gives legs and
  rungs different bond types), anisotropic and superlattice models come for free.
- **Disorder:** ALPS `DISORDER` gives every site/bond its own inhomogeneous
  type, but no values. Random on-site potentials need a value source: a
  per-site/per-bond file (HDF5 or text) named by a parameter, replacing the
  hard-coded `user_set_params()` block. Tables per site are `Nsites * (nmax+1)`.
- **Measurements gated at run time:** estimators indexed by translation
  (`relNumbering`, `G(k, omega_n)`) are switched off when the lattice is
  inhomogeneous or site/bond values come from a file; energies and densities
  stay on, and winding and stiffness stay on whenever the lattice has periodic
  directions (section 2.3).
- **Sign check:** the weights are used as probabilities. `XXZ` returns
  `0.5 * Jpm * ...` and Bose-Hubbard `t * ...` unchanged, so a negative
  off-diagonal coupling silently produces negative weights. On ALPS lattices
  (frustrated ones included) check at start-up: off-diagonal couplings must be
  sign-free, i.e. ferromagnetic `Jxy`/positive `t`, or the lattice bipartite
  (`g.is_bipartite()`, `g.parity(s)`) so the sign can be gauged away.

## 4. Parameters and build

- `Lx, Ly, Lz, pbcx, pbcy, pbcz` give way to `LATTICE`, `L`, `W`, `H` and the
  boundary conditions of the lattice description. Worm's names map to ALPS
  lattices as: `chain` -> `chain lattice`, `square` -> `square lattice`,
  `cubic` -> `simple cubic lattice`, `ladder` -> `ladder`,
  `triangular` -> `triangular lattice`, `honeycomb` -> `honeycomb lattice`.
- `PIMC_WORM_LATTICES` and `PIMC_WORM_UNIFORM` disappear; Eigen is no longer
  required.
- The test parameter files are rewritten accordingly.

## 5. Validation

1. A unit test that builds each of worm's six lattices and the ALPS
   counterpart and checks that the edge sets agree up to relabeling, and that
   winding increments agree on closed loops. This settles the honeycomb and
   triangular geometries before any physics runs.
2. Bit-for-bit comparison with the current code is lost once the neighbor
   order changes (the updates draw directions at random). Use the ED
   benchmarks in `test/` for all four test cases and all six lattices.
3. New cases the port enables: open boundaries, a two-type ladder, and a
   disordered chain against ED.

## 5a. What step 1 found

- ALPS numbers cells with the last dimension running fastest; worm runs x
  fastest. The adapter keeps ALPS's site numbering and computes
  `relNumbering` from cell offsets, so the correlation functions are indexed
  as before, while per-site output (`Density_Distribution`) comes in ALPS order.
- ALPS's `honeycomb lattice` uses the primitive vectors `a1`, `a2`; worm's
  `honeycomb` uses `a2` and `a2 - a1`. On an `L x L` torus they are the same
  lattice (a worm cell `(x, y)` is the ALPS cell `(-y, x + y)`), but the two
  winding components are taken along different vectors. For `L != W` the tori
  differ.
- The other five lattices agree with ALPS's `chain lattice`, `square lattice`,
  `simple cubic lattice`, `ladder` and `triangular lattice`, including the open
  variants of chain and square.
- `ZC` is a CMake option for `pimc_worm_alps` (`PIMC_WORM_ZCMAX`, default 6);
  a lattice with more neighbors is refused at start-up.

## 5b. Succeeding `worms`

The older worm code in `src/apps/qmc/worms` has no maintainer: its author has
left the collaboration, and only `pimc_worm` is to be maintained. The plan,
stated in the discussion of ALPSim/ALPS#254 (winding numbers on `GRAPH`
lattices, issue #43), is to deprecate `worms` once `pimc_worm` takes its
inputs, accepting small fixes to `worms` until then. `pimc_worm` can in
principle simulate every model `worms` reasonably can; what it must provide
before the switch:

- Any ALPS lattice, including a plain `GRAPH`, with the winding behaviour of
  section 2.3 (and the start-up note), tested.
- Couplings by bond and site type with the ALPS names (`t#`, `V#`, `U`, `mu`,
  ...), where an unset type falls back to the untyped value rather than to
  zero. ALPS's `boson Hubbard` defaults unset `t#` to zero, which on
  `honeycomb lattice` (three bond types) silently gives `t1 = t2 = 0`.
- Inhomogeneous and disordered couplings (section 3).
- The measurements `worms` reports (energies, densities, correlations,
  Green's function, winding and stiffness), with names mapped in the docs.
- ED benchmarks at a temperature that matches the reference: the honeycomb
  2x2 reference is a ground-state value and needs `beta >= 60` (at
  `beta = 20` both codes sit 3-5 sigma off in the same direction, N = 7.016).

## 6. Suggested order

1. Lattice adapter: build `nb`, `opposite_direction`, `bond_index`, winding
   table and distances from `graph_helper`, behind worm's current interface;
   keep `UNISYS` and one `ZC` per build. Validate (5.1, 5.2).
2. `ZC` dispatch, one executable.
3. Correlation functions on ALPS distance classes; update `compare.py`.
4. Type-indexed coupling tables; remove `UNISYS`.
5. Disorder input from a file.

## 7. Open questions

- Parameter names for type-dependent couplings (`t_hop0`/`t_hop1`, or ALPS's
  `t0`/`t1`, `J0`/`J1`)?
- Keep the old `Lx`/`pbcx` parameters as aliases for existing input files?
- Format of the disorder file?
- `ZC` dispatch (A) or a single compile option (B)?
- Is the single-reference-site `DensDens_CorrFun` intentional for cost, or
  should it average over pairs?
