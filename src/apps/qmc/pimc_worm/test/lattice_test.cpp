/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2021-2026 by Jonas Greitemann, Nicolas Sadoune and Lode Pollet
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Checks that alps_lattice, built from the ALPS lattice library, reproduces
// the hand-written lattices of lattice.hpp: the same sites and bonds, the same
// winding increment for every hop, consistent neighbor tables in both
// directions, and (for Bravais lattices) the same relative numbering.

#include "../lattice.hpp"
#include "../alps_lattice.hpp"

#include <gtest/gtest.h>

#include <functional>
#include <set>
#include <string>

namespace {

alps_lattice make_alps(std::string const& name, int L, int W = 0, int H = 0) {
    alps::Parameters p;
    p["LATTICE"] = name;
    p["L"] = std::to_string(L);
    if (W) p["W"] = std::to_string(W);
    if (H) p["H"] = std::to_string(H);
    return alps_lattice(p);
}

using Map = std::function<int(int)>;                       // worm site -> ALPS site
using Transform = std::function<alps_lattice::Direction(std::vector<double> const&)>;
using Cell = std::array<int, 3>;

// Maps worm's site numbering, cell (x, y, z) -> site (x + Lx*(y + Ly*z)) *
// n_basis + b, onto the ALPS site of cell cell_map(x, y, z) and basis b. ALPS
// numbers its cells in another order, so the map goes through cell offsets.
Map site_map(alps_lattice const& a, std::array<size_t, 3> const& L, std::function<Cell(Cell)> const& cell_map) {
    std::map<std::pair<Cell, int>, int> site;
    size_t const nb = a.get_nbasis();
    for (int s = 0; s < (int)a.get_Nsites(); ++s) site[{a.cell_offset(s), int(s % nb)}] = s;
    return [=](int s) {
        int const cell = s / nb;
        Cell c{int(cell % L[0]), int(cell / L[0] % L[1]), int(cell / (L[0] * L[1]))};
        Cell o = cell_map(c);
        for (size_t d = 0; d < 3; ++d) o[d] = (o[d] % int(L[d]) + int(L[d])) % int(L[d]);
        return site.at({o, int(s % nb)});
    };
}

Map identity_cells(alps_lattice const& a, std::array<size_t, 3> const& L) {
    return site_map(a, L, [](Cell c) { return c; });
}

alps_lattice::Direction pad(std::vector<double> const& w) {
    alps_lattice::Direction v = alps_lattice::Direction::Zero();
    for (size_t i = 0; i < w.size(); ++i) v[i] = w[i];
    return v;
}

// Checks the tables of a on their own: every neighbor relation is mutual, the
// opposite direction leads back, the bond index and the winding agree from
// both ends.
void check_consistent(alps_lattice const& a) {
    std::set<size_t> bonds;
    for (int s = 0; s < (int)a.get_Nsites(); ++s) {
        EXPECT_LE(a.get_zcoord(s), alps_lattice::zcmax);
        for (size_t d = 0; d < alps_lattice::zcmax; ++d) {
            if (!a.has_neighbor(s, d)) continue;
            int const t = a.nb(s, d);
            size_t const back = a.get_opposite_direction(s, t);
            ASSERT_TRUE(a.has_neighbor(t, back));
            EXPECT_EQ(a.nb(t, back), s);
            EXPECT_EQ(a.get_direction(s, t), d);
            EXPECT_EQ(a.get_bond_index(s, d), a.get_bond_index(t, back));
            EXPECT_EQ(a.winding_vector(s, d), -a.winding_vector(t, back));
            bonds.insert(a.get_bond_index(s, d));
        }
    }
    EXPECT_EQ(bonds.size(), a.get_Nbonds());
    if (!bonds.empty()) EXPECT_EQ(*bonds.rbegin() + 1, a.get_Nbonds());
}

template <class WORM>
void compare(WORM& w, alps_lattice const& a, Map const& map, Transform const& transform, bool same_numbering) {
    ASSERT_EQ(w.get_Nsites(), a.get_Nsites());
    ASSERT_EQ(w.get_Nbonds(), a.get_Nbonds());
    EXPECT_EQ(w.get_dim(), a.get_dim());
    EXPECT_EQ(w.get_nbasis(), a.get_nbasis());
    check_consistent(a);

    std::set<int> image;
    for (int s = 0; s < (int)w.get_Nsites(); ++s) image.insert(map(s));
    ASSERT_EQ(image.size(), a.get_Nsites()) << "the site map is not a bijection";

    auto const dirs = w.directions();
    for (int s = 0; s < (int)w.get_Nsites(); ++s) {
        int const as = map(s);
        EXPECT_EQ(w.get_zcoord(s), a.get_zcoord(as)) << "site " << s;
        for (size_t d = 0; d < WORM::zcmax; ++d) {
            if (!w.has_neighbor(s, d)) continue;
            int const at = map(w.nb(s, d));
            size_t e = 0;
            ASSERT_NO_THROW(e = a.get_direction(as, at)) << "bond " << s << " -- " << w.nb(s, d) << " missing";
            std::vector<double> wind(dirs[d].data(), dirs[d].data() + dirs[d].size());
            EXPECT_EQ(transform(wind), a.winding_vector(as, e))
                << "winding of hop " << s << " -> " << w.nb(s, d);
        }
    }

    if (same_numbering) {
        for (size_t d = 0; d < w.get_dim(); ++d) EXPECT_EQ(w.get_Ls(d), a.get_Ls(d));
        for (int s1 = 0; s1 < (int)w.get_Nsites(); ++s1)
            for (int s2 = 0; s2 < (int)w.get_Nsites(); ++s2)
                ASSERT_EQ(w.relNumbering(s1, s2), a.relNumbering(map(s1), map(s2)))
                    << "sites " << s1 << ", " << s2;
    }
}

Transform const same = [](std::vector<double> const& w) { return pad(w); };

} // namespace

TEST(PimcWormLattice, Chain) {
    chain w(6, 1, 1, true);
    alps_lattice const a = make_alps("chain lattice", 6);
    compare(w, a, identity_cells(a, {6, 1, 1}), same, true);
}

TEST(PimcWormLattice, OpenChain) {
    chain w(6, 1, 1, false);
    alps_lattice const a = make_alps("open chain lattice", 6);
    compare(w, a, identity_cells(a, {6, 1, 1}), same, true);
}

TEST(PimcWormLattice, Square) {
    square w(4, 5, 1, true, true);
    alps_lattice const a = make_alps("square lattice", 4, 5);
    compare(w, a, identity_cells(a, {4, 5, 1}), same, true);
}

TEST(PimcWormLattice, OpenSquare) {
    square w(4, 5, 1, false, false);
    alps_lattice const a = make_alps("open square lattice", 4, 5);
    compare(w, a, identity_cells(a, {4, 5, 1}), same, true);
}

TEST(PimcWormLattice, Cubic) {
    cubic w(3, 4, 5, true, true, true);
    alps_lattice const a = make_alps("simple cubic lattice", 3, 4, 5);
    compare(w, a, identity_cells(a, {3, 4, 5}), same, true);
}

TEST(PimcWormLattice, Ladder) {
    ladder w(5, 2, 1, true);
    alps_lattice const a = make_alps("ladder", 5);
    compare(w, a, identity_cells(a, {5, 2, 1}), same, true);
}

TEST(PimcWormLattice, Triangular) {
    triangular w(4, 5, 1, true, true);
    alps_lattice const a = make_alps("triangular lattice", 4, 5);
    compare(w, a, identity_cells(a, {4, 5, 1}), same, true);
}

// worm's honeycomb uses the primitive vectors a1' = a2 and a2' = a2 - a1 of
// the ALPS one. On an L x L torus both describe the same lattice; a worm cell
// (x, y) is the ALPS cell (-y, x + y), and windings transform alike.
TEST(PimcWormLattice, Honeycomb) {
    size_t const L = 4;
    honeycomb w(L, L, 1, true, true);
    alps_lattice const a = make_alps("honeycomb lattice", L);
    Map const map = site_map(a, {L, L, 1}, [](Cell c) { return Cell{-c[1], c[0] + c[1], 0}; });
    Transform const transform = [](std::vector<double> const& w) {
        return pad({-w[1], w[0] + w[1]});
    };
    compare(w, a, map, transform, false);
}

TEST(PimcWormLattice, RejectsParallelBonds) {
    EXPECT_THROW(make_alps("chain lattice", 2), std::runtime_error);
}

TEST(PimcWormLattice, NonBravaisLattice) {
    alps_lattice const a = make_alps("Kagome lattice", 3);
    EXPECT_EQ(a.get_nbasis(), 3u);
    EXPECT_EQ(a.get_Nsites(), 27u);
    for (int s = 0; s < (int)a.get_Nsites(); ++s) EXPECT_EQ(a.get_zcoord(s), 4u);
    check_consistent(a);
}
