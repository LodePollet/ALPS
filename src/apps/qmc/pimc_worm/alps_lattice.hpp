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

#pragma once

// A lattice from the ALPS lattice library (alps::graph_helper), presented
// through the interface of the hand-written lattices in lattice.hpp, so that
// worm uses it unchanged. See docs/alps_lattice_design.md.

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <alps/lattice/graph_helper.h>
#include <alps/ngs/make_deprecated_parameters.hpp>
#include <alps/parameter.h>
#include <alps/params.hpp>

#include <Eigen/Dense>

#ifndef PIMC_WORM_ZCMAX
#define PIMC_WORM_ZCMAX 6
#endif

struct alps_lattice {
    using Base = alps_lattice;

    // Capacity of a winding vector; the lattice's own dimension is get_dim().
    static constexpr size_t dim = 3;
    // Maximum coordination number, which sizes every diagram element.
    static constexpr size_t zcmax = PIMC_WORM_ZCMAX;

    using SiteType = int;
    using SiteIndex = int;
    using BondIndex = size_t;
    using BondType = std::pair<size_t, size_t>;
    using DirectionIndex = size_t;
    using Direction = Eigen::Matrix<double, dim, 1>;

    static void define_parameters(alps::params& params) {
        params.define<std::string>("LATTICE", "chain lattice", "name of the lattice in the ALPS lattice library")
              .define<int>("L", 4, "linear size of the lattice (W, H and LATTICE_LIBRARY are passed on when supplied)");
    }

    static bool parameters_supplied(alps::params const& params) {
        return params.supplied("LATTICE") && params.supplied("L");
    }

    explicit alps_lattice(alps::params const& params)
        : alps_lattice(to_parameters(params)) {}

    explicit alps_lattice(alps::Parameters const& p) {
        alps::graph_helper<> g(p);
        name_ = static_cast<std::string>(p["LATTICE"]);
        dim_ = g.dimension();
        if (dim_ < 1 || dim_ > dim)
            throw std::runtime_error("pimc_worm supports lattices of dimension 1 to 3, not " + std::to_string(dim_));
        const size_t Nsites = g.num_sites();

        // Cell offsets of every site. Sites are numbered cell by cell, as in
        // lattice.hpp: site = cell index * n_basis + basis index.
        size_t Ncells = 0;
        extent_.assign(dim, 1);
        periodic_.assign(dim, false);
        for (size_t d = 0; d < dim_; ++d) {
            extent_[d] = g.lattice().extent(d);
            periodic_[d] = g.lattice().boundary(d) == "periodic";
        }
        for (auto c = g.cells(); c.first != c.second; ++c.first) ++Ncells;
        if (Ncells == 0 || Nsites % Ncells != 0)
            throw std::runtime_error("pimc_worm cannot use a depleted lattice");
        n_basis_ = Nsites / Ncells;
        offset_.assign(Nsites, std::array<int, dim>{0, 0, 0});
        for (auto c = g.cells(); c.first != c.second; ++c.first) {
            auto const off = alps::offset(*c.first, g.lattice());
            size_t const cell = alps::index(*c.first, g.lattice());
            for (size_t b = 0; b < n_basis_; ++b)
                for (size_t d = 0; d < dim_; ++d)
                    offset_[cell * n_basis_ + b][d] = off[d];
        }

        // Neighbor tables in bond-index order: the direction index is the
        // slot in a site's neighbor list, filled from both ends of each bond.
        std::vector<std::pair<size_t, std::pair<size_t, size_t>>> bonds;
        std::vector<std::array<int, dim>> crossing;
        for (auto b = g.bonds(); b.first != b.second; ++b.first) {
            size_t const s = g.source(*b.first), t = g.target(*b.first);
            if (s == t)
                throw std::runtime_error("pimc_worm cannot use a lattice with a bond from a site to itself");
            alps::boundary_crossing const bc = boost::get(alps::boundary_crossing_t(), g.graph(), *b.first);
            std::array<int, dim> cr{0, 0, 0};
            for (size_t d = 0; d < dim_; ++d) cr[d] = bc.crosses(d);
            bonds.push_back({static_cast<size_t>(g.index(*b.first)), {s, t}});
            crossing.push_back(cr);
        }
        std::vector<size_t> order(bonds.size());
        for (size_t i = 0; i < order.size(); ++i) order[i] = i;
        std::sort(order.begin(), order.end(), [&](size_t i, size_t j) { return bonds[i].first < bonds[j].first; });

        sites_.assign(Nsites, std::vector<neighbor>());
        for (size_t const i : order) {
            size_t const s = bonds[i].second.first, t = bonds[i].second.second;
            BondType const key = s < t ? BondType(s, t) : BondType(t, s);
            if (!BondIndices_.emplace(key, BondIndices_.size()).second)
                throw std::runtime_error("pimc_worm cannot use a lattice with two bonds between sites "
                                         + std::to_string(s) + " and " + std::to_string(t)
                                         + " (is a periodic extent smaller than 3?)");
            // The winding of a hop s -> t is its displacement in units of the
            // primitive vectors, unwrapped across periodic boundaries: the
            // difference of cell offsets plus the extent for every crossing.
            Direction w = Direction::Zero();
            for (size_t d = 0; d < dim_; ++d)
                w[d] = offset_[t][d] - offset_[s][d] + crossing[i][d] * static_cast<int>(extent_[d]);
            add_neighbor(s, t, BondIndices_.at(key), w);
            add_neighbor(t, s, BondIndices_.at(key), -w);
        }
        for (size_t s = 0; s < Nsites; ++s)
            if (sites_[s].size() > zcmax)
                throw std::runtime_error("site " + std::to_string(s) + " has " + std::to_string(sites_[s].size())
                                         + " neighbors, more than PIMC_WORM_ZCMAX = " + std::to_string(zcmax));
    }

    size_t get_dim() const {return dim_;}
    size_t get_nbasis() const {return n_basis_;}
    size_t get_zcoord(SiteIndex s) const {return sites_.at(s).size();}
    size_t get_Nsites() const {return sites_.size();}
    size_t get_Nbonds() const {return BondIndices_.size();}
    size_t get_N_nb(SiteIndex s) const {return sites_.at(s).size();}
    size_t get_Ls(int d) const {
      if ((d < 0) || (d > 2)) throw std::runtime_error("call to get_Ls() with invalid argument");
      return extent_[d];
    }
    bool is_periodic(size_t d) const {return periodic_.at(d);}
    // Offset of the cell of site s in units of the primitive vectors.
    std::array<int, dim> const& cell_offset(SiteIndex s) const {return offset_.at(s);}

    SiteType operator[](SiteIndex s) const {return s;}
    SiteType get_site(SiteIndex s) const {return s;}

    bool has_neighbor(SiteIndex s, DirectionIndex dir) const {return dir < sites_.at(s).size();}
    SiteType nb(SiteIndex s, DirectionIndex dir) const {return sites_.at(s).at(dir).site;}

    DirectionIndex get_direction(SiteIndex s1, SiteIndex s2) const {
        auto const& n = sites_.at(s1);
        for (DirectionIndex d = 0; d < n.size(); ++d)
            if (n[d].site == s2) return d;
        throw std::out_of_range("sites " + std::to_string(s1) + " and " + std::to_string(s2) + " are not neighbors");
    }
    DirectionIndex get_opposite_direction(SiteIndex s1, SiteIndex s2) const {return get_direction(s2, s1);}
    BondIndex get_bond_index(SiteIndex s, DirectionIndex dir) const {return sites_.at(s).at(dir).bond;}

    // Winding of a hop from s along direction dir, as directions()[dir] in lattice.hpp.
    Direction winding_vector(SiteIndex s, DirectionIndex dir) const {return sites_.at(s).at(dir).winding;}

    // Relative numbering of the cell of s2 w.r.t. the cell of s1, as in lattice.hpp.
    // This is NOT symmetric in s1, s2, i.e. relNumbering(s1,s2) != relNumbering(s2,s1)
    size_t relNumbering(SiteType s1, SiteType s2) const {
        size_t r = 0, stride = 1;
        for (size_t d = 0; d < dim_; ++d) {
            int const L = static_cast<int>(extent_[d]);
            r += stride * static_cast<size_t>(((offset_[s2][d] - offset_[s1][d]) % L + L) % L);
            stride *= extent_[d];
        }
        return r;
    }

    void print_name(std::ostream& os) const {
        os << "# lattice: Name                                     : " << name_ << " (ALPS lattice library)\n";
    }

    void print() const {
        std::cout << "# Nsites, dim, zcoord " << sites_.size() << " " << dim_ << " "
                  << sites_.at(0).size() << "\n";
        for (size_t s = 0; s < sites_.size(); ++s) {
            std::cout << "site: " << s << "\t\t" << "has neighbors: ";
            for (size_t d = 0; d < sites_[s].size(); ++d)
                std::cout << sites_[s][d].site << "(dir " << d << ")\t";
            std::cout << std::endl;
        }
        for (auto const& b : BondIndices_){
            std::cout << "bond: " << b.first.first << "," << b.first.second << "\t\t" << "has index: "
                        << b.second << std::endl;
        }
    }

    void print_params(std::ostream& os) const {
        os << "# lattice: dim                                      : " << dim_ << "\n";
        for (size_t d = 0; d < dim_; ++d)
            os << "# lattice: periodic in direction " << d << "                : " << periodic_[d] << "\n";
        os << "# lattice: length                                   : " << "L, W, H  ";
        for (size_t d = 0; d < dim; ++d) os << extent_[d] << "\t";
        os << "num Sites  " << sites_.size() << "\t"
           << "num Bonds  "  << BondIndices_.size() << "\n";
    }

    // alps::params as the ALPS lattice library reads them: the supplied text,
    // with LATTICE and L also when they take their default.
    static alps::Parameters to_parameters(alps::params const& params) {
        alps::Parameters p = alps::make_deprecated_parameters(params);
        p["LATTICE"] = params["LATTICE"].as<std::string>();
        p["L"] = std::to_string(params["L"].as<int>());
        return p;
    }

private:
    struct neighbor {
        SiteIndex site;
        BondIndex bond;
        Direction winding;
    };

    void add_neighbor(size_t s, size_t t, BondIndex bond, Direction const& w) {
        sites_[s].push_back({static_cast<SiteIndex>(t), bond, w});
    }

    std::string name_;
    size_t dim_;
    size_t n_basis_;
    std::vector<size_t> extent_;
    std::vector<bool> periodic_;
    std::vector<std::array<int, dim>> offset_;
    std::vector<std::vector<neighbor>> sites_;
    std::map<BondType, BondIndex> BondIndices_; //without double counting
};
