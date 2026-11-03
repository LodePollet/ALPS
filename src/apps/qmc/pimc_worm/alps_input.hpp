/*****************************************************************************
*
* ALPS Project Applications
*
* Copyright (C) 2026 by Lode Pollet
*
* ALPS Project: https://alps.comp-phys.org/
* SPDX-License-Identifier: MIT
*
*****************************************************************************/

// Input in the conventions of ALPS's other programs, so that pimc_worm runs
// from pyalps as well as from its own .ini files:
//   - the arguments pyalps.runApplication adds (-T, --Tmin, --TMax,
//     --write-xml, --mpi),
//   - an input file written by pyalps.writeInputH5Files (a /parameters group
//     without ALPSCore's attributes), as opposed to a checkpoint,
//   - the ALPS parameter names (T, t, U, V, mu, Nmax, MODEL, LATTICE, L, W,
//     H, SWEEPS, THERMALIZATION, SEED), translated to pimc_worm's own.

#ifndef PIMC_WORM_ALPS_INPUT_HPP
#define PIMC_WORM_ALPS_INPUT_HPP

#include <alps/hdf5/archive.hpp>
#include <alps/params.hpp>
#include <alps/params/text_values.hpp>

#include <map>
#include <regex>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <vector>

namespace pimc_worm_input {

// True for an HDF5 file whose /parameters group lacks ALPSCore's attributes,
// as pyalps.writeInputH5Files writes it; a checkpoint has them.
inline bool is_plain_input(const std::string& name) {
    if (name.size() < 3 || name.compare(name.size() - 3, 3, ".h5") != 0) return false;
    try {
        alps::hdf5::archive ar(name, "r");
        return ar.is_group("/parameters") && !ar.is_attribute("/parameters@status");
    }
    catch (const std::exception&) {
        return false;
    }
}

// The command line with pyalps's additions rewritten for alps::params:
// "-T n" is the run time limit, "--Tmin n"/"--TMax n" become "--Tmin=n" and
// "--Tmax=n" (read by the MPI driver), --write-xml is dropped, and a plain
// input file is taken out to be supplied separately. --mpi, which
// runApplication(..., MPI=n) adds, is dropped by the MPI driver and refused
// by the serial one, whose copies would all write the same files.
class command_line {
public:
    command_line(int argc, char** argv, bool mpi_driver = false) {
        args_.push_back(argv[0]);
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            const bool has_value = i + 1 < argc;
            if (arg == "-T" && has_value) {
                args_.push_back("--runtimelimit=" + std::string(argv[++i]));
            }
            else if (arg == "--Tmin" && has_value) {
                args_.push_back("--Tmin=" + std::string(argv[++i]));
            }
            else if ((arg == "--TMax" || arg == "--Tmax") && has_value) {
                args_.push_back("--Tmax=" + std::string(argv[++i]));
            }
            else if (arg == "--mpi") {
                if (!mpi_driver)
                    throw std::runtime_error("--mpi: this program is serial; run pimc_worm_<lattice>_mpi "
                                             "under mpirun instead");
            }
            else if (arg == "--write-xml") {
                continue;
            }
            else if (input_file_.empty() && is_plain_input(arg)) {
                input_file_ = arg;
            }
            else {
                args_.push_back(arg);
            }
        }
        for (const std::string& a : args_) argv_.push_back(const_cast<char*>(a.c_str()));
        argv_.push_back(nullptr);
    }

    int argc() const { return int(args_.size()); }
    char** argv() { return argv_.data(); }
    const std::string& input_file() const { return input_file_; }

private:
    std::vector<std::string> args_;
    std::vector<char*> argv_;
    std::string input_file_;
};

namespace detail {

inline bool has(const alps::params& p, const std::string& name) {
    return p.supplied_text().count(name) != 0;
}

inline const std::string& text(const alps::params& p, const std::string& name) {
    return p.supplied_text().at(name);
}

// Supplies `native` from the ALPS name `alps` unless `native` is supplied
// too, in which case the two must not both be given.
inline void alias(alps::params& p, const std::string& alps, const std::string& native) {
    if (!has(p, alps)) return;
    if (has(p, native))
        throw std::runtime_error("Both " + alps + " (ALPS name) and " + native + " are given; give one of them");
    p.supply(native, text(p, alps));
}

// A bond coupling: `name` or the bond-type values name0, name1, ..., which
// pimc_worm accepts only when they are all equal (one coupling for every bond).
inline void bond_coupling(alps::params& p, const std::string& name, const std::string& native) {
    const std::regex typed(name + "[0-9]+");
    std::vector<std::string> values;
    if (has(p, name)) values.push_back(text(p, name));
    for (const auto& kv : p.supplied_text()) {
        if (std::regex_match(kv.first, typed)) values.push_back(kv.second);
    }
    if (values.empty()) return;
    for (const std::string& v : values) {
        if (std::stod(v) != std::stod(values.front()))
            throw std::runtime_error("Bond-type couplings " + name + "0, " + name + "1, ... must all be equal; "
                                     "pimc_worm has one " + native + " for every bond");
    }
    if (has(p, native))
        throw std::runtime_error("Both " + name + " or " + name + "0, " + name + "1, ... (ALPS names) and " + native
                                 + " are given; give one of them");
    p.supply(native, values.front());
}

} // namespace detail

// Supplies the values of a plain input file. Values given on the command
// line take precedence. The output file is named after the input,
// prefix.taskN.in.h5 -> prefix.taskN.out.h5, as pyalps.getResultFiles
// expects, and the checkpoint prefix.taskN.clone.h5, one per task (the
// default would drop every extension and give all tasks prefix.clone.h5).
inline void supply_input_file(alps::params& p, const std::string& input_file) {
    if (input_file.empty()) return;
    alps::params plain;
    {
        alps::hdf5::archive ar(input_file, "r");
        ar.set_context("/parameters");
        alps::params_ns::supply_plain_group(ar, plain);
    }
    for (const auto& kv : plain.supplied_text()) {
        if (!detail::has(p, kv.first)) p.supply(kv.first, kv.second);
    }
    std::string stem = input_file.substr(0, input_file.size() - 3);  // without .h5
    const std::string in_ext = ".in";
    if (stem.size() > in_ext.size() && stem.compare(stem.size() - in_ext.size(), in_ext.size(), in_ext) == 0)
        stem.erase(stem.size() - in_ext.size());
    if (!detail::has(p, "outputfile")) p.supply("outputfile", stem + ".out.h5");
    if (!detail::has(p, "checkpoint")) p.supply("checkpoint", stem + ".clone.h5");
}

// Translates the ALPS parameter names to pimc_worm's. Call before
// define_parameters() on parameters that are not restored from a checkpoint.
template <typename Lattice>
void supply_alps_names(alps::params& p) {
    using detail::has;
    using detail::text;

    if (has(p, "MODEL")) {
        if (text(p, "MODEL") != "boson Hubbard")
            throw std::runtime_error("MODEL=\"" + text(p, "MODEL") + "\" is not available with ALPS names; "
                                     "pimc_worm takes MODEL=\"boson Hubbard\", or model=\"XXZ\" with Jpm, Jzz, h");
        if (!has(p, "model")) p.supply("model", "BoseHubbard");
    }
    if (has(p, "T")) {
        if (has(p, "beta"))
            throw std::runtime_error("Both T (ALPS name) and beta are given; give one of them");
        p.supply("beta", alps::params_ns::to_text(1.0 / std::stod(text(p, "T"))));
    }
    detail::bond_coupling(p, "t", "t_hop");
    detail::bond_coupling(p, "V", "V_nn");
    detail::alias(p, "U", "U_on");
    detail::alias(p, "Nmax", "nmax");
    detail::alias(p, "SWEEPS", "sweeps");
    detail::alias(p, "THERMALIZATION", "thermalization");
    detail::alias(p, "SEED", "seed");

#ifndef PIMC_WORM_ALPS_LATTICE
    // pimc_worm_alps reads LATTICE, L, W, H itself; here they must describe
    // this executable's lattice, with periodic boundaries.
    if (has(p, "LATTICE") && text(p, "LATTICE") != Lattice::alps_name)
        throw std::runtime_error("LATTICE=\"" + text(p, "LATTICE") + "\" is not the lattice of this program (\""
                                 + Lattice::alps_name + "\"); use pimc_worm_alps for other ALPS lattices");
    if (has(p, "L")) {
        const std::string L = text(p, "L");
        const std::string W = has(p, "W") ? text(p, "W") : L;
        const std::string H = has(p, "H") ? text(p, "H") : L;
        if (std::is_same<Lattice, honeycomb>::value && W != L)
            throw std::runtime_error("pimc_worm_honeycomb agrees with the ALPS honeycomb lattice only for W = L; "
                                     "use pimc_worm_alps");
        std::map<std::string, std::string> extents{{"Lx", L}};
        if (Lattice::dim >= 2 && !std::is_same<Lattice, ladder>::value) extents["Ly"] = W;
        if (Lattice::dim >= 3) extents["Lz"] = H;
        for (const auto& e : extents) {
            if (has(p, e.first))
                throw std::runtime_error("Both L/W/H (ALPS names) and " + e.first + " are given; give one of them");
            p.supply(e.first, e.second);
        }
    }
#endif
}

} // namespace pimc_worm_input

#endif // PIMC_WORM_ALPS_INPUT_HPP
