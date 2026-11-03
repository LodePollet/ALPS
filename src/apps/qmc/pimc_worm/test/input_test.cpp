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

// Checks alps_input.hpp: the arguments pyalps.runApplication adds, pyalps
// input files, and the translation of ALPS parameter names.

#include "../lattice.hpp"
#include "../alps_input.hpp"

#include <alps/hdf5/archive.hpp>
#include <alps/utilities/temporary_filename.hpp>

#include <gtest/gtest.h>

#include <cstdio>
#include <string>
#include <vector>

namespace {

struct args {
    explicit args(std::vector<std::string> a) : strings(std::move(a)) {
        for (std::string& s : strings) pointers.push_back(&s[0]);
    }
    int argc() { return int(pointers.size()); }
    char** argv() { return pointers.data(); }
    std::vector<std::string> strings;
    std::vector<char*> pointers;
};

std::vector<std::string> as_vector(pimc_worm_input::command_line& cl) {
    return std::vector<std::string>(cl.argv(), cl.argv() + cl.argc());
}

// A /parameters group as pyalps.writeInputH5Files writes it
std::string write_plain_input(const std::string& name) {
    alps::hdf5::archive ar(name, "w");
    ar["/parameters/LATTICE"] << std::string("square lattice");
    ar["/parameters/MODEL"] << std::string("boson Hubbard");
    ar["/parameters/L"] << 6;
    ar["/parameters/T"] << 0.05;
    ar["/parameters/t"] << 0.055;
    ar["/parameters/U"] << 1.0;
    ar["/parameters/Nmax"] << 2;
    ar["/parameters/SWEEPS"] << 600000;
    return name;
}

} // namespace

TEST(CommandLine, RewritesWhatPyalpsAdds) {
    args a({"pimc_worm_square_mpi", "--mpi", "job.ini", "-T", "100", "--Tmin", "5", "--TMax", "60", "--write-xml"});
    pimc_worm_input::command_line cl(a.argc(), a.argv(), true);
    EXPECT_EQ(std::vector<std::string>({"pimc_worm_square_mpi", "job.ini", "--runtimelimit=100", "--Tmin=5",
                                        "--Tmax=60"}),
              as_vector(cl));
    EXPECT_TRUE(cl.input_file().empty());
}

TEST(CommandLine, TheSerialProgramRefusesMpi) {
    args a({"pimc_worm_square", "--mpi", "job.ini"});
    EXPECT_THROW(pimc_worm_input::command_line(a.argc(), a.argv()), std::runtime_error);
}

TEST(CommandLine, TakesOutAPlainInputFile) {
    const std::string input = write_plain_input(alps::temporary_filename("pimc_worm_input") + ".task1.in.h5");
    args a({"pimc_worm_square", input, "--seed=3"});
    pimc_worm_input::command_line cl(a.argc(), a.argv());
    EXPECT_EQ(input, cl.input_file());
    EXPECT_EQ(std::vector<std::string>({"pimc_worm_square", "--seed=3"}), as_vector(cl));
    std::remove(input.c_str());
}

TEST(CommandLine, LeavesACheckpointToAlpsParams) {
    const std::string checkpoint = alps::temporary_filename("pimc_worm_input") + ".clone.h5";
    {
        alps::params p;
        p.define<int>("nmax", 2, "maximum occupation");
        alps::hdf5::archive ar(checkpoint, "w");
        ar["/parameters"] << p;
    }
    args a({"pimc_worm_square", checkpoint});
    pimc_worm_input::command_line cl(a.argc(), a.argv());
    EXPECT_TRUE(cl.input_file().empty());
    EXPECT_EQ(2, cl.argc());
    std::remove(checkpoint.c_str());
}

TEST(InputFile, SuppliesValuesAndNamesTheOutputAfterTheTask) {
    const std::string stem = alps::temporary_filename("pimc_worm_input") + ".task2";
    const std::string input = write_plain_input(stem + ".in.h5");
    args a({"pimc_worm_square", "--L=8"});
    alps::params p(a.argc(), a.argv());
    pimc_worm_input::supply_input_file(p, input);
    EXPECT_EQ("8", p.supplied_text().at("L"));  // the command line wins
    EXPECT_EQ("0.055", p.supplied_text().at("t"));
    EXPECT_EQ(stem + ".out.h5", p.supplied_text().at("outputfile"));
    EXPECT_EQ(stem + ".clone.h5", p.supplied_text().at("checkpoint"));
    std::remove(input.c_str());
}

TEST(AlpsNames, TranslatesToPimcWormNames) {
    args a({"pimc_worm_square", "--LATTICE=square lattice", "--MODEL=boson Hubbard", "--L=6", "--T=0.05",
            "--t=0.055", "--U=1", "--V=0.1", "--Nmax=2", "--SWEEPS=600000", "--THERMALIZATION=1000", "--SEED=7"});
    alps::params p(a.argc(), a.argv());
    pimc_worm_input::supply_alps_names<square>(p);
    p.define<std::string>("model", "", "")
     .define<double>("beta", 1., "").define<double>("t_hop", 1., "").define<double>("U_on", 4., "")
     .define<double>("V_nn", 0., "").define<int>("nmax", 10, "").define<int>("sweeps", 1, "")
     .define<int>("thermalization", 1, "").define<size_t>("seed", 1, "")
     .define<size_t>("Lx", 4, "").define<size_t>("Ly", 1, "");
    EXPECT_EQ("BoseHubbard", p["model"].as<std::string>());
    EXPECT_DOUBLE_EQ(20., p["beta"].as<double>());
    EXPECT_DOUBLE_EQ(0.055, p["t_hop"].as<double>());
    EXPECT_DOUBLE_EQ(1., p["U_on"].as<double>());
    EXPECT_DOUBLE_EQ(0.1, p["V_nn"].as<double>());
    EXPECT_EQ(2, p["nmax"].as<int>());
    EXPECT_EQ(600000, p["sweeps"].as<int>());
    EXPECT_EQ(1000, p["thermalization"].as<int>());
    EXPECT_EQ(7u, p["seed"].as<size_t>());
    EXPECT_EQ(6u, p["Lx"].as<size_t>());
    EXPECT_EQ(6u, p["Ly"].as<size_t>());
}

TEST(AlpsNames, ExtentsFollowTheLattice) {
    args a({"pimc_worm_ladder", "--L=5"});
    alps::params p(a.argc(), a.argv());
    pimc_worm_input::supply_alps_names<ladder>(p);
    EXPECT_EQ("5", p.supplied_text().at("Lx"));
    EXPECT_EQ(0u, p.supplied_text().count("Ly"));

    args b({"pimc_worm_cubic", "--L=3", "--W=4"});
    alps::params q(b.argc(), b.argv());
    pimc_worm_input::supply_alps_names<cubic>(q);
    EXPECT_EQ("3", q.supplied_text().at("Lx"));
    EXPECT_EQ("4", q.supplied_text().at("Ly"));
    EXPECT_EQ("3", q.supplied_text().at("Lz"));
}

TEST(AlpsNames, EqualBondTypeCouplingsAreOneCoupling) {
    args a({"pimc_worm_honeycomb", "--t0=1", "--t1=1", "--t2=1.0"});
    alps::params p(a.argc(), a.argv());
    pimc_worm_input::supply_alps_names<honeycomb>(p);
    EXPECT_EQ("1", p.supplied_text().at("t_hop"));
}

TEST(AlpsNames, RefusesWhatItCannotRepresent) {
    auto refuses = [](std::vector<std::string> v) {
        args a(std::move(v));
        alps::params p(a.argc(), a.argv());
        EXPECT_THROW(pimc_worm_input::supply_alps_names<square>(p), std::runtime_error) << a.strings[1];
    };
    refuses({"pimc_worm_square", "--T=0.1", "--beta=10"});              // both names
    refuses({"pimc_worm_square", "--t=0.1", "--t_hop=0.1"});
    refuses({"pimc_worm_square", "--t0=0.1", "--t1=0.2"});              // different bond types
    refuses({"pimc_worm_square", "--LATTICE=chain lattice"});           // another lattice
    refuses({"pimc_worm_square", "--MODEL=spin"});                      // another model
    refuses({"pimc_worm_square", "--L=4", "--Lx=4"});

    args h({"pimc_worm_honeycomb", "--L=4", "--W=3"});                 // differs from ALPS's torus
    alps::params p(h.argc(), h.argv());
    EXPECT_THROW(pimc_worm_input::supply_alps_names<honeycomb>(p), std::runtime_error);
}
