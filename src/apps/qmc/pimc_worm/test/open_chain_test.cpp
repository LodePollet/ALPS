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

// Runs pimc_worm on an open four-site soft-core Bose-Hubbard chain and
// compares with exact diagonalization. The end sites have one neighbor and
// the others two, which exercises the coordination-number factors of
// INSERTKINK/DELETEKINK: before the DELETEKINK fix the kinetic energy here
// was -1.620 instead of -1.548, about 20 standard errors at this run length.
//
// H = -t sum_i (b_i^+ b_{i+1} + h.c.) + U/2 sum_i n_i (n_i - 1) - mu sum_i n_i
// on sites 0..3, n_i <= 3, t = 0.3, U = 1, mu = 0.5, T = 0.5 (the reproducer
// of ALPSim/ALPS#295); the exact grand-canonical averages below come from a
// full diagonalization of the 256 states.

#include "../worm.hpp"

#include <alps/mc/api.hpp>
#include <alps/mc/stop_callback.hpp>

#include <gtest/gtest.h>

#include <cmath>
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

} // namespace

TEST(OpenChain, SoftCoreBosonsAgreeWithExactDiagonalization) {
    args a({"pimc_worm_open_chain_test",
            "--model=BoseHubbard", "--Lx=4", "--Ly=1", "--Lz=1", "--pbcx=0", "--pbcy=0", "--pbcz=0",
            "--t_hop=0.3", "--U_on=1.0", "--V_nn=0.0", "--mu=0.5", "--nmax=3", "--beta=2.0",
            "--thermalization=100000", "--sweeps=3000000", "--runtimelimit=600", "--seed=2026",
            "--Ntest=1000000", "--Nsave=1000000000000"});  // no checkpoint files
    alps::params parameters(a.argc(), a.argv());
    worm::define_parameters(parameters);
    ASSERT_FALSE(parameters.has_missing());

    worm sim(parameters);
    sim.initialize();
    sim.run(alps::stop_callback(size_t(parameters["runtimelimit"])));
    ASSERT_GE(sim.fraction_completed(), 1.) << "the run stopped on its time limit";
    sim.test_conf();

    const auto results = alps::collect_results(sim);
    const struct { const char* name; double exact; } expected[] = {
        {"Kinetic_Energy",      -1.548107},
        {"Potential_Energy",    -0.792633},
        {"Total_Energy",        -2.340740},
        {"Number_of_particles",  4.915584},
    };
    for (const auto& e : expected) {
        const double mean = results[e.name].mean<double>();
        const double error = results[e.name].error<double>();
        EXPECT_GT(error, 0.) << e.name;
        EXPECT_LT(error, 0.02) << e.name << ": the run is too short to resolve the bias it guards against";
        EXPECT_LE(std::abs(mean - e.exact), 4 * error)
            << e.name << " = " << mean << " +/- " << error << ", exact " << e.exact;
    }
}
