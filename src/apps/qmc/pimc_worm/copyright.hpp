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
#include <iostream>
#include <alps/utilities/citations.hpp>

inline std::ostream& print_copyright(std::ostream& os) {
    os
    << "\n"
    << "# ALPS pimc_worm: path integral Monte Carlo with worm-type updates\n"
    << "# Copyright (C) 2021-2026  Jonas Greitemann, Nicolas Sadoune and Lode Pollet\n"
    << "# Part of the ALPS project, https://alps.comp-phys.org/\n"
    << "# Distributed under the MIT License; see LICENSE.txt in the ALPS sources.\n\n";
    alps::print_citations(os, "pimc_worm");

    return os;
}
