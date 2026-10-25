# Copyright (C) 2021-2026 by Jonas Greitemann, Nicolas Sadoune and Lode Pollet
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT

# Path integral Monte Carlo with worm-type updates for Bose-Hubbard and XXZ
# models, imported from https://github.com/LodePollet/worm; see README.md.

# The lattice sets the maximum coordination number that sizes every diagram
# element, so it is chosen at compile time: one executable per lattice.
set(PIMC_WORM_LATTICES chain square cubic ladder triangular honeycomb
  CACHE STRING "Lattices for which pimc_worm executables are built")
set(PIMC_WORM_DSTRUC LIST CACHE STRING
  "pimc_worm diagram data structure: AVL, LIST or LIST_STACK")
set_property(CACHE PIMC_WORM_DSTRUC PROPERTY STRINGS AVL LIST LIST_STACK)
option(PIMC_WORM_UNIFORM "pimc_worm: uniform site and bond parameters" ON)
option(PIMC_WORM_CWINDOW "pimc_worm: canonical window for the Green's function" OFF)
option(PIMC_WORM_MATSUBARA "pimc_worm: G(k,omega_n) and binned G(k=0,tau) measurements" OFF)
option(PIMC_WORM_DEBUG "pimc_worm: check the configuration after every update" OFF)

# lattice.hpp uses Eigen, which ALPS treats as optional (see src/alps/CMakeLists.txt).
find_package(Eigen3 CONFIG QUIET)
if(NOT Eigen3_FOUND OR Eigen3_VERSION VERSION_LESS 3.3)
  message(STATUS "Eigen 3.3 or newer not found: pimc_worm is not built")
  return()
endif()

if(NOT PIMC_WORM_DSTRUC MATCHES "^(AVL|LIST|LIST_STACK)$")
  message(FATAL_ERROR "Unknown PIMC_WORM_DSTRUC ${PIMC_WORM_DSTRUC}")
endif()
set(_pimc_worm_definitions DSTRUC_${PIMC_WORM_DSTRUC})
if(PIMC_WORM_UNIFORM)
  list(APPEND _pimc_worm_definitions UNISYS)
endif()
if(PIMC_WORM_CWINDOW)
  list(APPEND _pimc_worm_definitions CAN_WINDOW)
endif()
if(PIMC_WORM_MATSUBARA)
  list(APPEND _pimc_worm_definitions MATSUBARA_MEAS)
endif()
if(PIMC_WORM_DEBUG)
  list(APPEND _pimc_worm_definitions DEBUGMODE)
endif()

set(_pimc_worm_targets "")
foreach(lattice IN LISTS PIMC_WORM_LATTICES)
  if(NOT lattice MATCHES "^(chain|square|cubic|ladder|triangular|honeycomb)$")
    message(FATAL_ERROR "Unknown pimc_worm lattice ${lattice}")
  endif()
  set(impl pimc_worm_${lattice}_impl)
  add_library(${impl} STATIC worm.cpp worm.update.cpp worm.output.cpp)
  target_compile_definitions(${impl} PUBLIC LATTICE=${lattice} ${_pimc_worm_definitions})
  target_link_libraries(${impl} PUBLIC ALPS::alps Eigen3::Eigen)

  add_executable(pimc_worm_${lattice} worm.run.cpp)
  target_link_libraries(pimc_worm_${lattice} PRIVATE ${impl})
  list(APPEND _pimc_worm_targets pimc_worm_${lattice})

  if(ALPS_ENABLE_MPI)
    add_executable(pimc_worm_${lattice}_mpi worm.run_mpi.cpp)
    target_link_libraries(pimc_worm_${lattice}_mpi PRIVATE ${impl})
    list(APPEND _pimc_worm_targets pimc_worm_${lattice}_mpi)
  endif()
endforeach()
alps_install_applications(${_pimc_worm_targets})
