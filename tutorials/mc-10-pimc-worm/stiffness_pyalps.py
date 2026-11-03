# ****************************************************************************
#
# ALPS Project: Algorithms and Libraries for Physics Simulations
#
# Copyright (C) 2009-2010 by Matthias Troyer <troyer@phys.ethz.ch>
# Copyright (C) 2026 by Lode Pollet
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
#
# ****************************************************************************

"""Superfluid stiffness of the 2D Bose-Hubbard model with pimc_worm, run from pyalps.

tutorials/mc-05-bosons (tutorial5a.py, tutorial5b.py) with pimc_worm in place
of worm: the same ALPS parameters, written with writeInputH5Files instead of
writeInputFiles, since pimc_worm reads HDF5 rather than XML input.
stiffness_native.py runs the same with pimc_worm's own input files.

Usage: python stiffness_pyalps.py [a|b|ab]
"""

import sys

import matplotlib.pyplot as plt
import pyalps
import pyalps.plot


def parameters(Ls, ts, T, sweeps):
    return [{
        'LATTICE'        : 'square lattice',
        'MODEL'          : 'boson Hubbard',
        'T'              : T,
        'L'              : L,
        't'              : t,
        'U'              : 1.0,
        'mu'             : 0.5,
        'Nmax'           : 2,
        'THERMALIZATION' : 10000,
        'SWEEPS'         : sweeps,
        'runtimelimit'   : 3600,   # pimc_worm: the run ends after SWEEPS or this many seconds
        'SEED'           : 1000 + i,
    } for i, (L, t) in enumerate((L, t) for L in Ls for t in ts)]


def run(prefix, parms):
    input_files = pyalps.writeInputH5Files(prefix, parms)
    pyalps.runApplication('pimc_worm_square', input_files)
    data = pyalps.loadMeasurements(pyalps.getResultFiles(prefix=prefix), 'Stiffness')
    return pyalps.collectXY(data, x='t', y='Stiffness', foreach=['L'])


part = sys.argv[1] if len(sys.argv) > 1 else 'ab'

if 'a' in part:
    rhos = run('parm10a', parameters([4], [0.01 * i for i in range(1, 11)], T=0.1, sweeps=500000))
    plt.figure()
    pyalps.plot.plot(rhos)
    plt.xlabel('Hopping $t/U$')
    plt.ylabel('Superfluid density $\\rho _s$')
    plt.title('pimc_worm, $T = 0.1$, $4 \\times 4$')

if 'b' in part:
    rhos = run('parm10b', parameters([4, 6, 8], [0.045, 0.05, 0.0525, 0.055, 0.0575, 0.06, 0.065],
                                     T=0.05, sweeps=600000))
    # multiply with the system size for the scaling plot
    for s in rhos:
        s.y = s.y * float(s.props['L'])
    plt.figure()
    pyalps.plot.plot(rhos)
    plt.xlabel('Hopping $t/U$')
    plt.ylabel('$\\rho _sL$')
    plt.legend()
    plt.title('Scaling plot for Bose-Hubbard model (pimc_worm)')

plt.show()
