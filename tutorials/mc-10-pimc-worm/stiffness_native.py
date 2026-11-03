# ****************************************************************************
#
# ALPS Project: Algorithms and Libraries for Physics Simulations
#
# Copyright (C) 2026 by Lode Pollet
#
# ALPS Project: https://alps.comp-phys.org/
# SPDX-License-Identifier: MIT
#
# ****************************************************************************

"""Superfluid stiffness of the 2D Bose-Hubbard model with pimc_worm, run directly.

The physics of tutorials/mc-05-bosons, H = -t sum_<ij> b_i^+ b_j
+ U/2 sum_i n_i (n_i - 1) - mu sum_i n_i on a square lattice, with U = 1,
mu = 0.5 and at most two bosons per site, run with pimc_worm's own input files:
one .ini file per point, results read from the .out.h5 file with h5py.
stiffness_pyalps.py runs the same with pyalps.

  a: rho_s against t at T = 0.1 on a 4x4 lattice (tutorial 5a);
  b: rho_s L against t at T = 0.05 for L = 4, 6, 8, where the curves cross at
     the superfluid-Mott transition (tutorial 5b).

pimc_worm_square is taken from PIMC_WORM_BIN, else from PATH. Points run in
parallel, one per core.
"""

import argparse
import os
import shutil
import subprocess
from concurrent.futures import ThreadPoolExecutor

import h5py
import matplotlib.pyplot as plt

PROGRAM = 'pimc_worm_square'


def executable():
    directory = os.environ.get('PIMC_WORM_BIN')
    path = os.path.join(directory, PROGRAM) if directory else shutil.which(PROGRAM)
    if not path or not os.path.isfile(path):
        raise SystemExit(PROGRAM + ' not found: set PIMC_WORM_BIN or add it to PATH')
    return path


def ini_text(L, t, T, sweeps, seed):
    parms = {
        'model': '"BoseHubbard"',
        'Lx': L, 'Ly': L, 'Lz': 1,
        'pbcx': 1, 'pbcy': 1, 'pbcz': 0,
        't_hop': t, 'U_on': 1.0, 'V_nn': 0.0, 'mu': 0.5, 'nmax': 2,
        'beta': 1.0 / T,
        'thermalization': 10000, 'sweeps': sweeps,
        'runtimelimit': 3600,   # the run ends after `sweeps` or this many seconds
        'seed': seed,
    }
    return ''.join('%s = %s\n' % kv for kv in parms.items())


def run(task):
    """Runs one point in its own directory; returns (L, t, rho_s, error)."""
    program, directory, L, t, T, sweeps, seed = task
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, 'job.ini'), 'w') as f:
        f.write(ini_text(L, t, T, sweeps, seed))
    with open(os.path.join(directory, 'job.log'), 'w') as log:
        subprocess.run([program, 'job.ini'], cwd=directory, stdout=log, stderr=subprocess.STDOUT, check=True)
    with h5py.File(os.path.join(directory, 'job.out.h5'), 'r') as f:
        rho = f['simulation/results/Stiffness/mean']
        return L, t, rho['value'][()], rho['error'][()]


def scan(name, Ls, ts, T, sweeps):
    program = executable()
    tasks = [(program, os.path.join('runs', name, 'L%d_t%g' % (L, t)), L, t, T, sweeps, 1000 + i)
             for i, (L, t) in enumerate((L, t) for L in Ls for t in ts)]
    with ThreadPoolExecutor(max_workers=os.cpu_count()) as pool:
        return list(pool.map(run, tasks))


def plot(results, scale_by_L, ylabel, title):
    plt.figure()
    for L in sorted({r[0] for r in results}):
        points = sorted(r for r in results if r[0] == L)
        factor = [L if scale_by_L else 1 for _ in points]
        plt.errorbar([p[1] for p in points], [f * p[2] for f, p in zip(factor, points)],
                     yerr=[f * p[3] for f, p in zip(factor, points)], marker='o', label='L = %d' % L)
    plt.xlabel('Hopping $t/U$')
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument('part', nargs='?', choices=['a', 'b', 'ab'], default='ab')
    args = parser.parse_args()

    if 'a' in args.part:
        results = scan('5a', [4], [0.01 * i for i in range(1, 11)], T=0.1, sweeps=500000)
        for L, t, rho, err in results:
            print('a: L = %d  t = %-6g rho_s = %.5f +/- %.5f' % (L, t, rho, err))
        plot(results, False, 'Superfluid density $\\rho_s$', 'pimc_worm, $T = 0.1$, $4 \\times 4$')
    if 'b' in args.part:
        results = scan('5b', [4, 6, 8], [0.045, 0.05, 0.0525, 0.055, 0.0575, 0.06, 0.065], T=0.05, sweeps=600000)
        for L, t, rho, err in results:
            print('b: L = %d  t = %-6g rho_s L = %.4f +/- %.4f' % (L, t, rho * L, err * L))
        plot(results, True, '$\\rho_s L$', 'pimc_worm, scaling plot at $T = 0.05$')
    plt.show()


if __name__ == '__main__':
    main()
