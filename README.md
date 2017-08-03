# `pymbd` — Many-body dispersion method

Python 2/3 package for calculating [many-body dispersion](http://dx.doi.org/10.1063/1.4865104) energies. Also contains corresponding ScaLAPACK implementation of main routines.

Most functionality is implemented in the Fortran module in `src/mbd.f90`.

## Installation

Installation of Pymbd requires Numpy with properly configured Lapack libraries.

On macOS, this requirement is almost always satisfied with whatever Python with simple

```
python/python3 -m pip install numpy
```

This can be followed with

```
python/python3 -m pip install git+https://github.com/azag0/pymbd.git
```

You can run tests with

```
python/python3 -m unittest pymbd.tests -v
```

On Linux, the pip-installed Numpy sometimes doesn't report the correct Lapack [1]. To avoid this issue, use the Anaconda Python with the no-MKL Numpy (if you have Anaconda already installed, [make sure](https://www.continuum.io/blog/developer-blog/anaconda-25-release-now-mkl-optimizations) that you don't use the MKL Numpy):

```
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
[...]
<path to anaconda>/bin/conda install nomkl numpy
```

To install Pymbd, run

```
<path to anaconda>/bin/pip install git+https://github.com/azag0/pymbd.git
```

And finally test with

```
<path to anaconda>/bin/python -m unittest pymbd.tests -v
```

[1]: What really matters is whether `numpy.distutils.system_info.get_info('lapack_opt')` returns the correct information. On macOS, it always returns the Accelerate framework. With Anaconda and no-MKL Numpy, it always returns the internal openBLAS libraries. But as long as it returns a valid Lapack library, you can use whatever Python installation.

## Usage

Pymbd doesn't have any input files, it is called directly from Python scripts. 

For documentation, see the [basic examples](http://nbviewer.jupyter.org/github/azag0/mbd/blob/master/examples/basic.ipynb) or the [more advanced use](http://nbviewer.jupyter.org/github/azag0/mbd/blob/master/examples/advanced.ipynb).

## ScaLAPACK

The ScaLAPACK version currently uses a different build system, not compatible with setuptools. In order to make use of the ScaLAPACK version, run `make mbd` after adapting `system.example.mk` for your system, saving it as `system.mk` and making sure you have `f2py` installed.

ScaLAPACK routines have the suffix `_s`, e.g. `mbd.get_mbd_energy_s(...)` instead of `mbd.get_mbd_energy(...)`. However, original MPI/LAPACK routines (such as `get_mbd_energy`) are still available. Hence, `mbd.get_mbd_energy(...)` will call MPI/LAPACK routine and `mbd.get_mbd_energy_s(...)` calls the ScaLAPACK version of it.

Latest version support use of different (Sca)LAPACK eigensolvers ("QR", "Divide and Conquer", "MR3")

Example in python (after compiling corresponding MBD module):

. add path to mbd and pymbd module to your system's PYTHONPATH

`
from mbd(_scalapack) import mbd(_scalapack) as mbd
from pymbd import get_free_atom_data (also change import mbd in pymbd if neccessary)
from mpi4py import MPI

mbd.init_grid(<n_omega_SCS>)
mbd.param_vacuum_axis = <vacuum_axis(x,y,z)>
mbd.my_task = MPI.COMM_WORLD.Get_rank()
mbd.n_tasks = MPI.COMM_WORLD.Get_size()
mbd.eigensolver = <'qr'/'dandc'/'mrrr'>
alpha_free, C6_free, RvdW_free = get_free_atom_data(<symbols>)

... rescaling of dispersion parameters ~> alpha_0_eff, C6_eff, RvdW_eff ...

omega = mbd.omega_eff(C6_eff, alpha_0_eff)
alpha_dyn_eff = mbd.alpha_dynamic_ts_all('C', \
                                        mbd.n_grid_omega, \
                                        <alpha_0_eff>, \
                                        c6=<C6_eff>)

alpha_dyn_SCS = mbd.run_scs_s(<modus_scs>, \
                              <Coulomb_SCS>, \
                              <pos>, \
                              <alpha_dyn_eff>, \
                              r_vdw=<RvdW_eff>, \
                              beta=<damping_param_beta>, \
                              a=<damping_param_d>, \
                              unit_cell=<UC>)

alpha_0_SCS = alpha_dyn_SCS[0]
E_MBD = mbd.get_mbd_energy_s(<modus>, \
                             Coulomb_CFDM, \
                             <pos>, <alpha_0_SCS>, <omega_SCS>, \
                             supercell=<supercell>, \
                             k_grid=<kgrid>, \
                             unit_cell=<UC>, \
                             r_vdw=<rvdwAB>, \
                             beta=<damping_param_beta>, \
                             a=<damping_param_d>)

mbd.destroy_grid()
`

### Format of eigenmode/eigenenergy output
(calling top-level `get_mbd_energy` instance in LAPACK and all related routines in ScaLAPACK framework create output files)

. eigenenergies: formatted text file `mbd_eigenvalues(_reciprocal).out` sorted by k points (if applicable)

. eigenmodes: unformatted binary file `mbd_eigenmodes(_kptX).out` sorted by k points (if applicable)
    - `mbd_eigenmodes.out`: size 1st dim, size 2nd dim, elements of 1st eigenmode, elements of 2nd eigenmode, ...
    - `mbd_eigenmodes_kptX.out`: size 1st dim, size 2nd dim, kpointX(1), kpointX(2), kpointX(3), elements of 1st eigenmode on kpoint X, elements of 2nd eigenmode on kpoint X, ...


