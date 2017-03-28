# `mbd` — Many-body dispersion method

Implementation of the MBD method in Fortran, wrapped for Python with `f2py`.
Also contains corresponding ScaLAPACK implementation of main routines.

### Installation

Adapt `system.example.mk` for your system, save it as `system.mk` and make sure you have `f2py` installed. Then run `make`.

In order to make use of the ScaLAPACK version run `make mbd_scalapack`.

### Usage

See the Jupyter notebook in `playground.ipynb`.

ScaLAPACK routines have the suffix `_s`, e.g. `mbd_scalapack.get_mbd_energy_s(...)` instead of `mbd.get_mbd_energy(...)`. However, original MPI/LAPACK routines (such as `get_mbd_energy`) are still available in `mbd_scalapack`. Hence, `mbd_scalapack.get_mbd_energy(...)` will call MPI/LAPACK routine and `mbd_scalapack.get_mbd_energy_s(...)` calls the ScaLAPACK version of it.


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
(calling from top-level `get_mbd_energy` instance and in ScaLAPACK framework create output files)

. eigenenergies: formatted text file `mbd_eigenvalues(_reciprocal).out` sorted by k points (if applicable)

. eigenmodes: unformatted binary file `mbd_eigenmodes(_kptX).out` sorted by k points (if applicable)
    - `mbd_eigenmodes.out`: size 1st dim, size 2nd dim, elements of 1st eigenmode, elements of 2nd eigenmode, ...
    - `mbd_eigenmodes_kptX.out`: size 1st dim, size 2nd dim, kpointX(1), kpointX(2), kpointX(3), elements of 1st eigenmode on kpoint X, elements of 2nd eigenmode on kpoint X, ...


