
# 1D CAK radiation-driven wind

*Update 2023:* A base version of this user code is now included in MPI-AMRVAC as test problem ([see here](https://github.com/amrvac/amrvac/tree/master/tests/hd/CAKwind_spherical_1D)). The entire line-force computation is now also part of the source code. For the latter see the documentation on the webpage.

Make HD model of a CAK wind from a OB star using MPI-AMRVAC. The stellar wind is spherically symmetric and assumed to be isothermal. The difference with the MPI-AMRVAC test problem is that here the code has been extended to follow the model of [Poniatowski et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022A%26A...667A.113P/abstract). This accounts for spatially varying line-force parameters computed from local hydrodynamic conditions as well as an update to the lower boundary condition.

## Setup

After cloning go into the local directory and fetch AMRVAC makefile (assuming you have exported the `AMRVAC_DIR` variable in the `.bashrc` on Linux or `.bash_profile` on macOS)
```
$ setup.pl -d=1 -v=1
```
and do a make (if the AMRVAC source was not compiled yet, this will be done now as well)
```
$ make -j
```

## How-to

### Performing a simulation

Simulations can be run for a specific OB star using the .par file. The goal of the CAK simulations is to create a relaxed, steady-state wind model to be used as initial wind condition in various contexts. The current .par file has stellar parameters corresponding to a typical O-supergiant (Zeta Puppis) in the Galaxy.

Given that the CAK simulations are rather fast there are no options included in the .par file to do restarts or resumes. Depending on the stellar parameters the runtime can be best adjusted based on the `unit_time` variable that is printed to the screen at the start of the simulation.

## Additional user parameters

The meaning of the AMRVAC parameter lists and variables in the .par file can be looked up on the AMRVAC webpage. The .par files provided are 'basic' ones that work for the problem.

Additionally, a `star_list` is specified in the .par file containing variables specific for our problem. The ones relevant for computations are converted in the code to unitless quantities within the `make_dimless_vars` routine.

| Parameter| Explanation                                                       |
|----------|-------------------------------------------------------------------|
| mstar_sol    | stellar mass (solar units)                                    |
| rstar_sol    | stellar radius (solar units)                                  |
| twind_cgs    | wind temperature (Kelvin)                                     |
| rhosurf_cgs  | boundary density (g/cm^3)                                     |
| cak_alpha    | CAK line-force parameter (no units)                           |
| gayley_Qbar  | Gayley's line-force parameter (no units)                      |
| gayley_Qmax  | OCR's cut-off parameter (no units)                            |
| beta         | beta velocity law power index (no units)                      |
| ifrc         | wind option                                                   |
| tstat        | start time (in units of unit time) for invoking average quantity computation |
|use_lte_table | varying line-force parameters after Poniatowksi+ (2022)       |

## Notice

Tested with AMRVAC version 3.3, 3.2, 3.1, 3.0, and 2.2.

## Known issues

None
