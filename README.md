# CAK_1d

Make HD model of a CAK radiation-driven wind from a OB star using MPI-AMRVAC. Stellar wind is spherically symmetric and assumed to be isothermal. 

## Setup

After cloning go into the local directory and fetch AMRVAC makefile (assuming you have exported the AMRVAC_DIR variable in the .bashrc (linux) or .bash_profile (macos))
```
$ setup.pl -d=1
```

and do a make (if the AMRVAC source was not compiled yet, this will be done now as well)
```
$ make -j
```

## User options

Simulations can be run for a specific OB star using the .par file. The goal of the CAK simulations for our purposes is to create a relaxed, steady-state wind model to be used in various contexts (in this case the LDI). The current .par file has stellar parameters corresponding to a typical OSG in the Galaxy.

Given that the CAK simulations are rather trivial (i.e. runtime ~ 3 minutes on a modern laptop) there are no options included in the .par file to do restarts or resumes. Depending on the star/stellar parameters the runtime can be best adjusted based on the unit_time that is printed to the screen at the start of the simulation.

***NOTE:** Due to the nature of the CAK line-force we could perform without problem a parellel simulation in MPI-AMRVAC. However, as we use the .dat file as input file for our LDI simulations this is not possible. The reason being that the LDI line-force is not parallelized (and is also not easy to do so). Because for the restart mesh + block related information is stored in the .dat file, the LDI .par file has to have the **exact** same mesh settings as the CAK simulation (of course also the stellar parameters...).*

## User parameters

The meaning of the AMRVAC parameter lists and variables in the .par file can be looked up on the AMRVAC webpage. The .par files provided are 'basic' ones that work for the problem.

Additionally, a *star_list* is specified in the .par file containing variables relevant for our problem. The ones relevant for computations are converted in the code to unitless quantities.

- lstar = stellar luminosity (solar units) 
- mstar = stellar mass (solar units)
- rstar = stellar radius (solar units)
- twind = wind temperature (Kelvin)
- rhobound = boundary density (g/cm^3)
- alpha = CAK line-force parameter (no units)
- Qbar = Gayley's line-force parameter (no units)
- Qmax = OCR's cut-off parameter (no units)
- kappae = electron scattering opacity
- beta = beta velocity law power index (no units)
- ifrc = wind option

## Notice

Tested with AMRVAC version 2.3 (Fall 2019).

## Known issues

None
