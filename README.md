# CAK_1d

Make HD model of a CAK radiation-driven wind from a OB star using MPI-AMRVAC. The stellar wind is spherically symmetric and assumed to be isothermal. 

## Setup

After cloning go into the local directory and fetch AMRVAC makefile (assuming you have exported the AMRVAC_DIR variable in the .bashrc (linux) or .bash_profile (macos))
```
$ setup.pl -d=1
```

and do a make (if the AMRVAC source was not compiled yet, this will be done now as well)
```
$ make -j
```

## Additional user setup

In order to run this problem without problems we have to make some modifications to the AMRVAC source files. To do this without harming/messing up the AMRVAC source files we copy the required source files to our local CAK_1d directory. The included *local.make* file ensures that AMRVAC sees our modified, local source file and compiles this one without using the original source file in the compilation.

The only source file required is:
```
mod_input_output.t
```
living in the **{path to your amrvac}/amrvac/src/amrvacio** directory. Copy it to your local CAK_1d directory and rename it (e.g. *my_mod_input_output.t* as in the *local.make* file).

In this source file search for the following codeblock:
```
if(mod(domain_nx^D/block_nx^D,2)/=0) &
  call mpistop("number level 1 blocks in D must be even")
```
and comment these two lines (really one line in Fortran). This allows us to run the problem on just 1 block and not multiple, i.e. a serial simulation. Because AMRVAC is parallel coded it expects that the mesh is divided in multiple blocks via domain decomposition. These blocks are then normally distributed over processors to do a parallel simulation. However, as mentioned below, for LDI purposes we are forced to run the code in serial.

## User options

Simulations can be run for a specific OB star using the .par file. The goal of the CAK simulations is to create a relaxed, steady-state wind model to be used as initial wind condition in various contexts (in this case the LDI). The current .par file has stellar parameters corresponding to a typical OSG in the Galaxy.

Given that the CAK simulations are rather trivial (i.e. runtime ~ 3 minutes on a modern laptop/desktop) there are no options included in the .par file to do restarts or resumes. Depending on the star/stellar parameters the runtime can be best adjusted based on the unit_time that is printed to the screen at the start of the simulation.

***NOTE:** Due to the nature of the CAK line-force we could perform without problem a parellel simulation in MPI-AMRVAC. However, as we use the .dat file as input file for our LDI simulations this is not possible. The reason being that the LDI line-force is not parallelized (and is also not easy to do so). Because for the restart mesh + block related information is stored in the .dat file, the LDI .par file has to have the **exact** same mesh settings as the CAK simulation (of course also the stellar parameters...).*

## User parameters

The meaning of the AMRVAC parameter lists and variables in the .par file can be looked up on the AMRVAC webpage. The .par files provided are 'basic' ones that work for the problem.

Additionally, a *star_list* is specified in the .par file containing variables relevant for our problem. These parameters are converted in the code to unitless quantities by a subroutine.

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
