# 2017 IIb

**MESA version: 9575**
but you are welcome to try/adjust these files for future versions.

To easily download these files, simply use
```
svn export https://github.com/orlox/mesa_input_data/trunk/2017_IIb 2017_IIb
```

Input data to reproduce results in Sravan et al. (2017), submitted.

Throughout this README, it is assumed that the user already has MESA properly installed,
if not, follow in detail the instructions given in the [MESA website](http://mesa.sourceforge.net/prereqs.html).

![IIb](IIb.png)

## binary_template
This folder contains a template which was used to model binary systems.
Masses and initial orbital periods are specified in the *inlist_extra* file,
```
&binary_controls
m1 = 20d0
m2 = 15d0
initial_period_in_days = 100d0
/ ! end of binary_controls namelist

```
To adjust mass transfer efficiency, set mass_transfer_beta equal to 1-efficiency,
i.e., for an efficiency of 10%, set mass_transfer_beta=0.9d0 in inlist_project.

To run any of these models, cd into the corresponding directory, adjust inlist_extra,
compile and run
```
./clean && ./mk
./rn
```

## single_template
This folder contains a template which was used to model single stars.
Masses are specified in the *inlist_extra* file,
```
&controls
initial_mass = 20d0
/ ! end of controls namelist

```
To run any of these models, cd into the corresponding directory, adjust inlist_extra,
compile and run
```
./clean && ./mk
./rn
```
