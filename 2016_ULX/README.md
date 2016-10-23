# 2016 Double BH

**MESA version: 8845**
but you are welcome to try/adjust these files for future versions.

To easily download these files, simply use
```
svn export https://github.com/orlox/mesa_input_data/trunk/2016_ULX 2016_ULX
```

Input data to reproduce results in Marchant et al. (2016), to be submitted.

Throughout this README, it is assumed that the user already has MESA properly installed,
if not, follow in detail the instructions given in the [MESA website](http://mesa.sourceforge.net/prereqs.html).

![ULX](ULX.png)

## template
This folder contains a template which was used to model systems at all metallicities.
Masses and initial orbital periods are specified in the *inlist_extra* file,
```
&binary_controls
m1 = 70.7945784384d0
m2 = 14.1589156877d0
initial_period_in_days = 1.100d0
/ ! end of binary_controls namelist

```
and the composition is specified in the inlist_extra_sj file,
```
&star_job
relax_initial_Z = .true.
new_Z = 0.000316227766d0
relax_initial_Y = .true.
new_Y = 0.248300832755d0
/ ! end of star_job namelist

&controls
Zbase = 0.000316227766d0
/ ! end of controls namelist


```
which corresponds to the ULX model with q=0.2 from the paper discussed in section 3.

To run any of these models, cd into the corresponding directory, adjust inlist_extra,
and inlist_extra_sj, compile and run
```
./clean && ./mk
./rn
```

## data
This folder contains compressed data files. This includes summary tables with the outcomes of evolution for
each of the ~120000 models, and properties of the binary models during phases of mass transfer from the
secondary to a BH.

## scripts
These are mainly plotting tools and also scripts used to compute certain results. No
support provided for these, their here only for my convenience.
