# 2016 Double BH

**MESA version: 8115**
but you are welcome to try/adjust these files for future versions.

To easily download these files, simply use
```
svn export https://github.com/orlox/mesa_input_data/trunk/2016_double_bh 2016_double_bh
```

Input data to reproduce results in Marchant et al. (2016), A&A 588, A50. For convenience,
we also provide input files that work with MESA version 8845. A test_suite case implementing
this channel of evolution should remain updated in the code, and can be found in
$MESA_DIR/binary/test_suite/double_bh

![mass_period](mass_period.png)

## Z04, Z10, Z20, Z50
Each of these folders contains a work directory to model a binary system at different metallicities,
with Z04 corresponding to models with metallicity Zsun/4 (using Zsun=0.017 from Grevesse et al. 1996).
Masses and initial orbital periods are specified in the *inlist_extra* file,
```
&binary_controls
m1 = 70.0d0
m2 = 65.0d0
initial_period_in_days = 1.10d0
/ ! end of binary_controls namelist

```
For all templates except the one for Zsun/4, input parameters are the same.
The Zsun/4 models are different because due to the high metallicity, it was required to use
the reduction of superadiabaticity, aka MLT++ (see section 7.2 of Paxton et al. 2013).
This is specified in the files *inlist1* and *inlist2*
```
(...)
    okay_to_reduce_gradT_excess = .true.
    gradT_excess_age_fraction = 0.95d0
    gradT_excess_max_change = 0.001d0
(...)
```

The implementation of the custom wind scheme is done in *src/run_star_extras.f*. And the
condition for termination at helium depletion is set in *src/run_binary_extras.f*. The
reason to not specify this termination condition in the inlists, is because after one
component depletes helium, it is assumed to transform into a BH of the same mass, and
evolution follows assuming it is a point mass.
