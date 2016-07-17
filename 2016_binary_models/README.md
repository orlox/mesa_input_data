# 2016 binary models

MESA version: 8845

These models are intended to reproduce the single star LMC models from Brott et al. (2011)
(A&A, 530, A115), and also include the contribution from binary systems. Just in case,
input data is provided as well for the GAL and SMC models of Brott et al(2011), although
it has not been well tested.

A description of each folder is as follows.

- 1_composition: Generates input composition files for MESA. Unless you plan
on changing the nuclear network, you need only extract the input files
xa_GAL.data, xa_LMC.data, and xa_SMC.data.
- 2_opacity_tables: Custom opacity tables produced from the Opal website.
- 3_nets: Slightly modified basic and co_burn nuclear networks.
- 4_ZAMS_models: Produces sets of ZAMS models for each composition.
- 5_single: Template work directory to model single stars.
- 6_binary: Template work directory to model binary stars.

Users interested simply in reproducing the single and binary models (not on modifying the
composition or nuclear networks used), can jump directly to the description on 5_single
to setup single stellar models.

Throughout this README, it is assumed that the user already has MESA properly installed,
if not, follow in detail the instructions given in the [MESA website](http://mesa.sourceforge.net/prereqs.html)

## 1_composition
MESA work directory that uses the chem module to produce input composition
files for MESA with the custom GAL, LMC and SCM mixtures of Brott et al. (2011).
For our models, we use the simple basic and co_burn networks of MESA, with the inclusion
of iron and calcium. Iron is included because we scale winds in terms of iron abundance,
calcium is included as an isotope to pile up everything that is not directly included
in the nuclear reaction network.

The code contained in 1_composition/src/test_chem.f derives mass fractions in terms of abundances in (12+log10(n_i/n_h)) format,
assuming helium goes linearly from its primordial value Y=0.2477 (Peimbert et al. 2007)
at Z=0, to Y=0.28 at the solar metallicity fro Grevesse et al. (1996), Z=0.017. To do this, three equations must be satisfied:
```
1 = X+Y+Z
Y = 0.2477 + (0.28-0.2477)Z/0.017 = 0.2477 + 1.9 Z
Z = X * alpha = 0.7523 / (alpha^(-1)+2.9),
```
where alpha is given by
```
alpha = sum_i m_i/m_H * 10^(abundance_i-12)
```
and m_i, m_H are the atomic masses of isotope i and hydrogen. These equations can
be easily rearranged to give X in terms of alpha
```
X = 0.7523 / (1 + 2.9*alpha),
```
and after computing alpha, X, Y and Z, individual abundances are given by
```
X_i = X *  m_i/m_H * 10^(abundance_i-12)
```

To use, cd into the working direct, compile and run

```
./clean && ./mk
./rn
```

this will write down the three file xa_GAL.data, xa_LMC.data, and xa_SMC.data. 
