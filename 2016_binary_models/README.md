# 2016 binary models

**MESA version: 8845**

These models are intended to reproduce the single star LMC models from Brott et al. (2011)
(A&A, 530, A115), and also include the contribution from binary systems. Just in case,
input data is provided as well for the GAL and SMC models of Brott et al(2011), although
these have not been well tested.

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
files for MESA with the custom GAL, LMC and SMC mixtures of Brott et al. (2011).
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
Terminal output reports the values for X,Y and Z for each mixture, where the 
"final" values are after readjusting Y as 1-X-Z so that abundances add up to 1
within machine precision,
```
LMC X,Y,Z,sum  0.73826111690125806       0.25689791002408885        4.8410022998830000E-003   1.0000000292252300     
SMC X,Y,Z,sum  0.74598059315035770       0.25184032236653253        2.1791140256615994E-003   1.0000000295425517     
GAL X,Y,Z,sum  0.72559344897166456       0.26519741622157317        9.2091635112669136E-003   1.0000000287045048     
final LMC X,Y,Z,sum  0.73826111690125806       0.25689788079885895        4.8410022998829991E-003   1.0000000000000000     
final SMC X,Y,Z,sum  0.74598059315035770       0.25184029282398068        2.1791140256616011E-003   1.0000000000000000     
final GAL X,Y,Z,sum  0.72559344897166456       0.26519738751706851        9.2091635112669119E-003  0.99999999999999989
```
Note that these differ slightly from the values of Brott et al. 2011. This is due to two reasons.
First, we include all isotopes from the Asplund et al. (2005) solar values when computing Z, while
the Brott models only consider those present in the nuclear network they use. Second, when computing
mass abundances from abundances by number, they use atomic numbers instead of atomic masses.

## 2_opacity_tables
Brott et al. (2011) use opacity tables for solar scaled composition. As the custom mixtures
considered are not simply solar scaled, they compensate for this by using an effective metallicity
for the computation of opacities
```
Zeff=Zsun * Fe / Fe_sun
```
Instead of using this approach, we create custom opal tables. To do this, we use the OPAL website,
which allows the generation of [Type 1 tables](http://opalopacity.llnl.gov/type1inp.html) (i.e. fixed
metal fractions and X,Z variable) and
[Type 2 tables](http://opalopacity.llnl.gov/type2inp.html) (i.e. X,Z and C, O abundances variable).
The input required for the opal website are relative abundances by number of metals, and these values are
provided in the xa_GAL.data, xa_LMC.data, and xa_SMC.data files in case the user wants to regenerate the tables,
or create new ones for different mixtures.
