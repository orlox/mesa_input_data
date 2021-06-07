# 2016 binary models

**MESA version: 15140**
This is an update to the files to run them on version 15140.
The section on opacity tables and nets (2 and 3) are unchanged, for those the files
on the upper directory should be used.

To easily download these files, simply use
```
svn export https://github.com/orlox/mesa_input_data/trunk/2016_binary_models 2016_binary_models
```

These models are intended to reproduce the single star LMC models from Brott et al. (2011)
(A&A, 530, A115), and also include the contribution from binary systems. Just in case,
input data is provided as well for the GAL and SMC models of Brott et al. (2011), although
these have not been well tested.

A description of each directory is as follows.

- [1_composition](https://github.com/orlox/mesa_input_data/tree/master/2016_binary_models#1_composition): Generates input composition files for MESA. Unless you plan
on changing the nuclear network, you need only extract the input files
xa_GAL.data, xa_LMC.data, and xa_SMC.data.
- [2_opacity_tables](https://github.com/orlox/mesa_input_data/tree/master/2016_binary_models#2_opacity_tables): Custom opacity tables produced from the Opal website.
- [3_nets](https://github.com/orlox/mesa_input_data/tree/master/2016_binary_models#3_nets): Slightly modified basic and co_burn nuclear networks.
- [4_ZAMS_models](https://github.com/orlox/mesa_input_data/tree/master/2016_binary_models#4_zams_models): Produces sets of ZAMS models for each composition.
- [5_single](https://github.com/orlox/mesa_input_data/tree/master/2016_binary_models#5_single): Template work directory to model single stars.
- [6_binary](https://github.com/orlox/mesa_input_data/tree/master/2016_binary_models#6_binary): Template work directory to model binary stars.

Users interested simply in reproducing the single and binary models (not on modifying the
composition or nuclear networks used), can jump directly to the description on 5_single
to setup single stellar models.

Throughout this README, it is assumed that the user already has MESA properly installed,
if not, follow in detail the instructions given in the [MESA website](http://mesa.sourceforge.net/prereqs.html).

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
 LMC X,Y,Z,sum  0.73829020774503529       0.25687885050661802        4.8309709747724897E-003   1.0000000292264257     
 SMC X,Y,Z,sum  0.74599458477532132       0.25183115544024870        2.1742893275569371E-003   1.0000000295431271     
 GAL X,Y,Z,sum  0.72565584825315910       0.26515653393555483        9.1876465183556820E-003   1.0000000287070696     
 final LMC X,Y,Z,sum  0.73829020774503529       0.25687882128019224        4.8309709747724862E-003   1.0000000000000000     
 final SMC X,Y,Z,sum  0.74599458477532132       0.25183112589712175        2.1742893275569402E-003   1.0000000000000000     
 final GAL X,Y,Z,sum  0.72565584825315910       0.26515650522848522        9.1876465183556837E-003   1.0000000000000000
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

The tables included in this directory are already organized in a MESA-friendly way, and the
way to set them up is described in the 5_single section regarding single star models.

## 3_nets
As already mentioned, we use slightly modified basic and co_burn nets, these simply include the
isotopes fe56 and ca40 as a fillup element to account for the mass of other isotopes not included
in the network. Both of these are inert isotopes. 

## 4_ZAMS_models
This directory is to recreate the ZAMS models for each composition. You need to setup the opacity
tables and networks as described in the following section to use.

To use, modify the files
inlist_create_zams and inlist_zams_specification to use the desired composition (GAL, LMC or SMC),
cd into the directory, compile and run
```
./clean && ./mk
./rn
```
This will create the files with ZAMS models BROTT_GAL.data, BROTT_LMC.data and BROTT_SMC.data.
By default models are computed for the mass range 1-100 Msun, MESA can still start models outside
this mass range, but does so by adjusting the mass of the model at the beginning of the run. If
convenient, users can modify inlist_zams_specification to change the range of computed models.

## 5_single
And to finally compute some stellar models!

First, the custom opacity tables, nets, and ZAMS models need to be included into MESA.
Assuming $MESA_INPUT is the directory where these files have been downloaded, and
$MESA_DIR is the mesa directory to use, then the nets and ZAMS models can be easily imported
using
```
cp $MESA_INPUT/3_nets/*.net $MESA_DIR/data/net_data/nets/
cp $MESA_INPUT/4_ZAMS_models/*.data $MESA_DIR/data/star_data/zams_models
```
To include the opacity tables requires a bit more work, as these need to through the MESA
preprocessor first. This only parses the LMC opacity tables, if you want the LMC and GAL ones
you need to uncomment them on the rebuild all file.
```
cd $MESA_DIR/kap/preprocessor/
tar -Jxvf kap_input_data.tar.xz
cp $MESA_INPUT/2_opacity_tables/BROTT_* $MESA_DIR/kap/preprocessor/kap_input_data/opal/
cp -r $MESA_INPUT/2_opacity_tables/Type2_BROTT_* $MESA_DIR/kap/preprocessor/kap_input_data/opal/
cp $MESA_INPUT/2_opacity_tables/inlist* .
cp $MESA_INPUT/2_opacity_tables/rebuild_all .
./rebuild_all
./build_4_export
cd ..
./build_data_and_export
```
and now all neccesary data should be included into MESA. Note that we rewrite the rebuild_all script,
in case you have other custom opacities included into MESA, or have troubles getting this to work
in newer MESA versions, you might want to manually include the neccesary entries into rebuild_all.

And now to run a single star model, simply cd into 5_single, adjust inlist_extra to set
the initial mass, rotation velocity and inlist_project to specify
whether the model is GAL, LMC or SMC (you also need to copy locally the corresponding xa file,
by default xa_LMC.data is included). Then compile and run
```
./clean && ./mk
./rn
```
This work directory includes a run_star_extras that does multiple things:
- Define stellar winds as in Brott et al. (2011)
- Store profiles at specific points through H, He and C burning
- Include a timestep control in terms of absolute changes in central H
- Relax the model to the main sequence

The relaxation process is done to ignore loops at the beginning of the main sequence, caused due to the star being
off CN-equilibrium. The main steps of this relaxation are as follows:
- Relax the rotational velocity to the given initial value
- Model the star, without mass loss, until |log L - log L_nuc| < 0.005. At this point the star will essentially be
in CN equilibrium
- By this point there will be an offset in the rotational velocity from the desired one. To compensate for this,
re-relax the rotational velocity to the one we want.
- Proceed with normal evolution, turning on mass loss.

The output from the early, pre-relaxtion phase, should be ignored afterwards.

## 6_binary
Run a binary
