# 2016 binary models

MESA version: 8845

These models are intended to reproduce the single star LMC models from Brott etal. (2011)
(A&A, 530, A115), and also include the contribution from binary systems. Just in case,
input data is provided as well for the GAL and SMC models of Brott (2011).

A description of each folder is as follows.

- 1_composition: Generates input composition files for MESA. Unless you plan
on changing the nuclear network, you need only extract the input files
xa_GAL.data, xa_LMC.data, and xa_SMC.data.
- 2_opacity_tables: Custom opacity tables produced from the Opal website.
- 3_nets: Slightly modified basic and co_burn nuclear networks.
- 4_ZAMS_models: Produces sets of ZAMS models for each composition.
- 5_single: Template work directory to model single stars.
- 6_binary: Template work directory to model binary stars.

## 1_composition
MESA work directory that uses the chem module to produce input composition
files for MESA with the custom GAL, LMC and SCM mixtures of Brott etal. (2011).
For our models, we use the simple basic and co_burn networks of MESA, with the inclusion
of iron and calcium. Iron is included because we scale winds in terms of iron abundance,
calcium is included as an isotope to pile up everything that is not directly included
in the nuclear reaction network.

To use, cd into the working direct, compile and run

```
./clean && ./mk
./rn
```

this will write down the three file xa_GAL.data, xa_LMC.data, and xa_SMC.data. 
