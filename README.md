# **Alcott et al., Nat Geo 2024**

This five box ocean-atmosphere model is originally based on the model presented by Slomp and Van Cappellen 
(2007; Biogeosciences), which was later edited to consider additional aspects of the Earth system
by Tsandev et al. (2008; GBC) Tsandev and Slomp (2009; EPSL) and then tranlsated to consider geological timescales by
Alcott et al. (2019; Science). 

# **Contents:**

### Model Code:

Sensitivity_run.m // Model monte carlo simulation file // Standard sensitivity test includes 1000 model runs.

Alcott_et_al_2024_NatGeo_front.m // Model code that outlines starting parameters as well as ranges to be considered
in senstivity

Alcott_et_al_2024_NatGeo.m // Model equations file // code will not run directly // contains differential equations 
for reservoirs and flux calculations

CIplots.m // called after model has completed in order to plot results

### Data to validate model in plot:

AtmosO2proxy.mat // Atmospheric O2

CO2proxy.mat // Atmospheric CO2

Havig13C.mat // d13C data

ReinhardPdata.mat // Reinhard et al. 2017. P concentrations

### Data needed for model to run:

C_HW_2006.mat // carbonate fraction in the crust

C_Hayesandwaldbauer2006.csv // degassing rates



