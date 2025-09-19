%% Alcott et al., Nat Geo 2024
% Model based on Slomp and VanCappellen, 2007, Biogeosciences; Tsandev et
% al., 2008, GBC; Tsandev and Slomp, 2009, EPSL. Alcott et al., 2019
% Model front (runs a single model run)
function run = Alcott_et_al_2024_NatGeo_front(S)

%%% Set up Global parameters
global stepnumber
global pars
global workingstate
global when
global starting
global per
global present
global state
global sensparams

%% Options
when.start = -4.2e9 ; %Start time
when.end = 0; %End time
stepnumber = 1 ;

%%% Run parameters
per.sig_SCAV = 0; %Scavenging - percentage of upwelled that is scavenged
per.POP_deep_feedback = 0.1;%Percentage of POP burial that is redox dependent
per.CaP_deep_feedback = 0.1;%Percentage of Pauth burial that is redox dependent
pars.CPoxic = 106 ;

load('C_HW_2006.mat')
per.C_HW_2006_time = 1e9 * CHayesandwaldbauer2006(:,1);
per.C_HW_2006_data = CHayesandwaldbauer2006(:,2);

per.P =  1 ; %Riverine input of P relative to present
per.generic_red = 0 ; %Generic reduced gas input - 0 at present day

%% Initial Values

%% Hydrological Cycle
%%%%% Values from Slomp and Van Cappellen, 2007
starting.Water_P = 36e12 ;
starting.Water_D = 3600e12 ;       
starting.Water_S = 49830e12 ;
starting.Water_DP = 1.3e18 ;
%Proximal Coastal Zone
pars.y(1) = starting.Water_P; %m^3
%Distal Coastal Zone
pars.y(2) = starting.Water_D; %m^3
%Surface Layer Open Ocean
pars.y(3) = starting.Water_S; %m^3
%Deep Ocean
pars.y(4) = starting.Water_DP; %m^3

%% Marine Carbon Cycle
starting.POC_P = 4.5e12 ;
starting.POC_D = 243e12 ;
starting.POC_S = 3816e12 ;
starting.POC_DP = 5.6e16 ;
%Proximal Coastal Zone
pars.y(5) = starting.POC_P ; %mol
%Distal Coastal Zone
pars.y(6) = starting.POC_D ; %mol
%Surface Layer Open Ocean
pars.y(7) = starting.POC_S ; %mol
%Deep Ocean
pars.y(8) = starting.POC_DP ; %mol

%% Oxygen Cycle
starting.O2_P = 4.5e12 ;
starting.O2_D = 243e12 ;
starting.O2_S = 1.6145e16 ;
starting.O2_DP = 2.21e17;
starting.O2_A = 3.7e19 ;
%Proximal Coastal Zone
pars.y(9) = starting.O2_P ; %mol
%Distal Coastal Zone
pars.y(10) = starting.O2_D ; %mol
%Surface Layer Open Ocean
pars.y(11) = starting.O2_S ; %mol 
%Deep Ocean
pars.y(12) = starting.O2_DP ; %mol 
% Atmospheric O2 ;
pars.y(21) = starting.O2_A ; %mol

%% Phosphorous Cycle
starting.SRP_P = 9.7e9 ;
starting.OP_P = 4.3e10 ;
starting.SRP_D = 5e12 ;
starting.OP_D = 2.3e12 ;
starting.SRP_S = 47e12 ;
starting.OP_S = 36e12 ;
starting.SRP_DP = 2790e12 ;
starting.OP_DP = 530e12 ;
%SRP proximal coastal zone
pars.y(13) = starting.SRP_P ; %mol
%POP proximal coastal zone
pars.y(14) = starting.OP_P ; %mol
%SRP distal coastal zone
pars.y(15) = starting.SRP_D ; %mol
%POP distal coastal zone
pars.y(16) = starting.OP_D ; %mol
%SRP surface open ocean
pars.y(17) = starting.SRP_S ; %mol
%POP surface open ocean
pars.y(18) = starting.OP_S ; %mol
%SRP deep ocean
pars.y(19) = starting.SRP_DP ; %mol
%POP deep ocean
pars.y(20) = starting.OP_DP ; %mol

%% COPSE Carbon
starting.A_0 = 3.193e18 ; 
starting.Aiso_0 = 0 ; 
starting.GAST_0 = 288 ;

pars.y(22) = starting.A_0 ;
pars.y(23) = starting.Aiso_0 ;
%% Present day values to use for normalization
present.SRP_P = 9.7e9 ;
present.SRP_D = 5e12 ;
present.SRP_S = 47e12 ;
present.SRP_DP = 2790e12 ;
present.OP_P = 4.3e10 ;
present.OP_D = 2.3e12 ;
present.OP_S = 36e12 ;
present.OP_DP = 530e12 ;
present.POC_P = 4.5e12 ;
present.POC_D = 243e12 ;
present.POC_S = 3816e12 ;
present.POC_DP = 5.6e16 ; 
present.O2_D = 243e12 ;
present.O2_DP = 2.21e17 ;
present.O2_A = 3.7e19 ;
present.O2_S =1.6145e16 ;
present.O2_P = 4.5e12 ;

%% Starting concentrations 
% Calculate concentrations of SRP and OP for flux calculations
starting.SRP_Pconc = starting.SRP_P / starting.Water_P ;
starting.OP_Pconc = starting.OP_P / starting.Water_P ;
starting.SRP_Dconc = starting.SRP_D / starting.Water_D ;
starting.OP_Dconc = starting.OP_D / starting.Water_D ;
starting.SRP_Sconc = starting.SRP_S / starting.Water_S ;
starting.OP_Sconc = starting.OP_S / starting.Water_S ;
starting.SRP_DPconc = starting.SRP_DP / starting.Water_DP ;
starting.OP_DPconc = starting.OP_DP / starting.Water_DP ;
starting.O2_Sconc = (starting.O2_S / starting.Water_S) ;
starting.O2_DPconc = (starting.O2_DP/starting.Water_DP) ;
        
%% Calculate constants for circulation

% Water fluxes from Slomp and VC, 2007.
Water_S_DP_0 = 3780e12 + 378e12 ;
Deep_Water_S_0 = 3780e12 ;
Deep_Water_D_0 = 378e12 ;
Water_P_D_0 = 37e12 ;
Water_D_S_0 = 415e12 ;

% Constants for circulation of water masses.
pars.kWF4 = 4158e12 / starting.Water_S; 
pars.kWF5 = 3780e12 / starting.Water_DP;
pars.kWF6 = 378e12 / starting.Water_DP;


%% Calculations for Carbon Cycle

% Primary production in Water_P
pars.Redfield_CP = 106 ; %C:P ratio redfield
pars.Redfield_CO2 = 106/138; %C:O2 ratio redfield

starting.Prox_Prod_Photo = 3.975e13 ; 
pars.kPhotoprox = starting.Prox_Prod_Photo / pars.Redfield_CP  ;

% POC export from Water_P to Water_D
OP_P_D_0 = starting.OP_Pconc * Water_P_D_0 ;
XP_P_D_0 = OP_P_D_0 * pars.Redfield_CP ;

% Proximal sediment POC burial
starting.POC_P_Burial_0 = 2.3e12 ;
pars.Prox_C_Bur = starting.POC_P_Burial_0 / starting.Prox_Prod_Photo ;

% POC mineralisation in Prox_Water
POC_P_Min_0 = starting.Prox_Prod_Photo - starting.POC_P_Burial_0 - XP_P_D_0 ;
pars.kminprox = POC_P_Min_0 ;

% Primary Production in Dist_Water
starting.PP_D_0 = 5.6e14 ;
pars.kPhotodist = starting.PP_D_0 / pars.Redfield_CP  ;

% POC Export from Distal to Surface
Dist_OP_S_0 = starting.OP_Dconc * Water_D_S_0 ; 
XP_D_S_0 = Dist_OP_S_0 * pars.Redfield_CP ;

% Distal sediment POC burial
starting.POC_D_Burial_0 = 1.7e12 ;
pars.Dist_C_Bur = starting.POC_D_Burial_0 / ( XP_P_D_0 + starting.PP_D_0 ) ;

% POC mineralisation in Dist_Water
POC_Min_D_0 = XP_P_D_0 - XP_D_S_0 - starting.POC_D_Burial_0 + starting.PP_D_0 ;
pars.kmindist = POC_Min_D_0 ;

% Primary Production in Surf_Water
PP_S_0 = 3.8688e15 ;
pars.kPhotosurf = PP_S_0 /  pars.Redfield_CP ;

% monod constant for oxic respiration 
pars.KmO2 = 0.0001 ; 

% O2 Downwelling
O2_S_DP_0 = Water_S_DP_0 * starting.O2_Sconc ; 

% O2 coastal upwelling
O2_DP_D_0 =  Deep_Water_D_0 * starting.O2_DPconc; 

% O2 oceanic upwelling
O2_DP_S_0 = Deep_Water_S_0 * starting.O2_DPconc; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Reduced Fe flux %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pars.FeO = 0 ;
% respiration
Respiration_O2_0 = O2_S_DP_0 - O2_DP_D_0 - O2_DP_S_0 - pars.FeO ; 

% Monod relationship for respiration
Mon_O2_deep_0 = starting.O2_DPconc / ( pars.KmO2 + starting.O2_DPconc ) ; 

% respiration constant considering both steady state and monod relationship
Respiration_O21_0 = Respiration_O2_0 / Mon_O2_deep_0 ;

% Organic carbon respiration. CF12 is the name used in Slomp and VC, 2007
pars.kCF12 = Respiration_O21_0 * pars.Redfield_CO2 ;

% Present day respiration
POC_DP_Resp_0 = pars.kCF12 ;

%%% Present day value of organic carbon burial in deep and distal boxes
starting.POC_DP_Burial = 1e12 ;
starting.POC_D_Burial = 1.7e12 ;

% Required export for steady state.
XP_S_DP_0 = starting.POC_DP_Burial + POC_DP_Resp_0 ;

% Constant for surface - deep export
pars.Surf_Deep_XP = XP_S_DP_0 / ( XP_D_S_0 + PP_S_0 ) ;

% Carbon mineralisation in Surface ocean
POC_Min_S_0 = XP_D_S_0 - XP_S_DP_0 + PP_S_0 ;
pars.kminsurf = POC_Min_S_0 ;

%Deep Sea Sediment POC burial for steady state
starting.POC_DP_Burial_0 = XP_S_DP_0 - POC_DP_Resp_0 ; 

% Phosphorus Primary Production in Surface ocean
P_PP_S_0 = PP_S_0 / pars.Redfield_CP; 

% Slomp and VC, 2007
pars.kCF11 = ( 496.6 / ( 3600 + 28.0125 ) ); 
OP_S_DP_0 = pars.kCF11 * ( P_PP_S_0 + Dist_OP_S_0 ) ;

% Present day iron-bound P burial in deep ocean
% pars.kFeP_Deep = 6.75e9 ;
pars.kFeP_Deep = 1e12/106 ;

% Steady state of organic P cycle in deep ocean 
OP_DP_Min_0 = OP_S_DP_0 - pars.kFeP_Deep ;

% Remineralisation of Org-P in the deep ocean
pars.kPrel_deep = OP_DP_Min_0 / starting.OP_DP ; 

% OrgC burial in deep ocean
pars.Deep_C_Bur = starting.POC_DP_Burial_0 / ( pars.kPrel_deep * starting.POC_DP ) ; 


%% Phosphorus cycle constants

%% Proximal Coastal

% Present day fanoxic value used in proximal and distal boxes (Watson et al., 2017)
starting.fanoxic = 0.0025 ; 

% Riverine P input.
pars.River_SRP_0 = 0.09e12 ; %Slomp and VC, 2007; Berner and Rao,1994

% Primary Production in Water_P
P_PP_P_0 = starting.Prox_Prod_Photo / pars.Redfield_CP ; 

% SRP export from Water_P to Water_D
SRP_P_D_0 = starting.SRP_Pconc * Water_P_D_0 ; 

% POP export from Water_P to Water_D
OP_P_D_0 = starting.OP_Pconc * Water_P_D_0 ; 

% Proximal FeP burial
pars.kFePprox = 0.925 ;
P_FeP_P_0 = pars.kFePprox * starting.SRP_P ; 

% Proximal sediment POP burial
pars.korgP_prox = 1 ;
OP_P_Burial_0 = pars.Prox_C_Bur * starting.Prox_Prod_Photo * ( ( ( 1-starting.fanoxic ) / pars.CPoxic ) + ( starting.fanoxic / 1100 ) ) ;

% POP mineralisation in Proximal 
OP_P_Min_0 = P_PP_P_0 - OP_P_Burial_0 - OP_P_D_0 ; 
pars.kPrel_prox = OP_P_Min_0 ;

% Proximal CaP burial for steady state
P_AuthP_P_0 = pars.River_SRP_0 + OP_P_Min_0 - SRP_P_D_0 - P_FeP_P_0 - P_PP_P_0 ; 
pars.kCaP_P = P_AuthP_P_0 / ( OP_P_Min_0 * ( 1 - starting.fanoxic ) ) ; %NEW
%% Distal P

% SRP open ocean upwelling. Deep to Surface transport
Deep_Surf_P_0 = starting.SRP_DPconc * Deep_Water_S_0 ; 

% SRP coastal upwelling. Deep to Distal transport
Deep_Dist_P_0 =  starting.SRP_DPconc * Deep_Water_D_0 ;

% Primary Production in Distal box
P_PP_D_0 = starting.PP_D_0 / pars.Redfield_CP ; 

% SRP export. Distal to Surface
SRP_D_S_0 = starting.SRP_Dconc * Water_D_S_0 ;

% POP export from Distal to Surface waters
Dist_OP_S_0 = starting.OP_Dconc * Water_D_S_0 ;

% Distal sediment FeP burial
pars.kPF9 = 0.00135238 ; % As used in Slomp and VC, 2007. 
P_FeP_D_0 = starting.SRP_D * pars.kPF9 ;
pars.kFePDOADist = P_FeP_D_0 / ( starting.SRP_D * (1 - starting.fanoxic) ) ;

% Distal sediment POP burial
OP_D_Burial_0 = starting.POC_D_Burial_0 / pars.CPoxic ;
pars.kPOPDOADist = OP_D_Burial_0 / ( (starting.PP_D_0 + XP_P_D_0 ) * ( ( ( 1- starting.fanoxic ) / pars.CPoxic ) + ( starting.fanoxic / 1100 ) ) ) ;

% POP mineralisation in Distal
OP_D_Min_0 = P_PP_D_0 + OP_P_D_0 - OP_D_Burial_0 - Dist_OP_S_0 ; 
pars.kPrel_dist = OP_D_Min_0 ;

% Distal Sediment CaP burial
P_AuthP_D_0 = SRP_P_D_0 - P_PP_D_0 + OP_D_Min_0 - P_FeP_D_0 - SRP_D_S_0 + Deep_Dist_P_0; 
pars.kCaPDOADist = P_AuthP_D_0 / ( OP_D_Min_0 *(1-starting.fanoxic) ) ;


%% Surface Ocean P

% SRP downwelling from Surface to Deep ocean
SRP_S_DP_0 = starting.SRP_Sconc * Water_S_DP_0 ;

% POP export from Surface to Deep Ocean
OP_S_DP_0 = pars.kCF11 * ( P_PP_S_0 + Dist_OP_S_0 ) ; 

% POP mineralisation in Water_S for steady state
OP_S_Min_0 = Dist_OP_S_0 - OP_S_DP_0 + P_PP_S_0 ; 
pars.kPrel_surf = OP_S_Min_0 ;

%% Deep Ocean

%Deep FeP Burial
P_FeP_DP_0 = pars.kFeP_Deep ;

% POP mineralisation in Water_DP
OP_DP_Min_0 = OP_S_DP_0 - pars.kFeP_Deep ; 
pars.kPrel_deep = OP_DP_Min_0 / starting.OP_DP ;

% Deep POP Burial
OP_DP_Burial_0 = OP_S_DP_0 - OP_DP_Min_0 ; 
pars.kPOP_Bur_Deep = OP_DP_Burial_0 / (XP_S_DP_0 / pars.CPoxic) ;

% Deep sediment CaP burial for steady state
P_AuthP_DP_0 = SRP_S_DP_0 - P_FeP_DP_0 - Deep_Surf_P_0 - Deep_Dist_P_0 + OP_DP_Min_0 ;

% Flux PF34, as used in Slomp and VC, 2007.
pars.fPF34 = P_AuthP_DP_0 / OP_DP_Min_0 ; 


%% Weathering O2
% Total OrgC burial is equal to oxidative weathering in order for steady state Atmospheric O2.
pars.rgf = 0.4e12 ;
pars.Flocb_0 = 5e12 ;
pars.O2_A_Weathering = ( pars.Flocb_0 + starting.POC_P_Burial_0 + starting.POC_D_Burial_0 + starting.POC_DP_Burial_0 - pars.rgf)  ;

%% Inorg C
%Constants taken from COPSE, COPSE reloaded
%Set to represent present day to calculate present day rates

pars.fbiota = 1 ; %Realtive change in terrestrial biota. Keeping it at 1 means land plants were always there!
pars.kTgran = 0.0724 ;  %Granite weathering constant
pars.kTbas = 0.0608 ; %Basalt weathering constant
pars.kTsfw = 0.0608 ; %Sea floor weathering constant

%%%% org C cycle
    
    pars.Fmocb_0 = starting.POC_P_Burial_0 + starting.POC_D_Burial_0 + starting.POC_DP_Burial_0 ;
    pars.Fccdeg_0 = 25e12 ; 
    pars.Fcarbw_0 = 16e12 ;
    pars.Fsfw_0 = 1.75e12 ;
    pars.Fmccb_0 = pars.Fcarbw_0 + pars.Fccdeg_0 - pars.Fsfw_0 ;
    pars.Fsilw_0 = pars.Fmccb_0 - pars.Fcarbw_0 ;
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pars.Focdeg_0 = 0 ; % Organic carbon degassing

basfrac = 0.3 ;
    pars.k_granw = pars.Fsilw_0 * (1-basfrac) ;
    pars.k_basw = pars.Fsilw_0 * basfrac ;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%   Sensitivity analysis   %%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %%%%%%% parameter space to test
        sensparams.D = 0.25 + rand*(2-0.25); 

        sensparams.fbiota =  0.15 + rand*(1-0.15);

        sensparams.C = 0.25 + rand*(3-0.25) ; 
        
        sensparams.CP = pars.CPoxic + rand*(1100-pars.CPoxic) ; 
        
        sensparams.CP_1point4_oxic = 40 + rand*10 ; %%%% 40-50 
        
        sensparams.CP_1point4_anoxic = 800 + rand*100 ; %%%% 800-900

        sensparams.EXP = 0.1 + rand* (0.75-0.1) ;

        sensparams.EXP2 = sensparams.EXP + rand* (1-sensparams.EXP) ;

        sensparams.EXPtiming = 1.7 + rand*(3-1.7);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options = (odeset('maxstep', 1e8)) ;
tic
[rawoutput.T,rawoutput.Y] = ode15s(@Alcott_et_al_2024_NatGeo,[when.start when.end] , pars.y, options) ;

%%%% size of output 
pars.output_length = length(rawoutput.T) ;
%%%%%%%%%% model finished output to screen
fprintf('Integration finished \t') ; fprintf('Total steps: %d \t' , stepnumber ) ; fprintf('Output steps: %d \n' , pars.output_length ) 

%%%%%% assemble output state vectors
[sharedvals,trecords] = intersect(workingstate.time,rawoutput.T,'stable') ;
field_names = fieldnames(workingstate) ;
for numfields = 1:length(field_names)
    eval([' state.' char( field_names(numfields) ) ' = workingstate.' char( field_names(numfields) ) '(trecords) ; '])
end
       
%%%%%% save state
    run.state = state ;
    run.pars = pars ;
%     run.forcings = forcings ;
    
%%%%%% done message
fprintf('Done: ')
endtime = toc ;
fprintf('time (s): %d \n', endtime )
