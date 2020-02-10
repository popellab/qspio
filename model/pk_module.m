% Pharmacokinetic Module
%
% Models pharmacokinetics in four compartments
%
% Inputs: model        -- SimBiology model object with four compartments: C,P,T,LN
%         species_name -- name of drug [must be unique]
%         params       -- object containing model parameter
%                         - q_P--rate of diffusive transport C<->P
%                         - q_T--rate of diffusive transport C<->T
%                         - q_LN--rate of diffusive transport C<->LN
%                         - q_LD--rate of convective transport T->LN 
%                         - k_cl--clearence rate from central
%                         - gamma_C--volume fraction in C
%                         - gamma_P--volume fraction in P
%                         - gamma_T--volume fraction in T
%                         - gamma_LN--volume fraction in LN
%        
% Outputs: model -- SimBiology model object with new pk module
%
% Created: Oct 22, 2018 (Richard Sové)
% Last Modified: Dec 13, 2018 (RJS)

function model = pk_module(model,species_name,params)

% Get Compartments
comp_C = model.Compartment(1);
comp_P = model.Compartment(2);
comp_T = model.Compartment(3);
comp_LN = model.Compartment(4);

% Add Species
% Central
A_C = addspecies(comp_C,'A',0.0,'InitialAmountUnits','molarity'); 
    set(A_C,'Notes',['Concentration of ' species_name ' in central compartment']);
% Peripheral
A_P = addspecies(comp_P,'A',0.0,'InitialAmountUnits','molarity'); 
    set(A_P,'Notes',['Concentration of ' species_name ' in peripheral compartment']);
% Tumour
A_T = addspecies(comp_T,'A',0.0,'InitialAmountUnits','molarity'); 
    set(A_T,'Notes',['Concentration of ' species_name ' in tumour compartment']);
% Lymph Node
A_LN = addspecies(comp_LN,'A',0.0,'InitialAmountUnits','molarity'); 
    set(A_LN,'Notes',['Concentration of ' species_name ' in lymph node compartment']);

% Add Parameters
q_P = addparameter(model,'q_P',params.q_P.Value,'ValueUnits',params.q_P.Units); 
q_T = addparameter(model,'q_T',params.q_T.Value,'ValueUnits',params.q_T.Units); 
q_LN = addparameter(model,'q_LN',params.q_LN.Value,'ValueUnits',params.q_LN.Units); 
q_LD = addparameter(model,'q_LD',params.q_LD.Value,'ValueUnits',params.q_LD.Units); 
k_cl = addparameter(model,'k_cl',params.k_cl.Value,'ValueUnits',params.k_cl.Units);
gamma_C = addparameter(model,'gamma_C',params.gamma_C.Value,'ValueUnits',params.gamma_C.Units);
gamma_P = addparameter(model,'gamma_P',params.gamma_P.Value,'ValueUnits',params.gamma_P.Units);
gamma_T = addparameter(model,'gamma_T',params.gamma_T.Value,'ValueUnits',params.gamma_T.Units);
gamma_LN = addparameter(model,'gamma_LN',params.gamma_LN.Value,'ValueUnits',params.gamma_LN.Units);

% Add Reactions
% Diffusive Transport: Central to Peripheral
reaction = addreaction(model,'V_C.A <-> V_P.A'); % q_P = perm_CP_Ab*P_ratio_Nivo*S_CP*Peri, A = Nivo/f_vol
    set(reaction,'ReactionRate','q_P*(V_C.A/gamma_C - V_P.A/gamma_P)'); 
    set(reaction,'Notes',[species_name ' diffusive transport to peripheral compartment']);
% Diffusive Transport: Central to Tumour
reaction = addreaction(model,'V_C.A <-> V_T.A'); % q_T = perm_CT_Ab*S_CT*Tum
    set(reaction,'ReactionRate','q_T*(V_C.A/gamma_C - V_T.A/gamma_T)');
    set(reaction,'Notes',[species_name ' diffusive transport to tumor compartment']);
% Diffusive Transport: Central to Lymph Node
reaction = addreaction(model,'V_C.A <-> V_LN.A');
    set(reaction,'ReactionRate','q_LN*(V_C.A/gamma_C - V_LN.A/gamma_LN)');
    set(reaction,'Notes',[species_name ' diffusive transport to lymph node compartment']);
% Convective Transport: Tumour to Lymph Node
reaction = addreaction(model,'V_T.A -> V_LN.A');
    set(reaction,'ReactionRate','q_LD*V_T.A/gamma_T');
    set(reaction,'Notes',[species_name ' convective transport from tumor to lymph node']);
% Convective Transport: Lymph Node to Central
reaction = addreaction(model,'V_LN.A -> V_C.A');
    set(reaction,'ReactionRate','q_LD*V_LN.A/gamma_LN');
    set(reaction,'Notes',[species_name ' convective transport from lymph node to central']);
% Clearence from Central
reaction = addreaction(model,'V_C.A -> null');
    set(reaction,'ReactionRate','k_cl*V_C.A');
    set(reaction,'Notes','Drug clearance from central compartment');

% Rename Objects with 'species_name'
rename(A_C,species_name);
rename(A_P,species_name);
rename(A_T,species_name);
rename(A_LN,species_name);
rename(q_P,['q_P_' species_name]); 
rename(q_T,['q_T_' species_name]); 
rename(q_LN,['q_LN_' species_name]); 
rename(q_LD,['q_LD_' species_name]); 
rename(k_cl,['k_cl_' species_name]); 
rename(gamma_C,['gamma_C_' species_name]); 
rename(gamma_P,['gamma_P_' species_name]); 
rename(gamma_T,['gamma_T_' species_name]); 
rename(gamma_LN,['gamma_LN_' species_name]); 