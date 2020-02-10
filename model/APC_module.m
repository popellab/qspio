% APC Module
%
% Models APC activation and transport
%
% Inputs: model        -- simbio model object with four compartments 
%         species_name -- name of Tcells [must be unique]
%         params       -- object containing model parameter Values, Units, and Notes:
%                         - k_APC_mat--rate of APC maturation
%                         - k_APC_mig--rate of APC migration
%                         - k_APC_dth--rate of APC death
%                         - k_mAPC_dth--rate of mAPC death
%                         - APC0_T--baseline APC density in tumour
%                         - APC0_LN--baseline APC density in LN
%                         - k_c--cytokine time constant 
%                         - c0--baseline cytokine concentration
%                         - c50--cytokine concentration for half-maximal APC maturation
%                         - DAMPS--concentration of cytokines released per dying cancer cell 
%        
% Outputs: model -- SimBiology model object with new APC module
%
% Created: Oct 29, 2018 (Richard Sové)
% Last Modified: Jan 06, 2019 (RJS)

function model = APC_module(model,params)

% Add Species
APC = addspecies(model.Compartment(3),'APC',0,'InitialAmountUnits','cell');
    set(APC,'Notes','Number of naive antigen presenting cells in the tumour');
APC = addspecies(model.Compartment(4),'APC',0,'InitialAmountUnits','cell');
    set(APC,'Notes','Number of naive antigen presenting cells in the lymph node');
mAPC = addspecies(model.Compartment(3),'mAPC',0,'InitialAmountUnits','cell');
    set(mAPC,'Notes','Number of mature antigen presenting cells in the tumour');
mAPC = addspecies(model.Compartment(4),'mAPC',0,'InitialAmountUnits','cell');
    set(mAPC,'Notes','Number of mature antigen presenting cells in the lymph node');
c = addspecies(model.Compartment(3),'c',0,'InitialAmountUnits','molarity');
    set(c,'Notes','Concentration of maturation cytokines in the tumour');

% Add Parameters 
k_APC_mat = addparameter(model,'k_APC_mat',params.k_APC_mat.Value,'ValueUnits',params.k_APC_mat.Units);
    set(k_APC_mat,'Notes','Maximum rate of APC maturation');
k_APC_mig = addparameter(model,'k_APC_mig',params.k_APC_mig.Value,'ValueUnits',params.k_APC_mig.Units);
    set(k_APC_mig,'Notes','Rate of APC migration');
k_APC_death = addparameter(model,'k_APC_death',params.k_APC_death.Value,'ValueUnits',params.k_APC_death.Units);
    set(k_APC_death,'Notes','Rate of APC death');
k_mAPC_death = addparameter(model,'k_mAPC_death',params.k_mAPC_death.Value,'ValueUnits',params.k_mAPC_death.Units);
    set(k_mAPC_death,'Notes','Rate of mAPC death');
APC0_T = addparameter(model,'APC0_T',params.APC0_T.Value,'ValueUnits',params.APC0_T.Units);
    set(APC0_T,'Notes','APC density in the tumour');
APC0_LN = addparameter(model,'APC0_LN',params.APC0_LN.Value,'ValueUnits',params.APC0_LN.Units);
    set(APC0_LN,'Notes','APC density in the LN');
k_c = addparameter(model,'k_c',params.k_c.Value,'ValueUnits',params.k_c.Units);
    set(k_c,'Notes','Cytokine rate constant');
c0 = addparameter(model,'c0',params.c0.Value,'ValueUnits',params.c0.Units);
    set(c0,'Notes','Baseline cytokine concentration');
c50 = addparameter(model,'c50',params.c50.Value,'ValueUnits',params.c50.Units);
    set(c50,'Notes','Cytokine concentration for half-maximal APC maturation'); 
DAMPs = addparameter(model,'DAMPs',params.DAMPs.Value,'ValueUnits',params.DAMPs.Units);
    set(DAMPs,'Notes',['Concentration of cytokines released per dying cancer cell ' params.DAMPs.Notes]);
n_sites = addparameter(model,'n_sites_APC',params.n_sites_APC.Value,'ValueUnits',params.n_sites_APC.Units);
    set(n_sites,'Notes',['Maxium number of T Cells an APC can interact with ' params.n_sites_APC.Notes]);

% Add Reactions
% APC Reacruitment/Death in the Tumour 
reaction = addreaction(model,'null -> V_T.APC');
    set(reaction,'ReactionRate','k_APC_death*(APC0_T*V_T-V_T.APC)');
    set(reaction,'Notes','APC recruitment/death in the tumour');
% APC Recruitment/Death in LN
reaction = addreaction(model,'null -> V_LN.APC');
    set(reaction,'ReactionRate','k_APC_death*(APC0_LN*V_LN-V_LN.APC)');
    set(reaction,'Notes','APC recruitment/death in LN');
% APC Maturation 
reaction = addreaction(model,'V_T.APC -> V_T.mAPC');
    set(reaction,'ReactionRate','k_APC_mat*c/(c+c50)*V_T.APC');
    set(reaction,'Notes','APC maturation in the tumour');
% APC Migration to LN 
reaction = addreaction(model,'V_T.mAPC -> V_LN.mAPC');
    set(reaction,'ReactionRate','k_APC_mig*V_T.mAPC');
    set(reaction,'Notes','APC migration to the lymph node');
% mAPC Death in T 
reaction = addreaction(model,'V_T.mAPC -> null');
    set(reaction,'ReactionRate','k_mAPC_death*V_T.mAPC');
    set(reaction,'Notes','mAPC death in the tumour');    
% mAPC Death in LN 
reaction = addreaction(model,'V_LN.mAPC -> null');
    set(reaction,'ReactionRate','k_mAPC_death*V_LN.mAPC');
    set(reaction,'Notes','mAPC death in the lymph node');
% Baseline Cytokine Secretion/Degradation
reaction = addreaction(model,'null -> V_T.c');
    set(reaction,'ReactionRate','k_c*(c0-V_T.c)');
    set(reaction,'Notes','Baseline cytokine secretion/degradation');
% Cytokine Release in Response to the Tumour 
reaction = addreaction(model,'null -> V_T.c');
    set(reaction,'ReactionRate','R_Tcell*DAMPs');
    set(reaction,'Notes','Cytokine release in response to the tumour');

% Update Hill Functions for APCs
addrule(model,'H_APC = n_sites_APC*V_LN.APC/(n_sites_APC*V_LN.APC+T_total_LN+cell)','repeatedAssignment');
addrule(model,'H_mAPC = n_sites_APC*V_LN.mAPC/(n_sites_APC*V_LN.mAPC+T_total_LN+cell)','repeatedAssignment');