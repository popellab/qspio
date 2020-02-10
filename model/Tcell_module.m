% Tcell Module
%
% Models T Cell transport and activation by APCs [Use before antigen corresponding module]
%
% Inputs: model        -- SimBiology model object with four compartments 
%         ID           -- T cell-antigen ID number [must be unique]
%         params       -- object containing model parameter Values, Units, and Notes:
%                         - Q_in--rate of naive T cells into the LN
%                         - q_out--rate of naive T cells out of the LN
%                         - k_act--rate of naive T cell activation
%                         - k_pro--rate of activated T cell proliferation
%                         - k_death--rate of mature T cell death
%                         - q_P_in--rate of T cell transport C->P
%                         - q_P_out--rate of T cell transport P->C
%                         - q_T_in--rate of T cell transport C->T
%                         - q_LN_out--rate of T cell transport LN->C
%                         - k_IL2_deg--rate of IL2 degradtion
%                         - k_IL2_cons--rate of IL2 consumption
%                         - k_IL2_sec--rate of IL2 secretion
%                         - IL2_50--half-maximal IL2 concentration for T cell activation
%                         - IL2_50_Treg--half-maximal IL2 concentration for Treg activation
%                         - N0--minimum number of activated T cell generations (no IL2)
%                         - N_IL2--maximum number of generations added due to IL2 
%                         - n_clones--number of T cell clones
%         cancer_types -- cell array of strings containing names of cancer types Teff kill
%                         <leave empty for Tregs>
%        
% Outputs: model -- SimBiology model object with new Tcell module
%
% Created: Oct 29, 2018 (Richard Sové)
% Last Modified: Sep 17, 2019 (RJS)

function model = Tcell_module(model,ID,params,cancer_types)

% Species Names
species_name = ['T' ID];
antigen = ['P' ID];

% Determine if Treg
isTreg = false;
if (strcmp(ID,'0'))
    isTreg = true;
end

% Add Species
% Naive T cells
nT_C  = addspecies(model.Compartment(1),'nT',params.nT_C.Value/params.div.Value,'InitialAmountUnits','cell');
    set(nT_C,'Notes',['Number of naive ' species_name ' cells in the central compartment']);
nT_P  = addspecies(model.Compartment(2),'nT',params.nT_P.Value/params.div.Value,'InitialAmountUnits','cell');
    set(nT_P,'Notes',['Number of naive ' species_name ' cells in the peripheral compartment']);
nT_LN = addspecies(model.Compartment(4),'nT',params.nT_LN.Value/params.div.Value,'InitialAmountUnits','cell');
    set(nT_LN,'Notes',['Number of naive ' species_name ' cells in the lymph node']);
% Activated T cells    
aT = addspecies(model.Compartment(4),'aT',0,'InitialAmountUnits','cell');
    set(aT,'Notes',['Number of activated ' species_name ' cells in the lymph node']);
% Mature T cells
if (isTreg)
    T_C = addspecies(model.Compartment(1),'T',params.nTreg_C.Value,'InitialAmountUnits','cell');
        set(T_C,'Notes',['Number of ' species_name ' cells in the central compartment']);
    T_P = addspecies(model.Compartment(2),'T',params.nTreg_P.Value,'InitialAmountUnits','cell');
        set(T_P,'Notes',['Number of ' species_name ' cells in the peripheral compartment']);
    T_T = addspecies(model.Compartment(3),'T',0,'InitialAmountUnits','cell');
        set(T_T,'Notes',['Number of ' species_name ' cells in the tumour compartment']);
    T_LN = addspecies(model.Compartment(4),'T',params.nTreg_LN.Value,'InitialAmountUnits','cell');
        set(T_LN,'Notes',['Number of ' species_name ' cells in the lymph node compartment']);
else
    T_C = addspecies(model.Compartment(1),'T',0,'InitialAmountUnits','cell');
        set(T_C,'Notes',['Number of ' species_name ' cells in the central compartment']);
    T_P = addspecies(model.Compartment(2),'T',0,'InitialAmountUnits','cell');
        set(T_P,'Notes',['Number of ' species_name ' cells in the peripheral compartment']);
    T_T = addspecies(model.Compartment(3),'T',0,'InitialAmountUnits','cell');
        set(T_T,'Notes',['Number of ' species_name ' cells in the tumour compartment']);
    T_LN = addspecies(model.Compartment(4),'T',0,'InitialAmountUnits','cell');
        set(T_LN,'Notes',['Number of ' species_name ' cells in the lymph node compartment']);
end
% Determine if first call    
first_call = true;
try % add IL2 if it does not exist yet
IL2 = addspecies(model.Compartment(4),'IL2',1e-18,'InitialAmountUnits','molarity');
    set(IL2,'Notes','Concentration of IL2 in the lymph node compartment');
catch 
    first_call = false;
end
% Add Hill Functions for APC/mAPC
if (isTreg)
    H_APC = 'H_APC';
    % Add Treg Variable
    addrule(model,['Tregs_ = V_T.' species_name],'repeatedAssignment'); 
else
    H_APC = 'H_mAPC';
end
    
% Add Parameters
div = addparameter(model,'div',params.div.Value,'ValueUnits',params.div.Units);
    set(div,'Notes',['T Cell Diversity ',params.div.Notes]);
n_clones = addparameter(model,'n_clones',params.n_clones.Value,'ValueUnits',params.n_clones.Units);
    set(n_clones,'Notes',['Number of T cell clones ' params.n_clones.Notes]);
q_LN_in = addparameter(model,'q_LN_in',params.q_LN_in.Value,'ValueUnits',params.q_LN_in.Units);
    set(q_LN_in,'Notes',['Rate of naive T cell transport into the LN ' params.q_LN_in.Notes]);
q_LN_out = addparameter(model,'q_LN_out',params.q_LN_out.Value,'ValueUnits',params.q_LN_out.Units);
    set(q_LN_out,'Notes',['Rate of naive T cell transport out of the LN ' params.q_LN_out.Notes]);
k_act = addparameter(model,'k_act',params.k_act.Value,'ValueUnits',params.k_act.Units);
    set(k_act,'Notes',[species_name ' activation rate']);
k_pro = addparameter(model,'k_pro',params.k_pro.Value,'ValueUnits',params.k_pro.Units);
    set(k_pro,'Notes',[species_name ' proliferation rate']);
k_death = addparameter(model,'k_death',params.k_death.Value,'ValueUnits',params.k_death.Units);
    set(k_death,'Notes',[species_name ' death rate']);
q_P_in = addparameter(model,'q_P_in',params.q_P_in.Value,'ValueUnits',params.q_P_in.Units);
    set(q_P_in,'Notes',['rate of ' species_name ' tranport into the peripheral compartment']);   
q_P_out = addparameter(model,'q_P_out',params.q_P_out.Value,'ValueUnits',params.q_P_out.Units);
    set(q_P_out,'Notes',['rate of ' species_name ' tranport out of the peripheral compartment']); 
q_T_in = addparameter(model,'q_T_in',params.q_T_in.Value,'ValueUnits',params.q_T_in.Units);
    set(q_T_in,'Notes',['rate of ' species_name ' tranport into the tumour compartment']);  
Q_nT_thym = addparameter(model,'Q_nT_thym',params.Q_nT_thym.Value,'ValueUnits',params.Q_nT_thym.Units);
    set(Q_nT_thym,'Notes',['Thymic output of n',species_name,' ', params.Q_nT_thym.Notes]);
k_nT_pro = addparameter(model,'k_nT_pro',params.k_nT_pro.Value,'ValueUnits',params.k_nT_pro.Units);
    set(k_nT_pro,'Notes',['rate of n',species_name,' proliferation ', params.k_nT_pro.Notes]);
K_nT_pro = addparameter(model,'K_nT_pro',params.K_nT_pro.Value,'ValueUnits',params.K_nT_pro.Units);
    set(K_nT_pro,'Notes',['half-maximal peripheral proliferation of of n',species_name,' ', params.K_nT_pro.Notes]);
k_nT_death = addparameter(model,'k_nT_death',params.k_nT_death.Value,'ValueUnits',params.k_nT_death.Units);
    set(k_nT_death,'Notes',['rate of n',species_name,' death ', params.k_nT_death.Notes]);
    
% IL2 Parameters
if (first_call)    
    p = addparameter(model,'k_IL2_deg',params.k_IL2_deg.Value,'ValueUnits',params.k_IL2_deg.Units);
        set(p,'Notes',['rate of IL2 degradation ' params.k_IL2_deg.Notes]);
    p = addparameter(model,'k_IL2_cons',params.k_IL2_cons.Value,'ValueUnits',params.k_IL2_cons.Units);
        set(p,'Notes',['rate of IL2 consumption by T cells ' params.k_IL2_cons.Notes]);
    p = addparameter(model,'k_IL2_sec',params.k_IL2_sec.Value,'ValueUnits',params.k_IL2_sec.Units);
        set(p,'Notes',['rate of IL2 secretion from T cells ' params.k_IL2_sec.Notes]);
    p = addparameter(model,'IL2_50',params.IL2_50.Value,'ValueUnits',params.IL2_50.Units);
        set(p,'Notes',['T cell activation half-maximal IL2 concentration ' params.IL2_50.Notes]);
    p = addparameter(model,'IL2_50_Treg',params.IL2_50_Treg.Value,'ValueUnits',params.IL2_50_Treg.Units);
        set(p,'Notes',['Treg activation half-maximal IL2 concentration ' params.IL2_50_Treg.Notes]);
    p = addparameter(model,'N0',params.N0.Value,'ValueUnits',params.N0.Units);
        set(p,'Notes',['numer of activated T cell generation with no IL2 ' params.N0.Notes]);
    p = addparameter(model,'N_IL2',params.N_IL2.Value,'ValueUnits',params.N_IL2.Units);
        set(p,'Notes',['maximum number of activated T cell generations due to IL2 ' params.N_IL2.Notes]);
end
if (~isTreg)
    k_Tcell = addparameter(model,'k_Tcell',params.k_Tcell.Value,'ValueUnits',params.k_Tcell.Units);
        set(k_Tcell,'Notes',['Rate of T cell exhaustion by cancer cells ' params.k_Tcell.Notes]);
    k_C_Tcell = addparameter(model,'k_C_Tcell',params.k_C_Tcell.Value,'ValueUnits',params.k_C_Tcell.Units);
        set(k_C_Tcell,'Notes',['Rate of cancer cell death by T cells ' params.k_C_Tcell.Notes]);
end
try % only add once
    k_Treg = addparameter(model,'k_Treg',params.k_Treg.Value,'ValueUnits',params.k_Treg.Units);
        set(k_Treg,'Notes',['Rate of T cell death by Tregs ' params.k_Treg.Notes]);
catch
end
% Antigen Default Hill Function
addparameter(model,['H_' antigen],1,'ValueUnits','dimensionless','ConstantValue',false);
    
% Add Reactions
% Thymic output of naive T cell
reaction = addreaction(model,'null -> V_C.nT');
    set(reaction,'ReactionRate','Q_nT_thym/div');
    set(reaction,'Notes','Thymic output of naive T cell to blood');
% Naive T cell proliferation in the peripheral compartment and the central compartment
reaction = addreaction(model,'null -> V_P.nT');
    set(reaction,'ReactionRate','k_nT_pro/div*V_P.nT/(K_nT_pro/div+V_P.nT)');
    set(reaction,'Notes','Naive T cell proliferation in the peripheral compartment');
reaction = addreaction(model,'null -> V_C.nT');
    set(reaction,'ReactionRate','k_nT_pro/div*V_C.nT/(K_nT_pro/div+V_C.nT)');
    set(reaction,'Notes','Naive T cell proliferation in the central compartment');
% Naive T cell death in the peripheral compartment and the central compartment
reaction = addreaction(model,'V_P.nT -> null');
    set(reaction,'ReactionRate','k_nT_death*V_P.nT');
    set(reaction,'Notes','Naive T cell death in the peripheral compartment');
reaction = addreaction(model,'V_C.nT -> null');
    set(reaction,'ReactionRate','k_nT_death*V_C.nT');
    set(reaction,'Notes','Naive T cell death in the central compartment');
% Naive T cell transport into and out of the peripheral compartment
reaction = addreaction(model,'V_C.nT -> V_P.nT');
    set(reaction,'ReactionRate','q_P_in*V_C.nT');
    set(reaction,'Notes','Naive T cell entry into the peripheral compartment');
reaction = addreaction(model,'V_P.nT -> V_C.nT');
    set(reaction,'ReactionRate','q_P_out*V_P.nT');
    set(reaction,'Notes','Naive T cell exit from the peripheral compartment');
% Naive T cell transport into and out of the lymph node
reaction = addreaction(model,'V_C.nT -> V_LN.nT');
    set(reaction,'ReactionRate','q_LN_in*V_C.nT');
    set(reaction,'Notes','Naive T cell entry into the lymph node');
reaction = addreaction(model,'V_LN.nT -> V_C.nT');
    set(reaction,'ReactionRate','q_LN_out*V_LN.nT');
    set(reaction,'Notes','Naive T cell exit from the lymph node');
    
% Naive T cell activation
reaction = addreaction(model,'V_LN.nT -> null');
    set(reaction,'ReactionRate',['k_act*' H_APC '*H_' antigen '*V_LN.nT']); 
    set(reaction,'Notes','Naive T cell activation');
    reaction = addreaction(model,'null -> V_LN.aT');
    set(reaction,'ReactionRate',['k_act*' H_APC '*H_' antigen '*V_LN.nT*n_clones']);
    set(reaction,'Notes','Naive T cell activation');
% Activated T cell proliferation <verify with Hanwen>
reaction = addreaction(model,'V_LN.aT -> null');
    set(reaction,'ReactionRate','k_pro*V_LN.aT');
    set(reaction,'Notes',['a' species_name ' cell proliferation']);
reaction = addreaction(model,'null -> V_LN.T');
    set(reaction,'ReactionRate','k_pro*2^N_aT*V_LN.aT'); 
    set(reaction,'Notes',['a' species_name ' cell proliferation']);
% T cell Death
reaction = addreaction(model,'V_C.T -> null');
    set(reaction,'ReactionRate','k_death*V_C.T');
    set(reaction,'Notes','T cell death in the central compartment');
reaction = addreaction(model,'V_P.T -> null');
    set(reaction,'ReactionRate','k_death*V_P.T');
    set(reaction,'Notes','T cell death in the peripheral compartment');
reaction = addreaction(model,'V_T.T -> V_T.T_exh');
    set(reaction,'ReactionRate','k_death*V_T.T');
    set(reaction,'Notes','T cell death in the tumour compartment');
reaction = addreaction(model,'V_LN.T -> null');
    set(reaction,'ReactionRate','k_death*V_LN.T');
    set(reaction,'Notes','T cell death in the lymph node compartment');
if (~isTreg)    
% T cell death from Treg
reaction = addreaction(model,'V_T.T -> V_T.T_exh');
    set(reaction,'ReactionRate','k_Treg*V_T.T*Tregs_/(C_total+T_total+cell)');
    set(reaction,'Notes','T cell death from Tregs');
% T cell death from cancer
reaction = addreaction(model,'V_T.T -> V_T.T_exh');
    set(reaction,'ReactionRate','k_Tcell*V_T.T*C_total/(C_total+T_total+cell)*H_PD1_C1');
    set(reaction,'Notes','T cell death from cancer');
end
% T cell transport
% Central & Peripheral
reaction = addreaction(model,'V_C.T -> V_P.T');
    set(reaction,'ReactionRate','q_P_in*V_C.T');
    set(reaction,'Notes','T cell transport into the peripheral compartment');
reaction = addreaction(model,'V_P.T -> V_C.T');
    set(reaction,'ReactionRate','q_P_out*V_P.T');
    set(reaction,'Notes','T cell transport out of the peripheral compartment');
% Central & Tumour
reaction = addreaction(model,'V_C.T -> V_T.T');
    set(reaction,'ReactionRate','q_T_in*V_T*V_C.T');
    set(reaction,'Notes','T cell transport into the tumour compartment');
% Central & LN
reaction = addreaction(model,'V_LN.T -> V_C.T');
    set(reaction,'ReactionRate','q_LN_out*V_LN.T');
    set(reaction,'Notes','T cell transport out of the lymph node compartment'); 
% IL2 Reactions
if (first_call)
    % IL2 Degradation
    reaction = addreaction(model,'V_LN.IL2 -> null');
        set(reaction,'ReactionRate','k_IL2_deg*V_LN.IL2*V_LN');
        set(reaction,'Notes','IL2 degradation');
    % IL2 Consumption
    reaction = addreaction(model,'V_LN.IL2 -> null');
        set(reaction,'ReactionRate','k_IL2_cons*T_total_LN*IL2/(IL2_50+IL2)');
        set(reaction,'Notes','IL2 consumption by T cells');
end
if (isTreg)
    % IL2 Consumption by Tregs
    reaction = addreaction(model,'V_LN.IL2 -> null');
        set(reaction,'ReactionRate','k_IL2_cons*V_LN.T*IL2/(IL2_50_Treg+IL2)');
        set(reaction,'Notes','IL2 consumption by Tregs');
else
    % IL2 Secretion by Activated T Cells
    reaction = addreaction(model,'null -> V_LN.IL2');
       set(reaction,'ReactionRate','k_IL2_sec*V_LN.aT');
       set(reaction,'Notes','IL2 secretion from activated T cells');
end 

% Add Rules
if (first_call)
    % Set Number of Activated T Cell Generations
    addparameter(model,'N_aT',1,'ValueUnits','dimensionless','ConstantValue',false);
    addrule(model,'N_aT = N0 + N_IL2*IL2/(IL2_50+IL2)','repeatedAssignment');
end

% Get Model Rules for Updating
model_rules = get(model,'Rules');

% Update Tumour Volume (Rule 1)
volume_rule = model_rules(1);
rule = get(volume_rule,'Rule');
set(volume_rule,'Rule',[rule '+vol_Tcell*V_T.' species_name]);

% Update Total T Cells in Tumour (Rule 3)
Tcell_rule = model_rules(3);
rule = get(Tcell_rule,'Rule');
set(Tcell_rule,'Rule',[rule '+V_T.' species_name]);

% Update Total T Cells in LN (Rule 4)
if (~isTreg)
    Tcell_rule = model_rules(4);
    rule = get(Tcell_rule,'Rule');
    set(Tcell_rule,'Rule',[rule '+V_LN.' species_name]);
end

% Update Cancer Killing by T Cells (Rule 5)
if (exist('cancer_types','var') && ~isTreg)
    for i = 1:length(cancer_types)
        reaction = addreaction(model,['V_T.' cancer_types{i} ' -> V_T.C_x']);
            set(reaction,'ReactionRate',['k_C_Tcell*V_T.T*V_T.' cancer_types{i} '/(C_total+T_total+cell)*(1-H_PD1_C1)']);
            set(reaction,'Notes','Cancer cell killing by T cells');
        rule = get(model_rules(5),'Rule'); % Jan 8: edited Rule 5 k_Tcell -> k_C_Tcell
            set(model_rules(5),'Rule',[rule,'+k_C_Tcell*V_T.T*V_T.' cancer_types{i} '/(C_total+T_total+cell)*(1-H_PD1_C1)']);
    end
end
    
% Rename Objects with 'species_name'
rename(div,['div_' species_name]);
rename(n_clones,['n_' species_name '_clones']);
rename(nT_C,['n' species_name]);
rename(nT_P,['n' species_name]);
rename(nT_LN,['n' species_name]);
rename(aT,['a' species_name]);
rename(T_C,species_name);
rename(T_P,species_name);
rename(T_T,species_name);
rename(T_LN,species_name);
rename(q_LN_in,['q_LN_' species_name '_in']);
rename(q_LN_out,['q_LN_' species_name '_out']);
rename(k_act,['k_' species_name '_act']);
rename(k_pro,['k_' species_name '_pro']);
rename(k_death,['k_' species_name '_death']);
rename(q_P_in,['q_' species_name '_P_in']);
rename(q_P_out,['q_' species_name '_P_out']);
rename(q_T_in,['q_' species_name '_T_in']);
rename(q_LN_out,['q_' species_name '_LN_out']);
rename(Q_nT_thym,['Q_n' species_name '_thym']);
rename(k_nT_pro,['k_n' species_name '_pro']);
rename(K_nT_pro,['K_n' species_name '_pro']);
rename(k_nT_death,['k_n' species_name '_death']);
if (~isTreg)
    rename(k_Tcell,['k_' species_name]);
    rename(k_C_Tcell,['k_C_' species_name]);
end