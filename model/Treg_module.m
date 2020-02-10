% Treg Module
%
% Models Treg transport and activation by APCs [Use before antigen corresponding module]
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
%        
% Outputs: model -- SimBiology model object with new Tcell module
%
% Created: Sep 23, 2019 (Richard Sové)
% Last Modified: Sep 23, 2019 (RJS)

function model = Treg_module(model,params)

% Get General T Cell Parameters
params_Treg.q_LN_in = params.q_LN_in;
params_Treg.q_LN_out = params.q_LN_out;
params_Treg.k_Treg = params.k_Treg;
params_Treg.nT_C = params.nT_C;
params_Treg.nT_P = params.nT_P;
params_Treg.nT_LN = params.nT_LN;
params_Treg.nTreg_C = params.nTreg_C;
params_Treg.nTreg_P = params.nTreg_P;
params_Treg.nTreg_LN = params.nTreg_LN;
params_Treg.K_nT_pro = params.K_nT_pro;
params_Treg.k_nT_death = params.k_nT_death;
params_Treg.k_IL2_sec = params.k_IL2_sec;
params_Treg.k_IL2_cons = params.k_IL2_cons;
params_Treg.IL2_50 = params.IL2_50;
params_Treg.IL2_50_Treg = params.IL2_50_Treg;
params_Treg.N0 = params.N0;
params_Treg.N_IL2 = params.N_IL2;
params_Treg.k_Treg = params.k_Treg;

% Get Rename Treg-specific Parameters
params_Treg.div = params.div_Treg;
params_Treg.n_clones = params.n_clones_Treg;
params_Treg.k_act = params.k_act_Treg;
params_Treg.k_pro = params.k_pro_Treg;
params_Treg.k_death = params.k_death_Treg;
params_Treg.q_P_in = params.q_P_in_Treg;
params_Treg.q_P_out = params.q_P_out_Treg;
params_Treg.q_T_in = params.q_T_in_Treg;
params_Treg.Q_nT_thym = params.Q_nT_thym_Treg;
params_Treg.k_nT_pro = params.k_nT_pro_Treg;

% Call Tcell Module
model = Tcell_module(model,'0',params_Treg,[]);