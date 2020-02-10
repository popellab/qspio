% TCR Kintetic Proofreading with limited signaling Module
%
% Sub-module used by antigen module
%
% Inputs: model   -- SimBiology model object with four compartments 
%         ID      -- T cell-antigen ID number [must be unique]
%         params  -- object containing model parameter Values, Units, and Notes
%        
% Outputs: model -- SimBiology model object with new antigen module
%
% Created: May 08, 2019 (Mohammad Jafarnejad)
% Last Modified: May 08, 2019 (MJ)

function model = TCR_KPR_module(model,ID,params)

% Names
antigen_name = ['P' ID];
epitope_name = ['p' ID];
Tcell_name = ['T' ID];
i = 1;
Mp = ['M' num2str(i) epitope_name];

% add TCR parameters
k_TCR_on = addparameter(model,'k_TCR_on',params.k_TCR_on.Value,'ValueUnits',params.k_TCR_on.Units);
    set(k_TCR_on,'Notes',['Rate of TCR binding to MHC-peptide complex ' params.k_TCR_on.Notes]);
k_TCR_off = addparameter(model,'k_TCR_off',params.k_TCR_off.Value,'ValueUnits',params.k_TCR_off.Units);
    set(k_TCR_off,'Notes',['Rate of TCR unbinding from MHC-peptide complex ' params.k_TCR_off.Notes]);
k_TCR_p = addparameter(model,'k_TCR_p',params.k_TCR_p.Value,'ValueUnits',params.k_TCR_p.Units);
    set(k_TCR_p,'Notes',['Rate of MHC-peptide-TCR complex modification ' params.k_TCR_p.Notes]);
phi_TCR = addparameter(model,'phi_TCR',params.phi_TCR.Value,'ValueUnits',params.phi_TCR.Units);
    set(phi_TCR,'Notes',['Rate of MHC-peptide-TCR complex with maximal modification that leads to non-signaling' params.phi_TCR.Notes]);
N_TCR = addparameter(model,'N_TCR',params.N_TCR.Value,'ValueUnits',params.N_TCR.Units);
    set(N_TCR,'Notes',['Number of modifcation steps for MHC-peptide-TCR complex ' params.N_TCR.Notes]);
TCR_tot = addparameter(model,'TCR_tot',params.TCR_tot.Value,'ValueUnits',params.TCR_tot.Units);
    set(TCR_tot,'Notes',['Total number of TCR molecules per naive T cell ' params.TCR_tot.Notes]);
pTCR_MHC_tot = addparameter(model,'pTCR_MHC_tot',0,'ValueUnits',params.TCR_tot.Units,'ConstantValue',false);
    set(pTCR_MHC_tot,'Notes',['Total number of MHC-' epitope_name '-TCR complexes of all different activation levels']);

% This only works with 1 MHC at the moment
P_T = ['A_s.' Mp '/ n_' Tcell_name '_clones'];
K_D = 'k_TCR_off/k_TCR_on';
C_T   = ['(0.5*(' P_T ' + TCR_tot + ' K_D ' - TCR_tot*(((' P_T ' + TCR_tot + ' K_D ')/TCR_tot)^2 - 4*' P_T '/TCR_tot)^0.5))'];
alpha = '(k_TCR_p /( k_TCR_off + k_TCR_p))';
beta  = 'k_TCR_off /( k_TCR_off + phi_TCR)';

addrule(model,['pTCR_MHC_tot = ' beta ' * ' alpha '^N_TCR * ' C_T],'repeatedAssignment');
addrule(model,['H_' antigen_name ' = pTCR_MHC_tot/(pTCR_MHC_tot+' epitope_name '_50)'],'repeatedAssignment');

% Rename Objects 
rename(k_TCR_on,['k_' Mp '_TCR_on']);
rename(k_TCR_off,['k_' Mp '_TCR_off']);
rename(k_TCR_p,['k_' Mp '_TCR_p']);
rename(phi_TCR,['phi_' Mp '_TCR']);
rename(N_TCR,['N_' Mp '_TCR']);
rename(TCR_tot,['TCR_' epitope_name '_tot']);
rename(pTCR_MHC_tot,['pTCR_' epitope_name '_MHC_tot']);