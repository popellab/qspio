% TCR Receptor Occupancy Module
%
% Sub-module used by antigen module
%
% Inputs: model   -- SimBiology model object with four compartments 
%         ID      -- T cell-antigen ID number [must be unique]
%         params  -- object containing model parameter Values, Units, and Notes
%        
% Outputs: model -- SimBiology model object with new TCR module
%
% Created: May 13, 2019 (Mohammad Jafarnejad)
% Last Modified: May 13, 2019 (MJ)

function model = TCR_RO_module(model,ID,params)

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
TCR_tot = addparameter(model,'TCR_tot',params.TCR_tot.Value,'ValueUnits',params.TCR_tot.Units);
    set(TCR_tot,'Notes',['Total number of TCR molecules per naive T cell ' params.TCR_tot.Notes]);
pTCR_MHC_tot = addparameter(model,'pTCR_MHC_tot',0,'ValueUnits',params.TCR_tot.Units,'ConstantValue',false);
    set(pTCR_MHC_tot,'Notes',['Total number of MHC-' epitope_name '-TCR complexes of all different activation levels']);

% This only works with 1 MHC at the moment
P_T = ['A_s.' Mp '/n_' Tcell_name '_clones'];
K_D = 'k_TCR_off/k_TCR_on';
addrule(model,['pTCR_MHC_tot = 0.5*(' P_T ' + TCR_tot + ' K_D ' - TCR_tot *(((' P_T ' + TCR_tot + ' K_D ')/TCR_tot)^2 - 4*' P_T '/TCR_tot)^0.5)'],'repeatedAssignment');

addrule(model,['H_' antigen_name ' = pTCR_MHC_tot/(pTCR_MHC_tot+' epitope_name '_50)'],'repeatedAssignment');

% Rename Objects 
rename(k_TCR_on,['k_' Mp '_TCR_on']);
rename(k_TCR_off,['k_' Mp '_TCR_off']);
rename(TCR_tot,['TCR_' epitope_name '_tot']);
rename(pTCR_MHC_tot,['pTCR_' epitope_name '_MHC_tot']);