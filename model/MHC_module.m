% MHC Module
%
% Sub-module used by antigen module
%
% Inputs: model -- SimBiology model object with four compartments 
%         ID    -- ID number for MHC [must be unique]
%         ICs   -- array containing fraction of MHC in
%                  - edosomal compartment
%                  - surface compartment
%         MHC_T -- total amount of MHC [parameter]
%        
% Outputs: model -- SimBiology model object with new MHC module
%
% Created: Nov 21, 2018 (Richard Sové)
% Last Modified: Feb 06, 2019 (RJS)

function model = MHC_module(model,ID,MHC_T)
    
% Add MHC
M_e = addspecies(model.Compartment(6),'M',MHC_T.Value,'InitialAmountUnits',MHC_T.Units);
    set(M_e,'Notes','Amount of MHC per area on the endosomal surface');
M_s = addspecies(model.Compartment(7),'M',1e-6,'InitialAmountUnits',MHC_T.Units);
    set(M_e,'Notes','Amount of MHC per area on the cell surface');

% % Add Reactions
reaction = addreaction(model,'A_e.M -> A_s.M');
    set(reaction,'ReactionRate','kout*A_e.M*A_e-kin*A_s.M*A_s');
    set(reaction,'Notes','MHC translocation');

% Rename MHC Molecules
rename(M_e,['M' num2str(ID)]);
rename(M_s,['M' num2str(ID)]);