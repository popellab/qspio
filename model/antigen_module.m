% Antigen Module
%
% Models antigen presentation on APCs 
%
% Requirements: cancer_module, Tcell_module and APC_module
%
% Inputs: model   -- SimBiology model object with four compartments 
%         ID      -- T cell-antigen ID number [must be unique]
%         params  -- object containing model parameter Values, Units, and Notes:
%                    - N_MHC--number of types of MHC molecules
%                    - MHC_T--total amount of MHC per area
%                    - kin--rate of MHC internalization
%                    - kout--rate of MHC exernalization
%                    - V_e--endosomal volume
%                    - A_e--endosomal surface area
%                    - A_s--APC surface area
%                    - k_up--rate of antigen uptake by APCs
%                    - k_xP_deg--rate of extracellular antigen degradation
%                    - k_P_deg--rate of endosomal antigen degradation
%                    - k_p_deg--rate of endosomal epitope degradation
%                    - k_on--rate of antigen_MHC binding
%                    - p_50--epitope concentration for half-maximal T cell activation
%         antigen -- antigen structure
%        
% Outputs: model -- SimBiology model object with new antigen module
%
% Created: Oct 29, 2018 (Richard Sové)
% Last Modified: May 08, 2019 (MJ)

function model = antigen_module(model,ID,params,antigen)

% Names
antigen_name = ['P' ID];
epitope_name = ['p' ID];

% Get Number of T Cell Clones
nTcells = howManyClones(model);

% Add MHCs on First Call 
N_MHC = params.N_MHC.Value;
if (length(model.Compartment)<5) % there should be 4 compartments before this module is added
    first_call = true;
else
    first_call = false;
end
if (first_call)
    % Add Endosomal And Surface Compartments
    comp = addcompartment(model,'V_e',params.V_e.Value,'CapacityUnits',params.V_e.Units); 
        set(comp,'Notes',['APC endosomal comparment ' params.V_e.Notes]);
    comp = addcompartment(model,'A_e',params.A_e.Value,'CapacityUnits',params.A_e.Units);
        set(comp,'Notes',['APC endosomal surface comparment ' params.A_e.Notes]);
    comp = addcompartment(model,'A_s',params.A_s.Value,'CapacityUnits',params.A_s.Units); 
        set(comp,'Notes',['APC surface comparment ' params.A_s.Notes]);
     
    % Add Translocation Rates
    kin = addparameter(model,'kin',params.kin.Value,'ValueUnits',params.kin.Units);
        set(kin,'Notes',['Rate of MHC internalization ' params.kin.Notes]);
    kout = addparameter(model,'kout',params.kout.Value,'ValueUnits',params.kout.Units);
        set(kout,'Notes',['Rate of MHC externalization ' params.kout.Notes]);
    
    % Add MHCs
    for i = 1:N_MHC
        model = MHC_module(model,i,params.MHC_T);
    end
end

% Add Antigen and Epitope
Pf = addspecies(model.Compartment(3),'P',1e-18,'InitialAmountUnits','molarity');
    set(Pf,'Notes',['Concentration of  free antigen (' antigen_name ') in the T compartment']);
Pe = addspecies(model.Compartment(5),'P',1e-18,'InitialAmountUnits','molarity');
    set(Pe,'Notes',['Concentration of antigen ' antigen_name ' in the APC endosomes']);
p = addspecies(model.Compartment(5),'p',1e-18,'InitialAmountUnits','molarity');
    set(p,'Notes',['Concentration of epitope ' antigen_name ' in the APC endosomes']);
    
% Add Parameters
k_up = addparameter(model,'k_up',params.k_up.Value,'ValueUnits',params.k_up.Units);
    set(k_up,'Notes',['Rate of antigen uptake by APCs ' params.k_up.Notes]);
k_xP_deg = addparameter(model,'k_xP_deg', params.k_xP_deg.Value,'ValueUnits',params.k_xP_deg.Units);
    set(k_xP_deg,'Notes',['Rate of extracellular antigen degradation ' params.k_xP_deg.Notes]);  
k_P_deg = addparameter(model,'k_P_deg',params.k_P_deg.Value,'ValueUnits',params.k_P_deg.Units);
    set(k_P_deg,'Notes',['Rate of endosomal antigen degradation ' params.k_P_deg.Notes]);
k_p_deg = addparameter(model,'k_p_deg',params.k_p_deg.Value,'ValueUnits',params.k_p_deg.Units);
    set(k_p_deg,'Notes',['Rate of endosomal epitope degradation ' params.k_p_deg.Notes]);
k_on = addparameter(model,'k_on',params.k_on.Value,'ValueUnits',params.k_on.Units);
    set(k_on,'Notes',['Rate of antigen-MHC binding ' params.k_on.Notes]);
for i = 1:N_MHC
    k_d = addparameter(model,['k_' antigen_name '_d' num2str(i)],antigen.kd(i).Value,'ValueUnits',antigen.kd(i).Units);
        set(k_d,'Notes','Antigen-MHC kd');
end
p_50 = addparameter(model,[epitope_name '_50'],params.p_50.Value,'ValueUnits',params.p_50.Units);
    set(p_50,'Notes',params.p_50.Notes);
k_dep = antigen_rate(model,ID,antigen.cancers,antigen.concentration,nTcells);

% Add Reactions
% Antigen 
% Antigen Source
reaction = addreaction(model,'null -> V_T.P');
    set(reaction,'ReactionRate',k_dep);
    set(reaction,'Notes','Antigen deposition from dying cancer cells');
% Free Antigen Degradation
reaction = addreaction(model,'V_T.P -> null');
    set(reaction,'ReactionRate','k_xP_deg*V_T.P*V_T');
    set(reaction,'Notes','Free antigen degradation');
% Antigen Uptake by mAPCs 
reaction = addreaction(model,'V_T.P -> null');
    set(reaction,'ReactionRate','k_up*V_LN.mAPC*V_T.P*V_T');
    set(reaction,'Notes','Antigen uptake by mature antigen presenting cells');
reaction = addreaction(model,'null -> V_e.P');
    set(reaction,'ReactionRate','k_up*cell*V_T.P*V_e');
    set(reaction,'Notes','Antigen uptake by mature antigen presenting cells');
% Antigen Degradation in APC Endosomes    
reaction = addreaction(model,'V_e.P -> V_e.p');
    set(reaction,'ReactionRate','k_P_deg*V_e.P*V_e');
    set(reaction,'Notes','Antigen degradation in APC endosomes');
% Epitope
reaction = addreaction(model,'V_e.p -> null');
    set(reaction,'ReactionRate','k_p_deg*V_e.p*V_e');
    set(reaction,'Notes','Epitope degradation in APC endosomes');
% MHC-Epitope Complexes
for i = 1:N_MHC
    Mp_e = addspecies(model.Compartment(6),'Mp',1e-6,'InitialAmountUnits',params.MHC_T.Units);
        set(Mp_e,'Notes','Antigen-MHC complex');
    Mp_s = addspecies(model.Compartment(7),'Mp',1e-6,'InitialAmountUnits',params.MHC_T.Units);
        set(Mp_s,'Notes','Antigen-MHC complex');
    % Binding in Endosome
    reaction = addreaction(model,['V_e.p + A_e.M' num2str(i) ' -> A_e.Mp']);
        set(reaction,'ReactionRate',['k_on*V_e.p*A_e.M' num2str(i) '*A_e']);
        set(reaction,'Notes','Antigen-MHC binding in endosome');
    % Unbinding in Endosome
    reaction = addreaction(model,['A_e.Mp -> V_e.p + A_e.M' num2str(i)]);
        set(reaction,'ReactionRate',['k_' antigen_name '_d' num2str(i) '*k_on*A_e.Mp*A_e']);
        set(reaction,'Notes','Antigen-MHC unbinding in endosome');
    % Unbinding on Surface
    reaction = addreaction(model,['A_s.Mp -> A_s.M' num2str(i)]);
        set(reaction,'ReactionRate',['k_' antigen_name '_d' num2str(i) '*k_on*A_s.Mp*A_s']);
        set(reaction,'Notes','Antigen-MHC unbinding on APC surface');
    % Antigen-MHC Translocation
    reaction = addreaction(model,'A_e.Mp -> A_s.Mp');
        set(reaction,'ReactionRate','kout*A_e.Mp*A_e');
%         set(reaction,'ReactionRate','kout*A_e.Mp*A_e - kin*A_s.Mp*A_s');
        set(reaction,'Notes','Antigen-MHC translocation');
    rename(Mp_e,['M' num2str(i) epitope_name]);
    rename(Mp_s,['M' num2str(i) epitope_name]);
end

% Determine if checkpoint sizes have been defined     
first_call = true;
try % see if synapse exist
    parameter = addparameter(model,'A_syn' ,params.A_syn.Value ,'ValueUnits',params.A_syn.Units);
    set(parameter,'Notes',['Surface area of the synapse ' params.A_syn.Notes]);
catch 
    first_call = false;
end
if first_call
% Add surface areas
parameter = addparameter(model,'A_Tcell' ,params.A_Tcell.Value ,'ValueUnits',params.A_Tcell.Units);
    set(parameter,'Notes',['Surface area of the T cell ' params.A_Tcell.Notes]);
parameter = addparameter(model,'A_cell' ,params.A_cell.Value ,'ValueUnits',params.A_cell.Units);
    set(parameter,'Notes',['Surface area of the Cancer cell ' params.A_cell.Notes]);
% *** This parameter is double defined ***
parameter = addparameter(model,'A_APC' ,params.A_s.Value ,'ValueUnits',params.A_s.Units);
    set(parameter,'Notes',['Surface area of the APC ' params.A_s.Notes]);
end

% Update Hill Function for Antigens
if ID == '0'
    % Adds TCR receptor occupancy to the model for Treg
    model = TCR_RO_module(model,ID,params);
else
    % Adds TCR kinetic proofreading to the model for Teff
    %(it only works with 1 MHC at the moment)
    model = TCR_KPR_module(model,ID,params);
%     model = TCR_RO_module(model,ID,params);
end

% Rename Objects with 'species_name'
rename(Pf,antigen_name);
rename(Pe,antigen_name);
rename(p,epitope_name);
rename(k_up,['k_' antigen_name '_up']);
rename(k_xP_deg,['k_x' antigen_name '_deg']);
rename(k_P_deg,['k_' antigen_name '_deg']);
rename(k_p_deg,['k_' epitope_name '_deg']);
rename(k_on,['k_' antigen_name '_on']);