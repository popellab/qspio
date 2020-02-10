% Checkpoint Module
%
% Models PD1 Interactions
%
% Inputs: model        -- simbio model object with four compartments 
%         params       -- object containing the default parameters 
%         Tname        -- name of the T cell forming the checkpoint synapse 
%         Cname        -- name of the cancer or APC cell forming the 
%                        checkpoint synapse 
%         varargin     -- name-value pair:
%                         - drugName: string (or string array [max 2]) of drug
%                         name(s) for PD1 or PDL1 antibody 
%
% Outputs: model -- simbio model object with new PD1 module
%
% Note: This only works for Ti (i >= 1) interaction with APC and Cj. We
% would need to add a check for Tregs if we would like to generalize it to checkpoint
% model of Treg.
%
% Created: Jun 11, 2019 (Mohammad Jafarnejad)
% Last Modified: Oct 02 2019 (RJS)

function model = PD1_module(model,params,Tname,Cname,varargin)

% Parse varargin
ip = inputParser;
addParameter(ip,'drugName','');
parse(ip,varargin{:});
drugName = ip.Results.drugName;

% select the right compartment based on cancer or APC
if Cname(1)=='C'
    compDrug = model.Compartment(3);
    gamma = 'gamma_C';
elseif Cname(1)=='A'
    compDrug = model.Compartment(4);
    gamma = 'gamma_LN';
end

% Add the synapse compartment
comp = addcompartment(model,['syn_',Tname,'_',Cname],params.A_syn.Value,'CapacityUnits',params.A_syn.Units); 
    set(comp,'Notes',['synapse comparment between ',Tname,' and ',Cname,' ', params.A_syn.Notes]);

% Determine if checkpoint sizes have been defined before in antigen module    
first_call = true;
try % see if synapse exist
    p = addparameter(model,'A_syn' ,params.A_syn.Value ,'ValueUnits',params.A_syn.Units);
    set(p,'Notes',['Surface area of the synapse ' params.A_syn.Notes]);
catch 
    first_call = false;
end
if first_call
% Add surface areas
p = addparameter(model,'A_Tcell' ,params.A_Tcell.Value ,'ValueUnits',params.A_Tcell.Units);
    set(p,'Notes',['Surface area of the T cell ' params.A_Tcell.Notes]);
p = addparameter(model,'A_cell' ,params.A_cell.Value ,'ValueUnits',params.A_cell.Units);
    set(p,'Notes',['Surface area of the Cancer cell ' params.A_cell.Notes]);
p = addparameter(model,'A_APC' ,params.A_s.Value ,'ValueUnits',params.A_s.Units);
    set(p,'Notes',['Surface area of the APC ' params.A_s.Notes]);
end

% Determine if first call    
first_call = true;
if first_call
% Checkpoint Binding    
kon = addparameter(model,'kon_PD1_PDL1',params.kon_PD1_PDL1.Value,'ValueUnits',params.kon_PD1_PDL1.Units);
    set(kon,'Notes',['kon of PD1-PDL1 binding ' params.kon_PD1_PDL1.Notes]);

% Add Antibodies
p = addparameter(model,'aPD1',0,'ValueUnits','molarity','ConstantValue',false);
    set(p,'Notes','Concentration of PD1 antibody');
p = addparameter(model,'aPDL1',0,'ValueUnits','molarity','ConstantValue',false);
    set(p,'Notes','Concentration of PDL1 antibody'); 
% Add Pharmacokinetics 
% Set Default Gamma
gam.aPD1 = '1';
gam.aPDL1 = '1';
if iscell(drugName) % if drug names are given in cell array
    if length(drugName) > 2
        err('More than two drugs were specified in PD1 module.');
    elseif length(drugName) == 1
        params_pk = pk_parameters(drugName{1});
        model = pk_module(model,drugName{1},params_pk);
        addrule(model,[params_pk.type ' = ' compDrug.Name '.' drugName{1}],'repeatedAssignment');
        gam.(params_pk.type) = [gamma '_' drugName{1}];
    elseif length(drugName) == 2
        % Drug 1
        params_pk1 = pk_parameters(drugName{1});
        model = pk_module(model,drugName{1},params_pk1);
        addrule(model,[params_pk1.type ' = ' compDrug.Name '.' drugName{1}],'repeatedAssignment');
        gam.(params_pk1.type) = [gamma '_' drugName{1}];
        % Drug 2
        params_pk2 = pk_parameters(drugName{2});
        model = pk_module(model,drugName{2},params_pk2);
        addrule(model,[params_pk2.type ' = ' compDrug.Name '.' drugName{2}],'repeatedAssignment');
        gam.(params_pk2.type) = [gamma '_' drugName{2}];
    end
else % if drug name is given as single string
    params_pk = pk_parameters(drugName);
    model = pk_module(model,drugName,params_pk);
    addrule(model,[params_pk.type ' = ' compDrug.Name '.' drugName],'repeatedAssignment');
    gam.(params_pk.type) = [gamma '_' drugName];
end

% Add kon Values
kon = addparameter(model,'kon_PD1_PDL2',params.kon_PD1_PDL2.Value,'ValueUnits',params.kon_PD1_PDL2.Units);
    set(kon,'Notes',['kon of PD1-PDL2 binding ' params.kon_PD1_PDL2.Notes]);
kon = addparameter(model,'kon_PD1_aPD1',params.kon_PD1_aPD1.Value,'ValueUnits',params.kon_PD1_aPD1.Units);
    set(kon,'Notes',['kon of PD1-antiPD1 binding ' params.kon_PD1_aPD1.Notes]);
kon = addparameter(model,'kon_PDL1_aPDL1',params.kon_PDL1_aPDL1.Value,'ValueUnits',params.kon_PDL1_aPDL1.Units);
    set(kon,'Notes',['kon of PDL1-antiPDL1 binding ' params.kon_PDL1_aPDL1.Notes]);    
 
% Add koff Values
koff = addparameter(model,'koff_PD1_PDL1' ,params.koff_PD1_PDL1.Value ,'ValueUnits',params.koff_PD1_PDL1.Units);
    set(koff,'Notes',['koff of PD1-PDL1 binding ' params.koff_PD1_PDL1.Notes]);
koff = addparameter(model,'koff_PD1_PDL2' ,params.koff_PD1_PDL2.Value ,'ValueUnits',params.koff_PD1_PDL2.Units);
    set(koff,'Notes',['koff of PD1-PDL2 binding ' params.koff_PD1_PDL2.Notes]);
koff = addparameter(model,'koff_PD1_aPD1' ,params.koff_PD1_aPD1.Value ,'ValueUnits',params.koff_PD1_aPD1.Units);
    set(koff,'Notes',['koff of PD1-nivolumab binding ' params.koff_PD1_aPD1.Notes]);
koff = addparameter(model,'koff_PDL1_aPDL1',params.koff_PDL1_aPDL1.Value,'ValueUnits',params.koff_PDL1_aPDL1.Units);
    set(koff,'Notes',['koff of PDL1-durvalumab binding ' params.koff_PDL1_aPDL1.Notes]);

% Bivalent anibody parameters  
p = addparameter(model,'Chi_PD1' ,params.Chi_PD1.Value ,'ValueUnits',params.Chi_PD1.Units);
    set(p,'Notes',['Antibody cross-arm binding efficiency ' params.Chi_PD1.Notes]);
p = addparameter(model,'Chi_PDL1' ,params.Chi_PDL1.Value ,'ValueUnits',params.Chi_PDL1.Units);
    set(p,'Notes',['Antibody cross-arm binding efficiency ' params.Chi_PDL1.Notes]);

% PD1-related Hill parameters
p = addparameter(model,'PD1_50',params.PD1_50.Value,'ValueUnits',params.PD1_50.Units);
    set(p,'Notes',['PD1/PDL1 concentration for half-maximal T cell inactivation ' params.PD1_50.Notes]);
p = addparameter(model,'n_PD1',params.n_PD1.Value,'ValueUnits',params.n_PD1.Units);
    set(p,'Notes',['Hill coefficient for PD1/PDL1 half-maximal T cell inactivation ' params.n_PD1.Notes]);   
 
end

% Checkpoint Expressions
if first_call
    p = addparameter(model,'T_PD1_total',params.T_PD1_total.Value,'ValueUnits',params.T_PD1_total.Units,'ConstantValue',false); 
        set(p,'Notes',['Concentration of PD1 on T cells ' params.T_PD1_total.Notes]);
    p = addparameter(model,'C_PDL1_total',params.C_PDL1_total.Value,'ValueUnits',params.C_PDL1_total.Units,'ConstantValue',false);
        set(p,'Notes',['Number of PDL1 molecules per cancer cell ' params.C_PDL1_total.Notes]);
    p = addparameter(model,'C_PDL2_total',params.C_PDL2_total.Value,'ValueUnits',params.C_PDL2_total.Units,'ConstantValue',false);
        set(p,'Notes',['Number of PDL2 molecules per cancer cell ' params.C_PDL2_total.Notes]);
    p = addparameter(model,'APC_PDL1_total',params.APC_PDL1_total.Value,'ValueUnits',params.APC_PDL1_total.Units,'ConstantValue',false);
        set(p,'Notes',['Number of PDL1 molecules per APC ' params.C_PDL1_total.Notes]);
    p = addparameter(model,'APC_PDL2_total',params.APC_PDL2_total.Value,'ValueUnits',params.APC_PDL2_total.Units,'ConstantValue',false);
        set(p,'Notes',['Number of PDL2 molecules per APC ' params.APC_PDL2_total.Notes]);
end  
 
% Add Species
x = addspecies(comp,'PD1_PDL1',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PD1-PDL1 complex');
x = addspecies(comp,'PD1_PDL2',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PD1-PDL2 complex');
x = addspecies(comp,'PD1',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PD1 in synapse');
x = addspecies(comp,'PDL1',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PDL1 in synapse');
x = addspecies(comp,'PDL2',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PDL2 in synapse');
x = addspecies(comp,'PD1_aPD1',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PD1-nivolumab complex');
x = addspecies(comp,'PD1_aPD1_PD1',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PD1-nivolumab-PD1 complex');
x = addspecies(comp,'PDL1_aPDL1',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PDL1-durvalumab complex');
x = addspecies(comp,'PDL1_aPDL1_PDL1',0,'InitialAmountUnits','molecule/micrometer^2');
    set(x,'Notes','concentration of PDL1-durvalumab-PDL1 complex');
    
% Update Input Parameters
addrule(model,[comp.Name,'.PD1',' = T_PD1_total /A_Tcell'   ] ,'initialAssignment');
if Cname(1)=='C'
    addrule(model,[comp.Name,'.PDL1 = C_PDL1_total /A_cell'] ,'initialAssignment');
    addrule(model,[comp.Name,'.PDL2 = C_PDL2_total /A_cell'] ,'initialAssignment');
elseif Cname(1)=='A'
    addrule(model,[comp.Name,'.PDL1 = APC_PDL1_total /A_APC'] ,'initialAssignment');
    addrule(model,[comp.Name,'.PDL2 = APC_PDL2_total /A_APC'] ,'initialAssignment');
end
 
% Dynamics of PD1/PDL1/PDL2/aPD1/aPDL1a
% PD1-PDL1
R = addreaction(model,[comp.Name,'.PD1 + ',comp.Name,'.PDL1 <-> ',comp.Name,'.PD1_PDL1']);
    set (R, 'ReactionRate', ['kon_PD1_PDL1*(',comp.Name,'.PD1)*(',comp.Name,'.PDL1)  -  koff_PD1_PDL1*',comp.Name,'.PD1_PDL1']);
    set (R, 'Notes'       , 'binding and unbinding of PD1 PDL1 in synapse');
% PD1-PDL2
R = addreaction(model,[comp.Name,'.PD1 + ',comp.Name,'.PDL2 <-> ',comp.Name,'.PD1_PDL2']);
    set (R,'ReactionRate',['kon_PD1_PDL2*(',comp.Name,'.PD1)*(',comp.Name,'.PDL2)  -  koff_PD1_PDL2*',comp.Name,'.PD1_PDL2']);
    set (R,'Notes','binding and unbinding of PD1 PDL2 in synapse'); 
% PD1-aPD1 
R = addreaction(model,[comp.Name,'.PD1 <-> ',comp.Name,'.PD1_aPD1']);
    set (R,'ReactionRate',['2*kon_PD1_aPD1*(',comp.Name,'.PD1*aPD1/',gam.aPD1,') -  koff_PD1_aPD1*',comp.Name,'.PD1_aPD1']);
    set (R,'Notes',['binding and unbinding of PD1 to anti-PD1 on ',Tname,' surface in synapse']);
% PD1-aPD1-PD1
R = addreaction(model,[comp.Name,'.PD1_aPD1 + ',comp.Name,'.PD1 <-> ',comp.Name,'.PD1_aPD1_PD1']);
    set (R,'ReactionRate',['Chi_PD1*kon_PD1_aPD1*(',comp.Name,'.PD1 * ',comp.Name,'.PD1_aPD1) -  2*koff_PD1_aPD1*',comp.Name,'.PD1_aPD1_PD1']);
    set (R,'Notes',['binding and unbinding of PD1 to PD1-Nivo on ',Tname,' surface in synapse']); 
% PDL1-aPDL1 
R = addreaction(model,[comp.Name,'.PDL1 <-> ',comp.Name,'.PDL1_aPDL1']);
    set (R,'ReactionRate', ['2*kon_PDL1_aPDL1*(',comp.Name,'.PDL1 *aPDL1/',gam.aPDL1,') -  koff_PDL1_aPDL1*',comp.Name,'.PDL1_aPDL1']);
    set (R,'Notes',['binding and unbinding of PDL1 to anit-PDL1 on ',Cname,' surface in synapse']); 
% PDL1-aPDL1-PDL1 
R = addreaction(model,[comp.Name,'.PDL1_aPDL1 + ',comp.Name,'.PDL1 <-> ',comp.Name,'.PDL1_aPDL1_PDL1']);
    set (R, 'ReactionRate', ['Chi_PDL1*kon_PDL1_aPDL1*(',comp.Name,'.PDL1 * ',comp.Name,'.PDL1_aPDL1) -  2*koff_PDL1_aPDL1*',comp.Name,'.PDL1_aPDL1_PDL1']);
    set (R, 'Notes'       , ['binding and unbinding of PDL1 to PDL1-anti-PDL1 on ',Cname,' surface in synapse']);   

% Update PD1 Hill Function
addrule(model,['H_PD1_',Cname,' = ((',comp.Name,'.PD1_PDL1+',comp.Name,'.PD1_PDL2)/PD1_50)^n_PD1/(((',...
    comp.Name,'.PD1_PDL1+',comp.Name,'.PD1_PDL2)/PD1_50)^n_PD1 + 1)'],'repeatedAssignment');