% Function to initialize simbiology object
%
% Inputs: model_name -- string containing model name
%         time       -- vector of output times (in days)
%         params     -- object containing compartment parameters
%                       - V_C--central compartment volume
%                       - V_P--peripheral compartment volume
%                       - V_Tmin--cancer-free tumour compartment volume
%                       - V_LN--lymph node compartment volume
%                       - vol_cell--cancer cell volume
%                       - vol_Tcell--T cell volume
%                       - k_cell_clear--dead cell clearance rate
%         varargin   -- Name-Value pair
%                       - absoluteTolerance -- absolute tolerance of solver
%                       - relativeTolerance -- relative tolerance of solver
%
% Outputs: model -- simbiology model object containing user-specified
%                   configuration and four compartments: C,P,T,LN
%
% Created: Oct 22, 2018 (Richard Sové)
% Last Modified: Oct 01, 2019 (RJS)

function model = simbio_init(model_name,time,params,varargin)

% Parse varargin
p = inputParser;
addParameter(p,'absoluteTolerance',1e-12);
addParameter(p,'relativeTolerance',1e-6);
parse(p,varargin{:});
tol_abs = p.Results.absoluteTolerance;
tol_rel = p.Results.relativeTolerance;

% Model Object
model = sbiomodel(model_name);
config = getconfigset(model);
options = get(config,'CompileOptions');
set(options,'UnitConversion',true);
set(config,'TimeUnits','day');
set(config,'SolverType','ode15s'); 
set(config.SolverOptions,'OutputTimes',time);
set(config,'StopTime',time(end)); 
set(config.SolverOptions,'AbsoluteTolerance',tol_abs);
set(config.SolverOptions,'RelativeTolerance',tol_rel);

% Setup Compartments
comp_C = addcompartment(model,'V_C',params.V_C.Value,'CapacityUnits',params.V_C.Units); 
    set(comp_C,'Notes',['Central compartment (C) ' params.V_C.Notes]);
comp_P = addcompartment(model,'V_P',params.V_P.Value,'CapacityUnits',params.V_P.Units); 
    set (comp_P,'Notes',['Peripheral compartment (P) ' params.V_P.Notes]);
comp_T = addcompartment(model,'V_T',params.V_Tmin.Value,'CapacityUnits',params.V_Tmin.Units,'ConstantCapacity',false); 
    set (comp_T,'Notes',['Tumor compartment (T) ' params.V_Tmin.Notes]);
comp_LN = addcompartment(model,'V_LN',params.V_LN.Value,'CapacityUnits',params.V_LN.Units); 
    set(comp_LN,'Notes',['Lymph node (LN) compartment volume ' params.V_LN.Notes]);
    
% Dead Cells
P = addparameter(model,'k_cell_clear',params.k_cell_clear.Value,'ValueUnits',params.k_cell_clear.Units);
    set(P,'Notes',['Rate of dead cell clearance from tumour compartment ' params.k_cell_clear.Notes]);
S = addspecies(comp_T,'C_x',0,'InitialAmountUnits','cell');
    set(S,'Notes','Dead cancer cells in the tumour compartment');
S = addspecies(comp_T,'T_exh',0,'InitialAmountUnits','cell');
    set(S,'Notes','Exhausted T cells');
R = addreaction(model,'V_T.C_x -> null');
    set(R,'ReactionRate','k_cell_clear*V_T.C_x');
R = addreaction(model,'V_T.T_exh -> null');
    set(R,'ReactionRate','k_cell_clear*V_T.T_exh');

% Define Cell and Time    
addparameter(model,'cell',1.0,'ValueUnits','cell');
addparameter(model,'day',1.0,'ValueUnits','day');

% Set Tumour Volume (Rule 1) 
vol_cell = addparameter(model,'vol_cell',params.vol_cell.Value,'ValueUnits',params.vol_cell.Units);
    set(vol_cell,'Notes',['Average volume of cancer cell ' params.vol_cell.Notes]);
vol_Tcell = addparameter(model,'vol_Tcell',params.vol_Tcell.Value,'ValueUnits',params.vol_Tcell.Units);
    set(vol_Tcell,'Notes',['Average volume of T cells ' params.vol_Tcell.Notes]);
V_Tmin = addparameter(model,'V_Tmin',params.V_Tmin.Value,'ValueUnits',params.V_Tmin.Units);
    set(V_Tmin,'Notes',['Cancer-Free Tumour compartment volume' params.V_Tmin.Notes]);
addrule(model,'V_T = V_Tmin+vol_cell*C_x+vol_Tcell*T_exh','repeatedAssignment');

% Set Total Number of Cancer Cells (Rule 2)
addparameter(model,'C_total',0,'ValueUnits','cell','ConstantValue',false);
addrule(model,'C_total = 0*cell','repeatedAssignment');

% Set Total Number of T Cells in Tumour (Rule 3)
addparameter(model,'T_total',0,'ValueUnits','cell','ConstantValue',false);
addrule(model,'T_total = 0*cell','repeatedAssignment');

% Set Total Number of T Cells in LN (Rule 4) 
addparameter(model,'T_total_LN',0,'ValueUnits','cell','ConstantValue',false);
addrule(model,'T_total_LN = 0*cell','repeatedAssignment');

% Set Total Rate of Cancer Death by T Cells (Rule 5)
addparameter(model,'R_Tcell','ValueUnits','cell/day','ConstantValue',false);
addrule(model,'R_Tcell = 0*cell/day','repeatedAssignment');

% Set Default Number of Tregs
addparameter(model,'Tregs_',0,'ValueUnits','cell','ConstantValue',false);

% Set Default Hill Function for APCs 
addparameter(model,'H_APC',0.5,'ValueUnits','dimensionless','ConstantValue',false);

% Set Default Hill Function for mAPCs
addparameter(model,'H_mAPC',1,'ValueUnits','dimensionless','ConstantValue',false);

% Set Default Hill Function for PD1 Checkpoint
addparameter(model,'H_PD1_C1',0.90,'ValueUnits','dimensionless','ConstantValue',false);
addparameter(model,'H_PD1_APC',0.90,'ValueUnits','dimensionless','ConstantValue',false);

% Set Default Hill Function for CTLA4 Checkpoint
addparameter(model,'H_CD28_C1',0.1,'ValueUnits','dimensionless','ConstantValue',false);
addparameter(model,'H_CD28_APC',0.1,'ValueUnits','dimensionless','ConstantValue',false);