% Function to check the Hill functions connecting the modules and to define
% them  if they have not been already defined
%
% Inputs: model      -- simbiology model object
%         params     -- object containing model parameters
%
% Outputs: model -- simbiology model object containing the hill parameters
%
% Created: Jun 12, 2019 (Mohammad Jafarnejad)
% Last Modified: Jun 12, 2019 (MJ)

function model_out = simbio_connections(model,params)

% Set Default Hill Function for PD1 Checkpoint
addparameter(model,'H_PD1_C1',0.90,'ValueUnits','dimensionless','ConstantValue',false);
addparameter(model,'H_PD1_APC',0.90,'ValueUnits','dimensionless','ConstantValue',false);

% Set Default Hill Function for CTLA4 Checkpoint
addparameter(model,'H_CD28_C1',1.0,'ValueUnits','dimensionless','ConstantValue',false);
addparameter(model,'H_CD28_APC',1.0,'ValueUnits','dimensionless','ConstantValue',false);

model_out = copyobj(model);