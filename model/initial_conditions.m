% Initial Conditions Module
%
% Generates initial conditions
%
% Inputs: model           -- pre-initialized SimBiology model object 
%         varargin        -- variant object      
%        
% Outputs: model   -- initialized SimBiology model
%          success -- shows if it reached the initial tumor diameter
%          simData -- initial condition data
%
% future release: name/value pair inputs, IC output time resolution
%
% Created: Jan 28, 2019 (Richard Sove)
% Last Modified: Oct 31, 2019 (MJ)

function [model_out,success,simData,idx,w] = initial_conditions(model_in,varargin)

% Copy Model Object
model = copyobj(model_in);
  
% Optional Inputs
p = inputParser;
addParameter(p,'Variant','');
addParameter(p,'Debug','');
parse(p,varargin{:});
Variant = p.Results.Variant;
Debug = p.Results.Debug;
if (isempty(Debug))
    Debug = false;
end

% Get Initial Tumour Diameter
if (~isempty(Variant))
    variantWithInitialTumourDiameter = false;
    % Check if 'Variant Object' Contains Initial Tumour Diameter
    for i = 1:length(Variant.content)
        variantParam = Variant.content{i};
        if strcmp(variantParam(2),'initial_tumour_diameter')
            % Get Initial Tumour Diameter from 'variant'
            D_Ti = sbioselect(model,'Name','initial_tumour_diameter');
            tumour_diameter.Value = variantParam{4};
            tumour_diameter.Units = D_Ti.ValueUnits;
            tumour_diameter.Notes = D_Ti.Notes;
            variantWithInitialTumourDiameter = true;
        end
    end
    if ~variantWithInitialTumourDiameter
        % Get Initial Tumour Diameter from 'model'
        D_Ti = sbioselect(model,'Name','initial_tumour_diameter');
        tumour_diameter.Value = D_Ti.Value;
        tumour_diameter.Units = D_Ti.ValueUnits;
        tumour_diameter.Notes = D_Ti.Notes;
    end
else
    % Get Initial Tumour Diameter from 'model'
    D_Ti = sbioselect(model,'Name','initial_tumour_diameter');
    tumour_diameter.Value = D_Ti.Value;
    tumour_diameter.Units = D_Ti.ValueUnits;
    tumour_diameter.Notes = D_Ti.Notes;
end
% Calculate Target Tumour Volume
tumour_volume = 4/3*pi*(tumour_diameter/2)^3;

% Get User-Defined Confiuration Settings
config = getconfigset(model);
user_output_times = config.SolverOptions.OutputTimes; 

% Reset Output Times for IC Simulation (assumes time in days)
set(config,'StopTime',8000);
set(config.SolverOptions,'OutputTimes',[]);

% IC Simulation 
if ((Debug)&&(~isempty(Variant))) 
    simData = sbiosimulate(model,Variant);
elseif ((Debug)&&(isempty(Variant)))
    simData = sbiosimulate(model);
elseif ((~Debug)&&(isempty(Variant)))
    try
        simData = sbiosimulate(model);
    catch
        disp('There was an error while running the initialization - use Debug mode for more information'); 
        model_out = copyobj(model);
        success = false;
        simData = [];
        return
    end   
elseif ((~Debug)&&(~isempty(Variant)))
    try
        simData = sbiosimulate(model,Variant);
    catch
        disp('There was an error while running the initialization - use Debug mode for more information'); 
        model_out = copyobj(model);
        success = false;
        simData = [];
        return
    end
end

% Get Time to Reach Target Tumour Size
% Get Tumour Volume from Simulation
V_T = selectbyname(simData,'V_T','Format','struct');
V_T_value = V_T.CompartmentData.Data;
V_T_units = V_T.CompartmentData.Units{:};
% Convert Units of Target Volume to Units in Simulation
zero.Value = 0; % define 'zero' parameter object with same units as simulation
zero.Units = V_T_units;
tumour_volume = zero + tumour_volume; % operator converts to units of first input 
target_V_T = tumour_volume.Value; % tumour volume in same units as simulation
% Get Time of Target Tumour Size
idx = find((target_V_T-V_T_value)<0,1); % difference should cross zero at target time

% Re-Initialize Model with New ICs
if isempty(idx) % tumour did not reach target size
    success = false;
    w = [];
else
    success = true;
    % Reset User Defined Configuration Settings
    set(config,'StopTime',user_output_times(end));
    set(config.SolverOptions,'OutputTimes',user_output_times);
    
    % Calculate Interpolation Weight
    w = (target_V_T-V_T_value(idx-1))/(V_T_value(idx)-V_T_value(idx-1));
    
    % Set New ICs for Species
    for i = 1:length(model.Species)
        % It makes sure not to initialize the checkpoints because of
        % initial assignments
        if (isempty(strfind(simData.DataInfo{i}.Compartment,'syn' )))&&(isempty(strfind(simData.DataNames{i},'OX40')))
            model.Species(i).InitialAmount = w*simData.Data(idx,i)+(1-w)*simData.Data(idx-1,i);
        end
    end    
    % Set New ICs for varying Parameters
    for i = 1:length(model.Parameters)
        if ~(isempty(find(strcmp(simData.DataNames,model.Parameters(i).Name ))))
            k = find(strcmp(simData.DataNames, model.Parameters(i).Name));
            model.Parameters(i).Value =  simData.Data(idx,k);
        end
    end
    % Set New ICs for size-changing Compartment
    for i = 1:length(model.Compartments)
        if ~(isempty(find(strcmp(simData.DataNames,model.Compartments(i).Name ))))
            k = find(strcmp(simData.DataNames, model.Compartments(i).Name));
            model.Compartments(i).Capacity = simData.Data(idx,k);
        end
    end

end

model_out = copyobj(model);