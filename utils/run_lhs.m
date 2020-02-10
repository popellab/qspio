% Function to Run a Latin Hypercube Sampling of Parameter Space
%
% Inputs: model         -- SimBiology model object 
%         dose_schedule -- SimBiology dose object
%         params        -- parameter object created with get_params function
%         N             -- size of LHS (should be 10 the number of params)
%        
% Outputs: rho    -- parameter sensitivities
%          input  -- parameter sensitivity inputs
%          output -- parameter sensitivity outputs
%          fig    -- figure handle to tumour volume of LHS
%
% Created: Mar 27, 2019 (Richard J. Sové)
% Last Modified: Mar 27, 2019 (RJS)

function [rho,input,output,fig] = run_lhs(model,dose_schedule,params,N)

fig = figure; hold on;

% Get Parameter Information
names = fieldnames(params);
P = length(names);

% Create Parameter Sensitivity Inputs/Outputs
input = zeros(N,P);
output = zeros(N,5);

% Create Latin Hypercube Sampling
lhs = lhsdesign(N,P);

% Run Simulations
count = 0;
for k = 1:N
    % Pogress
    disp([num2str(k) '/' num2str(N)]);
    
    % Copy Model Object
    model_cpy = copyobj(model);
    
    % Create Variant
    variantObj = addvariant(model_cpy,['v',num2str(k,'%5.5i')]);
    for i = 1:P
        min_val = params.(names{i}).min;
        max_val = params.(names{i}).max;
        value = min_val+(max_val-min_val)*lhs(k,i);
        input(k,i) = value;
        addcontent(variantObj,{'parameter',names{i},'Value',value});
    end
  
    % Set Initial Conditions
    [model_cpy_ic,success,~] = initial_conditions(model_cpy,'Variant',variantObj);
    
    % Run Simulation
    if success
        simData = sbiosimulate(model_cpy_ic,[],variantObj,dose_schedule);
        [~,V_T,~] = selectbyname(simData,'V_T');
        output(k,1) = (V_T(round(end/2))-V_T(1))/V_T(1)*100;
        output(k,2) = (V_T(round(2*end/3))-V_T(1))/V_T(1)*100;
        output(k,3) = (V_T(round(3*end/4))-V_T(1))/V_T(1)*100;
        output(k,4) = (V_T(round(4*end/5))-V_T(1))/V_T(1)*100;
        output(k,5) = (V_T(end)-V_T(1))/V_T(1)*100;
        h = simbio_plot(simData,'V_T');
        h.UserData = input(k,:);
    else
        count = count+1;
        disp(['Failed to initialize (' num2str(count) ')']);
    end
end

% Calculate Partial Rank Correlation Coefficients
rho = partialcorr(input,output,ones(N,1),'type','Spearman','rows','complete');