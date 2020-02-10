% Example 2 - script for setting up the model of example 2 from Sove et al 2020

% Model Settings
model_name = 'Immune Oncology Model';
start_time = 0.0; % solution start time [days]
time_step = 1; % solution output interval [days] (0.01 days ~ 15 mins)
end_time = 3000; % solution end time [days]

% Load Model Parameters 
parameter_filename = 'parameters/example2_parameters.json';
params = load_parameters(parameter_filename);
params2 = params;
params2.k_C_growth = 0.5*params.k_C_growth;
params2.k_C_growth.Notes = '';
params2.k_C_Tcell = 0.1*params.k_C_Tcell;
params2.k_C_Tcell.Notes = '';

% Model Initialization
time = start_time:time_step:end_time;
model = simbio_init(model_name,time,params);

% Add Modules
% Cancer Module
model = cancer_module(model,'C1',params);
model = cancer_module(model,'C2',params2);
% T Cell Module
model = Tcell_module(model,'1',params,{'C1'});
model = Tcell_module(model,'2',params2,{'C1','C2'});
% PD1 Immune Checkpoint Module
model = PD1_module(model,params,'T1','C1','drugName','nivolumab'); 

% Dose Schedule
dose_schedule = schedule_dosing({'nivolumab'});