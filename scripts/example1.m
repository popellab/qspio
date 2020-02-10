% Example 1 - script for setting up the model of example 1 from Sove et al 2020

% Model Settings
model_name = 'Immune Oncology Model';
start_time = 0.0; % solution start time [days]
time_step = 0.01; % solution output interval [days] (0.01 days ~ 15 mins)
end_time = 400; % solution end time [days]

% Load Model Parameters 
parameter_filename = 'parameters/example1_parameters.json';
params = load_parameters(parameter_filename);

% Model Initialization
time = start_time:time_step:end_time;
model = simbio_init(model_name,time,params);

% Add Modules
% Cancer Module
model = cancer_module(model,'C1',params);
% T Cell Module
model = Tcell_module(model,'1',params,{'C1'});
% T Reg Module
model = Treg_module(model,params);
% APC Module
model = APC_module(model,params);
% Antigen Module
self_antigen = create_antigen({'C1'},1e-8,'antigenID',0);
neo_antigen = create_antigen({'C1'},1e-8,'antigenID',1);
model = antigen_module(model,'0',params,self_antigen);
model = antigen_module(model,'1',params,neo_antigen);
% PD1 Immune Checkpoint Module
model = PD1_module(model,params,'T1','C1','drugName','nivolumab'); 

% Dose Schedule
dose_schedule = schedule_dosing({'nivolumab'});