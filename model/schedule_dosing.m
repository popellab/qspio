% Dose Schedule
%
% Generates Donsing Schedule 
%
% Inputs: drugName -- drug name or character array of drug names
%         varargin -- Name-Value Pairs
%                     - drugName_dose (mg/kg)
%                     - drugName_schedule [start,interval,repeat]
%                     - patientWeight (kg)
%
%       vaild drugName values: nivolumab, durvalumab, ipilimumab
%        
% Outputs: dosing -- SimBiology model object with new antigen module
%
% Created: Feb 19, 2019 (Richard Sové)
% Last Modified: Feb 20, 2019 (RJS)

function dose_schedule = schedule_dosing(drugName,varargin)

% Check if drugName is cell array 
if (iscell(drugName))
    N = length(drugName);
else
    N = 1;
    drugName = {drugName};
end

% Optional Inputs
in = inputParser;
% Nivolumab
addParameter(in,'nivolumab_dose',3); % 3 mg/kg every two weeks
addParameter(in,'nivolumab_schedule',[0,14,30]); 
% Durvalumab
addParameter(in,'durvalumab_dose',10); % 10 mg/kg every two weeks
addParameter(in,'durvalumab_schedule',[0,14,30]); 
% Ipilimumab
addParameter(in,'ipilimumab_dose',1); % 1 mg/kg every three weeks
addParameter(in,'ipilimumab_schedule',[0,21,30]); 
% Cabozantinib
addParameter(in,'cabozantinib_dose',60); % 60 mg daily
addParameter(in,'cabozantinib_schedule',[0,1,500]); 
% Patient Weight
addParameter(in,'patientWeight',70); % 70 kg (standard man)
% Parse Inputs
parse(in,varargin{:});
% Nivolumab
dose_nivo = in.Results.nivolumab_dose;
schedule_nivo = in.Results.nivolumab_schedule;
% Durvalumab
dose_durv = in.Results.durvalumab_dose;
schedule_durv = in.Results.durvalumab_schedule;
% Ipilimumab
dose_ipil = in.Results.ipilimumab_dose;
schedule_ipil = in.Results.ipilimumab_schedule;
% Cabozantinib
dose_cabo = in.Results.cabozantinib_dose;
schedule_cabo = in.Results.cabozantinib_schedule;
% Patient Weight
patient_weight = in.Results.patientWeight;

% Nivolumab
MW_nivo = 1.436E8; % milligrams per mole
doseObj_nivo = sbiodose('nivolumab','Amount',patient_weight*dose_nivo/MW_nivo,'AmountUnits','mole','TargetName','V_C.nivolumab');
doseObj_nivo.StartTime = schedule_nivo(1);
doseObj_nivo.Interval = schedule_nivo(2);
doseObj_nivo.TimeUnits = 'day';
doseObj_nivo.RepeatCount = schedule_nivo(3);
doseObj_nivo.Active = true;

% Durvalumab 
MW_durv = 1.436E8; % milligrams per mole 
doseObj_durv = sbiodose('durvalumab','Amount',patient_weight*dose_durv/MW_durv,'AmountUnits','mole','TargetName','V_C.durvalumab');
doseObj_durv.StartTime = schedule_durv(1);
doseObj_durv.Interval = schedule_durv(2);
doseObj_durv.TimeUnits = 'day';
doseObj_durv.RepeatCount = schedule_durv(3);
doseObj_durv.Active = true;

% Ipilimumab 
MW_ipil = 1.486349E8; % milligrams per mole
doseObj_ipil = sbiodose('ipilulumab','Amount',patient_weight*dose_ipil/MW_ipil,'AmountUnits','mole','TargetName','V_C.ipililumab');
doseObj_ipil.StartTime = schedule_ipil(1);
doseObj_ipil.Interval = schedule_ipil(2);
doseObj_ipil.TimeUnits = 'day';
doseObj_ipil.RepeatCount = schedule_ipil(3);
doseObj_ipil.Active = true;

% Cabozantinib
MW_cabo = 5.051E5; % milligrams per mole
doseObj_cabo = sbiodose('cabozantinib','Amount',dose_cabo/MW_cabo,'AmountUnits','mole','TargetName','V_C.cabozantinib');
doseObj_cabo.StartTime = schedule_cabo(1);
doseObj_cabo.Interval = schedule_cabo(2);
doseObj_cabo.TimeUnits = 'day';
doseObj_cabo.RepeatCount = schedule_cabo(3);
doseObj_cabo.Active = true;

% Dose Schedule Array
dose_schedule(N) = sbiodose('empty'); % preallocate array
for i = 1:N
    switch drugName{i}
        case 'nivolumab'
            dose_schedule(i) = doseObj_nivo;
        case 'durvalumab'
            dose_schedule(i) = doseObj_durv;
        case 'ipilimumab'
            dose_schedule(i) = doseObj_ipil;
        case 'cabozantinib'
            dose_schedule(i) = doseObj_cabo;
        otherwise
            error('No match for drug name');
    end
end