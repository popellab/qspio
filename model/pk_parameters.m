% PK parameter database
%
% Inputs: name -- string containing drug name <optional>:
%                 - nivolumab
%                 - durvalumab *not calibrated*
%                 - ipilimumab *not calibrated*
%
% Output: params_out -- object containing model parameters
%                       - q_P--rate of diffusive transport C<->P
%                       - q_T--rate of diffusive transport C<->T
%                       - q_LN--rate of diffusive transport C<->LN
%                       - q_LD--rate of convective transport T->LN 
%                       - k_cl--clearence rate from central
%                       - gamma_C--volume fraction in C
%                       - gamma_P--volume fraction in P
%                       - gamma_T--volume fraction in T
%                       - gamma_LN--volume fraction in LN
%                       - type--string describing drug type
%
% Created: Oct 25, 2018 (Richard Sové)
% Last Modified: Feb 18, 2018 (RJS) 

function params_out = pk_parameters(name)

% Optional Inputs
if ~exist('name','var')
    name = '';
end

% Drug Database
% - Transport Rates (q_P,q_T,q_LN,q_LD)
% - Clearence Rate (k_cl)
% - Volume Fractions (gamma_C,gamma_P,gamma_T,gamma_LN)
switch (name)
    case 'nivolumab'
        % Central to Peripheral
        q_P.Value = 1.73655e-7;
        q_P.Units = '1/second';
        q_P.Notes = '(Finley 2012)';
        % Central to Tumour
        q_T.Value = 8.52e-6;
        q_T.Units = '1/second';
        q_T.Notes = '(Finley 2012)';
        % Central to LN
        q_LN.Value = 1.73655e-7;
        q_LN.Units = '1/second';
        q_LN.Notes = '(Padera 2017)';
        % Tumour to LN
        q_LD.Value = 0.0015;
        q_LD.Units = '1/minute';
        q_LD.Notes = '(Zhu Jain 1996)';
        % Clearence
        k_cl.Value = 7.1813e-07;
        k_cl.Units = '1/second';
        k_cl.Notes = '(fitted?)';
        % Volume Fractions
        gamma_C.Value = 0.77378;
        gamma_C.Notes = '(estimated?)';
        gamma_P.Value = 0.058847;
        gamma_P.Notes = '(estimated?)';
        gamma_T.Value = 0.718;
        gamma_T.Notes = '(Coughlin 2010)';
        gamma_LN.Value = 0.1;
        gamma_LN.Notes = '(Jafarnejad 2017)';
        params_out.type = 'aPD1';
    case 'cabozantinib'
        % Central to Peripheral
        q_P.Value = 4.9370e-05;
        q_P.Units = '1/second';
        q_P.Notes = '(estimated)';
        % Central to Tumour
        q_T.Value = 2.4222e-03;
        q_T.Units = '1/second';
        q_T.Notes = '(estimated)';
        % Central to LN
        q_LN.Value = 4.9370e-05;
        q_LN.Units = '1/second';
        q_LN.Notes = '(estimated)';
        % Tumour to LN
        q_LD.Value = 7.1075e-03;
        q_LD.Units = '1/second';
        q_LD.Notes = '(estimated)';
        % Clearence
        k_cl.Value = 6.2e-03;
        k_cl.Units = '1/hour';
        k_cl.Notes = '(Lacey 2017)';
        % Volume Fractions
        gamma_C.Value = 0.77378;
        gamma_C.Notes = '(estimated)';
        gamma_P.Value = 0.058847;
        gamma_P.Notes = '(estimated)';
        gamma_T.Value = 0.718;
        gamma_T.Notes = '(Coughlin 2010)';
        gamma_LN.Value = 0.1;
        gamma_LN.Notes = '(Jafarnejad 2017)';
        params_out.type = 'TKI';
    otherwise
        disp('drug requested not available, using parameters for nivolumab')
        % Central to Peripheral
        q_P.Value = 1.73655e-7;
        q_P.Units = '1/second';
        q_P.Notes = '(Finley 2012)';
        % Central to Tumour
        q_T.Value = 8.52e-6;
        q_T.Units = '1/second';
        q_T.Notes = '(Finley 2012)';
        % Central to LN
        q_LN.Value = 1.73655e-7;
        q_LN.Units = '1/second';
        q_LN.Notes = '(Padera 2017)';
        % Tumour to LN
        q_LD.Value = 0.0015;
        q_LD.Units = '1/minute';
        q_LD.Notes = '(Zhu Jain 1996)';
        % Clearence
        k_cl.Value = 7.1813e-07;
        k_cl.Units = '1/second';
        k_cl.Notes = '(fitted?)';
        % Volume Fractions
        gamma_C.Value = 0.77378;
        gamma_C.Notes = '(estimated?)';
        gamma_P.Value = 0.058847;
        gamma_P.Notes = '(estimated?)';
        gamma_T.Value = 0.718;
        gamma_T.Notes = '(Coughlin 2010)';
        gamma_LN.Value = 0.1;
        gamma_LN.Notes = '(Jafarnejad 2017)';
        params_out.type = 'aPD1';
end

% Diffusive Transport: C<->P
params_out.q_P = q_P;
% Diffusive Transport: C<->T 
params_out.q_T = q_T;
% Diffusive Transport: C<->LN 
params_out.q_LN = q_LN;
% Convective Transport: T->LN || LN->C 
params_out.q_LD = q_LD;
% Clearence from Central
params_out.k_cl = k_cl;
% Volume Fractions
gamma_C.Units = 'dimensionless';
gamma_P.Units = 'dimensionless';
gamma_T.Units = 'dimensionless';
gamma_LN.Units = 'dimensionless';
params_out.gamma_C = gamma_C; 
params_out.gamma_P = gamma_P; 
params_out.gamma_T = gamma_T;
params_out.gamma_LN = gamma_LN;
