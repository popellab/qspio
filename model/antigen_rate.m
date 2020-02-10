% Function to generate object with default model parameters
%
% Inputs: model   -- SimBiology model object
%         ID      -- T cell-antigen ID number
%         cancers -- cell array of string containing cancer names
%         P_C     -- parameter array (same length as 'cancers') of antigen concentrations per cancer cell
%         nTcells -- number of cytotoxic T cell clones in model
%
% Output: rate    -- string containing rate of antigen release from dying cancer cells
%
% Created: Nov 27, 2018 (Richard Sové)
% Last Modified: March 6, 2019 (RJS)

function rate = antigen_rate(model,ID,cancers,P_C,nTcells)

% Parse Inputs
P = ['P' ID];
C = cancers;

% Add Antigen Concentrations to Model
for i = 1:length(P_C)
    param = addparameter(model,[P '_' C{i}],P_C(i).Value,'ValueUnits',P_C(i).Units);
        set(param,'Notes',['Concentration of ' P ' in ' C{i}]);
end

rate = '';
for i = 1:length(C)
    kT = '';
    for j = 1:nTcells
        % kT = kT + k_Tj*Tj 
        kT = sprintf('%s+k_C_T%d*V_T.T%d',kT,j,j);
    end
    % rate = rate + P_Ci*(k_Ci_death+k_Ci_therapy+kT/(C_total+T_total)*(1-H_PD1))*Ci
    rate = sprintf('%s+%s_%s*(k_%s_death+k_%s_therapy+(%s)/(C_total+T_total+cell)*(1-H_PD1_C1))*%s',...
          rate,P,C{i},C{i},C{i},kT,C{i});
end
rate = ['n_T' ID '_clones*(' rate ')*V_T'];