% SimBiology Data Retriever
%
% Gets Vector of Data from SimBiology Simulation Output
%
% Inputs: simData - SimBiology output data
%         name    - string with species name
%
% Optional Name-Value Pairs
%         compartmentName -- string containing compartment name <required for species in multiple compartments>
%         
% Created: Feb 08, 2019 (Richard Sové)
% Last Modified: Feb 08, 2019 (RJS) 

function data = simbio_get(simData,name,varargin)

% Optional Inputs
p = inputParser;
addParameter(p,'compartmentName','');
parse(p,varargin{:});
comp = p.Results.compartmentName;

% Get Data Index
for i = 1:size(simData.Data,2)
    if (~isempty(comp))
        try
            comp_match = strcmp(simData.DataInfo{i}.Compartment,comp);
        catch
            break;
        end
    else
        comp_match = true;
    end
    if (strcmp(simData.DataInfo{i}.Name,name)&&comp_match)
        idx = i;
        break;
    end
end

% Plot Data
try
    data = simData.data(:,idx);
catch
    data = [];
    disp(['Species ' name ' not found']);
end