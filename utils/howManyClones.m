% Function to find the number of clones
%
% Inputs: dataObj    -- Object containing SimBiology model object or model outputs 
%        
% Outputs: numClones -- number of T cell clones in the simulation
%
% Created: Jan 23, 2019 (Mohammad Jafarnejad)
% Last Modified: Feb 08, 2019 (RJS)

function numClones = howManyClones(dataObj)

% Determine Object Type
try % model object
    data_names = cell(length(dataObj.Species),1);
    for i = 1:length(dataObj.Species)
        data_names{i} = dataObj.Species(i).Name;
    end
catch
    data_names = dataObj.DataNames;
end 

% Finds the number of clones in the simulation by parsing 'Ti's
Tclones = zeros(size(data_names));
for i = 1:length(data_names)
    temp = sscanf(data_names{i},'T%d');
    if ~isempty(temp)
        Tclones(i) = temp;
    else
        Tclones(i) = NaN;
    end
end
numClones = nanmax(Tclones);