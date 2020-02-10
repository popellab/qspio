% Function to load model parameters from file
%
% Inputs: filename   -- string containing path and filename to parameters file
%
% Output: params_out -- structure containing model parameters [Value,Units,Notes]
%
% Created: Aug 15, 2019 (Richard Sové)
% Last Modified: Aug 15, 2019 (RJS)

function params_out = load_parameters(filename)

% Add 'cell' unit to SimBiology and Symbolic Toolboxes
% SimBiology
if (isempty(sbioshowunits('cell')))
    cell_unit = sbiounit('cell','molecule');
    sbioaddtolibrary(cell_unit);
end
% Symbolic Unit
u = symunit;
try u.cell;
catch
    newUnit('cell',u.molecule);
end

% Load text from .json file
json = jsondecode(fileread(filename));

% Parse Text
for i = 1:length(json)
    params_out.(json{i}.name).Value = json{i}.value;
    params_out.(json{i}.name).Units = json{i}.units;
    params_out.(json{i}.name).Notes = json{i}.description;
end

% Calculate Derived Parameters
for i = 1:length(json)
    if (strcmp(json{i}.source,'derived'))
        n = length(json{i}.derived_from);
        for k = 1:n
            p(k) = params_out.(json{i}.derived_from{k});
        end
        params_out.(json{i}.name) = eval(json{i}.expression);
        params_out.(json{i}.name).Notes = json{i}.description;
    end
end