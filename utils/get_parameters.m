% Function to get model parameters from simbiology object
%
% Inputs: model -- SimBiology model object 
%        
% Outputs: params -- parameter object
%
% Created: Mar 27, 2019 (Richard J. Sové)
% Last Modified: Mar 27, 2019 (RJS)

function params = get_parameters(model)

% Loop Through Model Parameters
for i = 1:length(model.Parameters)
    if (model.Parameters(i).ConstantValue) 
        if ~((model.Parameters(i).Name(1)=='H')||(model.Parameters(i).Name(1)=='a')||...
                contains(model.Parameters(i).Name,'tot')||strcmp(model.Parameters(i).Name,'cell')||...
                strcmp(model.Parameters(i).Name,'day')||strcmp(model.Parameters(i).Name,'Treg_'))
            params.(model.Parameters(i).Name).val = model.Parameters(i).Value;
        end
    end
end