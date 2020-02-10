% Function to set Initial Conditions
%
% Inputs: model_in -- SimBiology model object 
%         ICs      -- vector of IC values 
%        
% Outputs: model_out -- SimBiology model object with ICs
%
% Created: Mar 27, 2019 (Richard J. Sové)
% Last Modified: Mar 27, 2019 (RJS)

function model_out = set_ICs(model_in,ICs)
  
model_out = copyobj(model_in);  

for i = 1:length(model_in.Species)
  model_out.Species(i).InitialAmount = ICs(i);
end