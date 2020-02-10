% Function to get Initial Conditions from SimBiology Output Object
%
% Inputs: simData -- SimBiology data object  
%        
% Outputs: ICs -- initial conditions
%
% Created: Mar 27, 2019 (Richard J. Sové)
% Last Modified: Mar 27, 2019 (RJS)

function ICs = get_ICs(simData)
  
if isempty(simData)
    ICs = [];
else 
    N = simData.DataCount.Species;  
    ICs = zeros(1,N);  
    for i = 1:N
      ICs(i) = simData.Data(1,i);
    end
end