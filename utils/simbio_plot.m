% SimBiology Plot Generator
%
% Plots SimBiology Data
%
% Inputs: simData - SimBiology output data
%         name    - string with species name
%
% Optional Name-Value Pairs
%         legend -- string containing legend entry
%         compartment -- string containing compartment name
%         axisLabel -- srting containing y-axis label, units are added by default
%         normalizeBy -- string or number containing normalization factor
%         
% Created: Nov 12, 2018 (Richard Sové)
% Last Modified: Oct 23, 2019 (RJS) 

function h = simbio_plot(simData,name,varargin)

% Optional Inputs
p = inputParser;
addParameter(p,'legend','');
addParameter(p,'compartment','');
addParameter(p,'axisLabel','');
addParameter(p,'normalizeBy','');
addParameter(p,'linespec','');
parse(p,varargin{:});
legend = p.Results.legend;
comp = p.Results.compartment;
axis_label = p.Results.axisLabel;
normFac_str = p.Results.normalizeBy;
linespec = p.Results.linespec;

% Determine Normalization
if ~strcmp(normFac_str,'')
    % Normalization
    if isnumeric(normFac_str)
        % Not a String
        normFac = normFac_str;
    else
        % String
        if ~strcmp(normFac_str(1:2),'1/')
            % Normal Normalization
            norm_name = normFac_str;
            [~,normFac,~] = selectbyname(simData,norm_name);
        else
            % Reciprocal Normalization
            norm_name = normFac_str(3:end);
            [~,temp,~] = selectbyname(simData,norm_name);
            normFac = 1./temp;
        end        
    end
else
    % No Normalization
    normFac = 1;
end

% Set Figure Defaults
figure_defaults;

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
    h = plot(simData.time,simData.data(:,idx)./normFac,linespec,'DisplayName',legend); % plot data
    units = simData.DataInfo{idx}.Units; % get units for y-axis label
    if (strcmp(units,'dimensionless')||strcmp(units,'cell'))
        units = ''; % do not display units for dimensionless quantities
    else
        units = symunit2str(str2symunit(units,'SimBiology'),'Simulink'); % convert to symbols
        units = [' ($' units '$)']; % enclose units in brackets
    end
    ystring = [axis_label units];
    xlabel('Time (days)'); ylabel(ystring);
catch
    h = [];
    disp(['Species ' name ' not found']);
end