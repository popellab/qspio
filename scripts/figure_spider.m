% Spider Plot Figure - script to generate tumour volumes of LHS sampled parameters
% (Figure 4 in Sove et al 2020)
%
% Set N_lhs to change the number of "patients" simulated
%
% Note: Each run will produce a different figure since the LHS is randomly
% generated. 

% Setup Model
example1;

% Number of Simulations
N_lhs = 5;

% Get Model Parameters
all_params = get_parameters(model);

% Set Parameter Subspace
input_params.k_C1_growth.min = 1e-3; input_params.k_C1_growth.max = 0.05;
input_params.C_max.min = 1e11; input_params.C_max.max = 1e13;
input_params.initial_tumour_diameter.min = 0.5; input_params.initial_tumour_diameter.max = 5;
input_params.n_T1_clones.min = 1; input_params.n_T1_clones.max = 1000;
input_params.Q_nT1_thym.min = 1e8; input_params.Q_nT1_thym.max = 1e10;
input_params.k_C_T1.min = 0.5; input_params.k_C_T1.max = 50;
input_params.k_T1_death.min = 0.01; input_params.k_T1_death.max = 10;
input_params.q_T1_LN_out.min = 0.1; input_params.q_T1_LN_out.max = 10;
input_params.k_cell_clear.min = 0.001; input_params.k_cell_clear.max = 0.04;
input_params.PD1_50.min = 1; input_params.PD1_50.max = 50;

% Run LHS of Parameter Subspace
[rho,in,out,fig] = run_lhs(model,dose_schedule,input_params,N_lhs);
set(gca,'yscale','log'); ylabel('Tumour Volume ($\mu$L)');

% Get Data
axObj = fig.Children;
dataObj = axObj.Children;
t = dataObj.XData;
V_T = zeros(length(dataObj),length(t));
for i = 1:length(dataObj)
    V_T(i,:) = dataObj(i).YData;
end
close(fig);

% Subsample Time
t = t(1:100:end);
V_T = V_T(:,1:100:end);

% Get Diameter
D_T = 2*(3*V_T/4/pi).^(1/3);

% Calculate % Change
perc_change = zeros(size(V_T));
for i = 1:length(dataObj)
    perc_change(i,:) = (D_T(i,:)-D_T(i,1))/D_T(i,1)*100;
end

% Calculate RECIST
PR = sum(perc_change<-30)/size(perc_change,1);
PD = sum(perc_change>20)/size(perc_change,1);
SD = 1-PR-PD;

% Waterfall Plot
t_ind = find(t==60);
perc_change_ordered = sort(perc_change(:,t_ind));

% Figure
figure;
subplot(2,2,1); plot(t,V_T/1e3); set(gca,'YScale','log');
xlabel('Time (days)'); ylabel('Volume (mL)');
subplot(2,2,2); plot(t,perc_change); ylim([-100 150]); xlim([0 500]);
hold on; plot([0 500],[-30 -30],'k--'); plot([0 500],[20 20],'k--');
xlabel('Time (days)'); ylabel('\% Change in Diameter');
subplot(2,2,3);
patch([t t(end) 0],[PR 0 0],[204/255 229/255 1]);
patch([t t(end) t(end:-1:1)],[PR+SD PR(end) PR(end:-1:1)],[224/255 224/255 224/255]);
patch([0 t(end) t(end:-1:1)],[1 1 PR(end:-1:1)+SD(end:-1:1)],[1 229/255 204/255]);
ylim([0 1]);
xlabel('Time (days)'); ylabel('RECIST');
subplot(2,2,4); bar(perc_change_ordered(end:-1:1),'BarWidth',1); ylim([-100 100]);
ylabel('\% Change in Diameter'); xlabel('Simulation Number'); xlim([0 N_lhs-1]);