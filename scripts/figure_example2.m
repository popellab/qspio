% Example 2 Figure - script to plot example 2 (Figure 5 in Sove et al 2020)

% Setup the Model for Example 1
example2;

% Run Simulation
[model,success] = initial_conditions(model);
simData_noTreat = sbiosimulate(model);
simData = sbiosimulate(model,dose_schedule);

figure; 
% Tumour Volume
subplot(2,1,1); hold on;
simbio_plot(simData_noTreat,'V_T','normalizeBy',1e3);
simbio_plot(simData,'V_T','normalizeBy',1e3);
plot([420 420],[0 100],'k:');
ylim([0 100]); ylabel('Volume (mL)'); axis square; legend('Control','Anti-PD-1');
% Cell Count
subplot(2,1,2); hold on;
simbio_plot(simData_noTreat,'C1');
simbio_plot(simData,'C1');  set(gca,'ColorOrderIndex',1);
simbio_plot(simData_noTreat,'C2','linespec','--');
simbio_plot(simData,'C2','linespec','--');
plot([420 420],[1e0 1e15],'k:');
ylabel('Cell Count'); set(gca,'YScale','log'); axis square;
legend('$C_1$ (Control)','$C_1$ (Anti-PD-1)','$C_2$ (Control)','$C_2$ (Anti-PD-1)');