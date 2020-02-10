% Example 1 Figure - script to plot example 1 (Figure 3 in Sove et al 2020)

% Setup the Model for Example 1
example1;

% Run Simulation
[model,success] = initial_conditions(model);
simData_noTreat = sbiosimulate(model);
simData = sbiosimulate(model,dose_schedule);

H = figure;
% Tumour Volume
subplot(3,3,1);
simbio_plot(simData_noTreat,'V_T','legend','Control','normalizeBy',1e3);
hold on; simbio_plot(simData,'V_T','legend','Anti-PD1 Treatment','normalizeBy',1e3);
xlabel(''); ylabel('Volume (mL)'); ylim([0 50]); title({'Tumor','Volume'}); yticks([0 25 50]);

% Number of Cancer Cells
subplot(3,3,4);
simbio_plot(simData_noTreat,'C1','legend','Control');
hold on; simbio_plot(simData,'C1','legend','Anti-PD1 Treatment');
xlabel(''); ylabel('Cell Count'); title({'Cancer','Cells'}); ylim([0 1e10]);

% Number of Naive T Cells
subplot(3,3,3); V_C = 5e6; % volume in uL
simbio_plot(simData_noTreat,'nT1','compartment','V_C','legend','Control','normalizeBy',V_C/1.11e6);
hold on; simbio_plot(simData,'nT1','compartment','V_C','legend','Anti-PD1 Treatment','normalizeBy',V_C/1.11e6);
xlabel(''); ylabel({'Cell Density','(cell/$\mu$L)'}); title({'Na\"ive T cells in','Blood'});

% Number of T Cells in the Blood
subplot(3,3,6);
simbio_plot(simData_noTreat,'T0','compartment','V_T','legend','Control','normalizeBy','V_T');
hold on; simbio_plot(simData,'T0','compartment','V_T','legend','Anti-PD1 Treatment','normalizeBy','V_T');
xlabel(''); ylabel({'Cell Density','(cell/$\mu$L)'}); title({'Activated Tregs in','Tumor'});

% Number of T Cells in the Tumour
subplot(3,3,9);
simbio_plot(simData_noTreat,'T1','compartment','V_T','legend','Control','normalizeBy','V_T');
hold on; simbio_plot(simData,'T1','compartment','V_T','legend','Anti-PD1 Treatment','normalizeBy','V_T');
xlabel(''); xlabel(''); ylabel({'Cell Density','(cell/$\mu$L)'}); title({'Activated T cells in','Tumor'});

% Immune Checkpoint Inhibitors
subplot(3,3,8); hold on;
simbio_plot(simData,'nivolumab','compartment','V_T','legend','nivolumab','normalizeBy',1e-6/1.436e2*0.55);
simbio_plot(simData,'nivolumab','compartment','V_T','legend','nivolumab','normalizeBy',1e-6/1.436e2*0.55);
ylabel({'Concentration','($\mu$g/mL)'}); title({'Anti-PD1 in','Blood'});

% APCs
subplot(3,3,7);
simbio_plot(simData_noTreat,'APC','compartment','V_T','linespec','k-');
hold on; simbio_plot(simData_noTreat,'APC','compartment','V_LN','linespec','k--');
hold on; simbio_plot(simData_noTreat,'mAPC','compartment','V_T','linespec','k-.');
hold on; h_ = simbio_plot(simData_noTreat,'mAPC','compartment','V_LN','linespec','k:');
set(h_,'LineWidth',4);
set(gca,'ColorOrderIndex',1);
hold on; simbio_plot(simData_noTreat,'APC','compartment','V_T','linespec','-');
hold on; simbio_plot(simData,'APC','compartment','V_T','linespec','-');
set(gca,'ColorOrderIndex',1);
hold on; simbio_plot(simData_noTreat,'APC','compartment','V_LN','linespec','--');
hold on; simbio_plot(simData,'APC','compartment','V_LN','linespec','--');
set(gca,'ColorOrderIndex',1);
hold on; simbio_plot(simData_noTreat,'mAPC','compartment','V_T','linespec','-.');
hold on; simbio_plot(simData,'mAPC','compartment','V_T','linespec','-.');
set(gca,'ColorOrderIndex',1);
hold on; h(1) = simbio_plot(simData_noTreat,'mAPC','compartment','V_LN','linespec',':');
hold on; h(2) = simbio_plot(simData,'mAPC','compartment','V_LN','linespec',':');
set(h,'LineWidth',4);
title({'Antigen Presenting','Cells'}); xlabel(''); ylabel('Cell Count');
% legend({'APC in Tumor','APC in LN','mAPC in Tumor','mAPC in LN'},'fontSize',10);
set(gca,'Yscale','log'); ylim([1e0 1e15]); yticks([1e0 1e3,1e6,1e9,1e12 1e15]);

% Antigen
subplot(3,3,2);
simbio_plot(simData_noTreat,'P1','compartment','V_T','legend','Control','normalizeBy',1e-6);
hold on; simbio_plot(simData,'P1','compartment','V_T','legend','Anti-PD1','normalizeBy',1e-6);
xlabel(''); ylabel({'Concentration','($\mu$M)'}); title('Free Antigen');
plot([],[]);
legend('Location','northWest','fontSize',10); ylim([0 3.5]);

% Antigen Presentation
subplot(3,3,5); A_s = 900; % APC surface area in um^2
simbio_plot(simData_noTreat,'M1p1','compartment','A_s','legend','Control','normalizeBy',1/A_s);
hold on; simbio_plot(simData,'M1p1','compartment','A_s','legend','Anti-PD1 Treatment','normalizeBy',1/A_s);
title('Antigen Presentation'); xlabel(''); ylabel({'Number of','Molecules'});