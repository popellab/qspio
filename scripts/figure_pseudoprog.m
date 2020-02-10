% Pseudo-Progression Figure - script to plot pseudo-progression figure
% (Figure 2 in Sove et al 2020)

% Create Model
example1;

% Modify Parameters
variantObj = addvariant(model,'v1');
addcontent(variantObj,{'parameter','k_cell_clear','Value',0.001});
addcontent(variantObj,{'parameter','k_C_T1','Value',8.7});

% Run Simulation
[model,success] = initial_conditions(model,'Variant',variantObj);
simData = sbiosimulate(model,[],variantObj,dose_schedule);

figure; 
% Tumour Volume
subplot(1,3,1); hold on; simbio_plot(simData,'V_T','normalizeBy',1e3);
set(gca,'ColorOrderIndex',4); simbio_plot(simData,'V_T','normalizeBy',1e3);  axis square;
title('Tumor Volume'); ylabel('Volume (mL)'); xlim([0,300]);
% Number of Cancer Cells
subplot(1,3,2); hold on; simbio_plot(simData,'C1'); simbio_plot(simData,'C_x'); axis square;
title('Cancer Cells'); ylabel('Cell Count'); xlim([0,300]);
% Number of T Cells
subplot(1,3,3); hold on; simbio_plot(simData,'T1','compartment','V_T','normalizeBy','V_T');
simbio_plot(simData,'T_exh','normalizeBy','V_T'); axis square; set(gca,'Yscale','log');
title({'T Cells in','the Tumor'}); ylabel('Cell Count (cell/$\mu$L)'); xlim([0,300]); yticks([1e1,1e2,1e3,1e4,1e5]);
legend('viable','dead','location','east');