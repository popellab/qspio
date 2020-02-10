function figure_defaults()

% Text Defaults
set(groot,'defaultTextInterpreter','Latex');
set(groot,'defaultTextFontSize',18);
set(groot,'defaultTextColor','k');
set(groot,'defaultTextHorizontalAlignment','left');
set(groot,'defaultTextVerticalAlignment','top');

% Axes Defaults
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultAxesFontSize',12); %20
set(groot,'defaultAxesLabelFontSizeMultiplier',1.5);
set(groot,'defaultAxesBox','on');
%set(groot,'defaultAxesColorOrder',[nx3]);
%set(groot,'defaultAxesLineStyleOrder',[]);

% Legend Defaults
set(groot,'defaultLegendInterpreter','Latex');
set(groot,'defaultLegendFontSizeMode','manual','defaultLegendFontSize',15);%15
set(groot,'defaultLegendBox','off');

% Colorbar Defaults
set(groot,'defaultColorbarTickLabelInterpreter','Latex');
set(groot,'defaultColorbarFontSizeMode','manual','defaultColorbarFontSize',15);

% Line Plot Defaults
set(groot,'defaultLineLineWidth',1.8);

% Surface Plot Defaults
set(groot,'defaultSurfaceEdgeColor','none');

% Contour Plot Defaults
set(groot,'defaultContourLineWidth',1.8);
set(groot,'defaultContourLineColor','k');
set(groot,'defaultContourLabelSpacing',3000);