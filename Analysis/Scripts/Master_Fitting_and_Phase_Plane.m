% Connor Tynan- Jan 2023 - Final Year MMath Project
% Script for visualising parameter estimates

% Figure(1): gLV Estimate Trajectories
% Figure(2): Comparing Non-Linear Co-Culture Trajectories
% Figure(3): Phase Plane and Nullclines of 2D gLV System

% inters = [(M12-M22)/((M12/M11)-(M22/M21)),(M21-M11)/((M21/M22)-(M11/M12))];
% % 132391.610038639
% % -3459.22902814445

% Miscellaneous
clc; clear; clf; close all; format short

% Importing mean data
[ms1,ms2,ms3,ms4] = importmeandata;

% Time parameters
t0 = 0; t1 = 13; dt = 0.001;

% Parameter estimates (Calculated via Bayesian Inferencing (R))

% % Case 1a - Mean Averaged
% % Parameter Estimates
mu1 = 1.01; M11 = 16897.31; M12 = 506.1;
mu2 = 0.87; M22 = 13789.37; M21 = 105840.3;
% Mean Standard Deviation (Error)
err = [3817.02, 1676.4, 38.3, 405474.3]; 

% Case 2
% Parameter Estimates
% mu1 = 0.97; M11 = 1716875.51; M12 = 3348.39;
% mu2 = 1.12; M22 = 2615720.70; M21 = 492545.63;
% % Mean Standard Deviation (Error)
% err = [0, 0, 40, 30000]; 

% Mono-Culture calculations
[sam,pam] = monoplot(t0,dt,t1,mu1,mu2,M11,M22);

% Co-Culture calculations
[sac,pac] = coplot(t0,dt,t1,mu1,mu2,M11,M22,M12,M21);

% Time Series Estimation - Estimate Time Series + Mean SD Error

% S. Aureus Mono-Culture Plot
sp1 = subplot(2,2,1); hold on; grid minor; grid on; box on;
% Transparent Error Bars
fill([(t0:dt:t1)';flipud((t0:dt:t1)')],[sam-err(1);flipud(sam+err(1))],[1 0.6 0.4], ...
    'LineStyle','none','FaceAlpha',0.3);
% Plotting Trajectory
plot(t0:dt:t1,sam,'LineWidth',1.7,'Color',[1 0.6 0.1],'DisplayName',' S. Aureus Mono Est.');
% Plotting Observed Data
plot(t0:1:t1,ms1/ms1(1),'xk'); 
% Plot Miscellaneous
ylim([0 3*10^4]);
legsp1 = legend('',"Mono-Culture",'Location','northwest');
legsp1.ItemTokenSize = [20,18];
title('S. Aureus');
sp1.TitleHorizontalAlignment = 'right'; 

% P. Aeruginosa Mono-Culture Plot
sp2 = subplot(2,2,2); hold on; grid minor; grid on; box on;
% Transparent Error Bars
fill([(t0:dt:t1)';flipud((t0:dt:t1)')],[pam-err(2);flipud(pam+err(2))],[.8 .7 .8], ...
    'LineStyle','none','FaceAlpha',0.3);
% Plotting Trajectory
plot(t0:dt:t1,pam,'LineWidth',1.7,'Color',[0.8 0.5 0.8],'DisplayName', ...
    ' P. Aeruginosa Mono Est.');
% Plotting Observed Data
plot(t0:1:t1,ms2/ms2(1),'xk'); 
% Plot Miscellaneous
ylim([0 15000]);
legsp2 = legend('',"Mono-Culture",'Location','northwest');
legsp2.ItemTokenSize = [20,18];
title('P. Aeruginosa');
sp2.TitleHorizontalAlignment = 'right'; 
ax = gca; % axes handle
ax.YAxis.Exponent = 4;

% S. Aureus Co-Culture Plot
sp3 = subplot(2,2,3); hold on; grid minor; grid on; box on;
% Transparent Error Bars
fill([(t0:dt:t1)';flipud((t0:dt:t1)')],[sac-err(3);flipud(sac+err(3))],[1 0.6 0.4], ...
    'LineStyle','none','FaceAlpha',0.3);
% Plotting Trajectory
plot(t0:dt:t1,sac,'--','LineWidth',1.7,'Color',[1 0.6 0.1],'DisplayName', ...
    ' S. Aureus Co Est.');
% Plotting Observed Data
plot(t0:1:t1,ms3/ms3(1),'xk'); 
% Plot Miscellaneous
ylim([0 600]);
legsp3 = legend('',"Co-Culture",'Location','northwest');
legsp3.ItemTokenSize = [20,18];
title('S. Aureus');
sp3.TitleHorizontalAlignment = 'right'; 

% P. Aeruginosa Co-Culture Plot
sp4 = subplot(2,2,4); hold on; grid minor; grid on; box on;
% Transparent Error Bars
fill([(t0:dt:t1)';flipud((t0:dt:t1)')],[pac-err(4);flipud(pac+err(4))],[.8 .7 .8], ...
    'LineStyle','none','FaceAlpha',0.3);
% Plotting Trajectory
plot(t0:dt:t1,pac,'--','LineWidth',1.7,'Color',[0.8 0.5 0.8],'DisplayName', ...
    ' P. Aeruginosa Co Est.');
% Plotting Observed Data
plot(t0:1:t1,ms4/ms3(1),'xk'); 
% Plot Miscellaneous
ylim([0 1.5*10^6]);
legsp4 = legend('',"Co-Culture",'Location','northwest');
legsp4.ItemTokenSize = [20,18];
title('P. Aeruginosa');
sp4.TitleHorizontalAlignment = 'right'; 
ax = gca; % axes handle
ax.YAxis.Exponent = 4;

% Figure Configuration
subfig = axes(gcf,'visible','off');
subfig.Title.Visible='on';
subfig.XLabel.Visible='on';
subfig.YLabel.Visible='on';
set(gca, 'fontname', 'Helvetica');
ylabel(subfig,'Bacteria count (mL^{-1}) - Mean Averaged')
xlabel(subfig,'Time (hr)')

% Export Settings
set(gcf, 'units', 'centimeters', 'position', [1 1 12.1 9]);
% export_fig('1a_MA_traj_fig', gcf, '-transparent', '-png', '-r500')

% Phase Plane and Nullclines of 2D gLV System

figure(2); hold on; box on; % Figure Miscellaneous

% Generating Phase Plane
% Function to Calculate Phase Arrows for
fphase = @(t,Y) [mu1*Y(1)*(1-Y(1)/M11-Y(2)/M12);
mu2*Y(2)*(1-Y(2)/M22-Y(1)/M21)];
% Vectors Outlining Phase Plane Space
y1 = linspace(0, 150000,20);
y2 = linspace(-.7*10^4, 1.7*10^4,20);
% Generating a Mesh Grid
[m,n] = meshgrid(y1,y2);
size(m)
size(n)

% Solving Function - Generating Magnitude and Direction of Flow
% Initial Conditions
u = zeros(size(m));
v = zeros(size(n));
t=0; 
% Main Computation
for i = 1:numel(m)
    Yprime = fphase(t,[m(i); n(i)]);
    u(i) = Yprime(1);
    v(i) = Yprime(2);
end

% Generating Phase Plane
scale = 1; q = quiver(m,n,u,v,'r'); q.AutoScaleFactor = 1.5;

% Dummy Vector for Nullclines
x = linspace(-1*10^10, 1.5*10^10, 1001);

% System Nullclines
plot(0*x,x,'LineWidth',1.7,'color',[0.8 0 0 1])
plot(x,0*x,'LineWidth',1.7,'color',[0 0 0.8 1])
plot(x, M12*(1-(x)/M11),'LineWidth', 2,'color',[0.8 0 0 1])
plot(x, M22*(1-(x)/M21),'LineWidth', 2,'color',[0 0 0.8 1])

% Plotting Estimate Trajectory (figure(1) and figure(2))
plot(sac,pac,'--k','LineWidth', 2)

% Finding Equilibrium Points (Intersection of Nullclines)
P1 = InterX([x; M12*(1-(x)/M11)],[x; M22*(1-(x)/M21)]); % Intersection Saddle
P2 = InterX([x; M12*(1-(x)/M11)],[x; x*0]); % Axis Saddle Node
P3 = InterX([x; M22*(1-(x)/M21)],[x*0; x]); % Axis Steady State

text(-9500,-1100,'Unstab.')
text(P1(1),P1(2)*0.65,'Stab.')
text(P2(1),P2(2)-1100,'Saddle')
text(P3(1)+2500,P3(2)+1000,'Stab.')

% Plotting Equilibrium Points
plot(P1(1),P1(2)','xg','MarkerSize',10, 'linewidth',2)
plot(P2(1),P2(2)','xg','MarkerSize',10, 'linewidth',2)
plot(P3(1),P3(2)','xg','MarkerSize',10, 'linewidth',2)
plot(0,0,'xg','MarkerSize',10, 'linewidth',2)

% Figure Configuration
% Legend
phaseleg = legend(' Phase Arrows','','',' Nullcline (S. Aureus)', ...
    ' Nullcline (P. Aeruginosa)', ' Co-Culture Trajectory', ...
    ' Equilibrium Points','location','northeast');
phaseleg.ItemTokenSize = [20,18];
% Title
title('System Nullclines'); ax = gca; ax.TitleHorizontalAlignment = 'right'; 
% Axis
xlim([-10000 1.5*10^5]); xlabel('S. Aureus Bacteria Count (y_{1})')
ylim([-.75*10^4 1.75*10^4]); ylabel('P. Aeruginosa Bacteria Count (y_{2})')

ax = gca; % axes handle
ax.YAxis.Exponent = 4;

% Export Settings
set(gcf, 'units', 'centimeters', 'position', [1 1 12.1 9]);
% export_fig('figure_for_today', gcf, '-transparent', '-png', '-r500')
