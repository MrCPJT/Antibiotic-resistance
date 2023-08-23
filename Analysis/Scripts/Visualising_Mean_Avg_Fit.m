% Connor Tynan - Visualising Mean-average Parameter Estimates
% We solve the gLV equations using MATLAB's in-built ODE45 solver.
% We can then visualise the time series and phase plane 

% Miscellaneous
clc; clf; close all; clear;

% Importing mono-culture data
[mono_sa, mono_pa] = importrawdata;
[mono_sa_mean, mono_pa_mean] = importmeandata;

% Time parameters
t0 = 0; t1 = 25;

% Arrays for parameter estimates (median values)
mu1 = [0.91, 0.86, 0.99, 1.17, 1.01, 1.08];
mu2 = [0.94, 0.96, 0.85, 0.85, 0.76, 0.74];
L11 = [28800.88, 33287.12, 16604.13, 15629.94, 12193.25, 8416.49];
L22 = [17663.69, 12626.88, 9098.26, 3310.00, 152428.87, 163773.34];
L12 = zeros(1,6);
L21 = zeros(1,6);

% Parameters fit to mean data
mu1m = 1.01;
mu2m = 0.87;
L11m = 16724.89;
L22m = 13790.53;
L12m = 0;
L21m = 0;

% Parameter Estimates using Mean-averaged Mono-Culture Data
figure(1)

% Customising Figure Layout
ha = tight_subplot(1,1,[.125 .065],[.15 .06],[.25 .25]);
axes(ha);

hold on

% Highlighting biologically infeasible region
p1 = patch([1e5 -1e5 -1e5 1e5], [-1e5 -1e5, 0 0], [0.8 0.8 0.8], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0 0 0]);
p2 = patch([0 -1e5 -1e5 0], [0 0, 1e5 1e5], [0.8 0.8 0.8], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0 0 0]);

% ODEs to solve
Fm = @(t,y) [mu1m.*y(1)-mu1m.*y(1).*y(1)./L11m;
             mu2m.*y(2)-mu2m.*y(2).*y(2)./L22m];

% Solving the ODEs in time
[t,xts] = ode45(Fm, [t0 t1], [1 1]);

% Quiver Initial Conditions
X=linspace(L11m*0.03,L11m*0.97,25); 
Y=linspace(L22m*0.03,L22m*0.97,25); 

% Plotting Phase Arrows (See bottom of script)
pp_plot(X,Y,mu1m,mu2m,L11m,L22m)

% Quiver Initial Conditions - Lower Right Quadrant
X=linspace(L11m*1.03,L11m*1.97,25); 
Y=linspace(L22m*0.03,L22m*0.97,25); 

% Plotting Phase Arrows (See bottom of script)
pp_plot(X,Y,mu1m,mu2m,L11m,L22m)

% Quiver Initial Conditions - Upper Left Quadrant
X=linspace(L11m*0.03,L11m*0.97,25); 
Y=linspace(L22m*1.03,L22m*1.97,25); 

% Plotting Phase Arrows (See bottom of script)
pp_plot(X,Y,mu1m,mu2m,L11m,L22m)

% Quiver Initial Conditions - Upper Right Quadrant
X=linspace(L11m*1.03,L11m*1.97,25); 
Y=linspace(L22m*1.03,L22m*1.97,25); 

% Plotting Phase Arrows (See bottom of script)
pp_plot(X,Y,mu1m,mu2m,L11m,L22m)

% Plotting Phase Diagram
plot(xts(:,1),xts(:,2),'r--','linewidth',2);    % Trajectory
plot(mono_sa_mean,mono_pa_mean,'xk','MarkerSize',10,'LineWidth',1.5) % Observed Data
plot(0,0,'color','#1E88E5','linewidth',2)       % Dummy Nullcline for Legend
plot(0,0,'color','#FFC107','linewidth',2)       % Dummy Nullcline for Legend

% Equilibrium points
plot(0,L22m,'o','color',[0.9 0 0],'MarkerFaceColor',[0.9 0 0],'MarkerSize',10)
plot(L11m,0,'o','color',[0.9 0 0],'MarkerFaceColor',[0.9 0 0],'MarkerSize',10)
plot(L11m,L22m,'o','color',[0.7 0 0.7],'MarkerFaceColor',[0.7 0 0.7],'MarkerSize',10)
plot(0,0,'o','color','#FFD1DC','MarkerFaceColor','#FFD1DC','MarkerSize',10)

% y1 and y2 Nullclines
xh = xline(0,'color','#1E88E5','linewidth',2,'DisplayName','Nullcline');
xh2 = xline(L11m,'color','#1E88E5','linewidth',2);
xh3 = yline(0,'color','#FFC107','linewidth',2); 
xh4 = yline(L22m,'color','#FFC107','linewidth',2);

% Formatting of Nullclines
Sxh = struct(xh);
Sxh.Edge.Layer = 'back';
Sxh = struct(xh2);
Sxh.Edge.Layer = 'back';
Sxh = struct(xh3);
Sxh.Edge.Layer = 'back';
Sxh = struct(xh4);
Sxh.Edge.Layer = 'back';

% Axis Aesthetics
Xlimit = xlim; Ylimit = ylim;
xlim([-L11m*0.1 L11m*1.1])
ylim([-L22m*0.1 L22m*1.1])

ax=gca; ax.XAxis.Exponent = 4; ax.YAxis.Exponent = 4;

% Legend
leg = legend('','','','','',' Phase Arrows',' Trajectory',' Observed Data',' y_{1} Nullcline',' y_{2} Nullcline', ...
    ' Saddle Point','',' Stable Node',' Unstable Node','location','west');
leg.ItemTokenSize = [20,18];

% Plot aesthetics
LimitsX = xlim; LimitsY = ylim;
box on; xlabel('y_{1}'); ylabel('y_{2}')
set(gca,'FontSize',12)
set(gcf, 'units', 'centimeters', 'position', [1 1 12.1*2 9*1.35]);
% export_fig('Mean_Avg_Param_Est_Phase_Plane', gcf, '-transparent', '-png', '-r125')

% Solving gLV equations in time with parameter estimated coefficients
% 1x2 Subplot showcasing S.aureus and P.aeruginosa trajectories against
% observed data

% Initialising new figure
fig2 = figure(2);

% Subplot 1: S.aureus
subplot(1,2,1); hold on; box on
% Plotting trajectory
plot(t,xts(:,1),'LineWidth',2,'color','#1E88E5')
% Plotting data for comparison/reference
plot(0:1:13,mono_sa_mean/mono_sa_mean(1),'kx','MarkerSize',8,'LineWidth',2)
% Axis
xlim([0 15])
xlabel('Time [hr]'); ylabel('y_{1}')
% Legend
leg = legend(' Trajectory',' Data','location','northwest');
leg.ItemTokenSize = [20,18];
% Plot Aesthetics
set(gca,'FontSize',12)
ax=gca; ax.YAxis.Exponent = 4;

% Subplot 2: P.aeruginosa
subplot(1,2,2); hold on; box on
% Plotting trajectory
plot(t,xts(:,2),'LineWidth',2,'color','#FFC107')
% Plotting data for comparison/reference
plot(0:1:13,mono_pa_mean/mono_pa_mean(1),'kx','MarkerSize',8,'LineWidth',2)
% Axis
xlim([0 15])
xlabel('Time [hr]'); ylabel('y_{2}')
% Legend
leg = legend(' Trajectory',' Data','location','northwest');
leg.ItemTokenSize = [20,18];
% Plot Aesthetics
set(gca,'FontSize',12)
ax=gca; ax.YAxis.Exponent = 4;
% Legend
leg = legend(' Trajectory',' Data','location','northwest');
leg.ItemTokenSize = [20,18];

% Figure Aesthetics
set(gcf, 'units', 'centimeters', 'position', [1 1 12.1*2 9]);
% export_fig('Mean_Avg_Param_Est_Time_Series', gcf, '-transparent', '-png', '-r125')

%%

% Function to plot phase plane
function q = pp_plot(X,Y,mu1m,mu2m,L11m,L22m)

    % ODEs to solve
    Fm = @(t,y) [mu1m.*y(1)-mu1m.*y(1).*y(1)./L11m;
                 mu2m.*y(2)-mu2m.*y(2).*y(2)./L22m];
    
    % Solving the ODEs in time
    [t,xts] = ode45(Fm, [0 13], [1 1]);

    % Initialising Grid for Quiver - Lower Left Quadrant
    [m,n] = meshgrid(X,Y);
    
    % Initialising Quiver Solution Vectors
    u = zeros(size(m));
    v = zeros(size(n));
    
    % Quiver Computation - Lower Left Quadrant
    
    for j = 1:numel(m)
    
        % Solving System of Equations
        Yprime = Fm(t,[m(j); n(j)]);
        
        % Storing the solutions
        u(j) = Yprime(1);
        v(j) = Yprime(2);
    end
    
    % Generating Quiver Plot --- Grid of points & Grid of solutions
    q = quiver(m,n,u,v,'color',[255/255,140/255,0,0.6]);

end