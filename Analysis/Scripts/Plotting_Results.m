% Connor Tynan - Visualising Mean-average Parameter Estimates
% Plotting the results of cases 1a-d

% Parameter Estimates using Mean-averaged Normalised Co-Culture Data

% Miscellaneous
clc; clf; close all; clear; hold on

% Importing mono and co-culture data
[mono_sa, mono_pa, co_sa, co_pa] = importrawdata;
[mono_sa_mean, mono_pa_mean, co_sa_mean, co_pa_mean] = importmeandata;

% Time parameters
t0 = 0;     % Initial time
t1 = 25;    % Final time

% Arrays for parameter estimates (median values)
mu1 = [1.01, 0.99, 1.01, 0.97];
mu2 = [0.87, 1.28, 0.87, 1.12];
L11 = [16724.89, 16724.89, 97130.57, 94155.04];
L22 = [13790.53, 13790.53, 121839.35, 2394258.29];
L12 = [460.64, 5002.28, 460.18, 2813.77];
L21 = [105897.27, -321331.69, 106665.19, 103731.30];

% Main computation
% - Plotting phase plane and time series of each population for each case

for i = 1:4
    
    % Initialising figure
    figure(i); hold on
    
    % Highlighting biologically infeasible region
    p1 = patch([1e8 -1e8 -1e8 1e8], [-1e8 -1e8, 0 0], [0.8 0.8 0.8], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0 0 0]);
    p2 = patch([0 -1e8 -1e8 0], [0 0, 1e8 1e8], [0.8 0.8 0.8], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0 0 0]);
    
    % ODEs to solve
    Fm = @(t,y) [mu1(i).*y(1)-mu1(i).*y(1).*y(1)./L11(i)-mu2(i).*y(1).*y(2)./L12(i);
                 mu2(i).*y(2)-mu2(i).*y(2).*y(2)./L22(i)-mu1(i).*y(1).*y(2)./L21(i)];
    
    % Solving the ODEs in time
    [t,xts] = ode45(Fm, [t0 t1], [1 1]);
    
    % Calculating non-zero equilibrium position
    xinter = (L12(i)-L22(i))*L11(i)*L21(i)/(L12(i)*L21(i)-L11(i)*L22(i));
    yinter = (L21(i)-L11(i))*L22(i)*L12(i)/(L12(i)*L21(i)-L11(i)*L22(i));

    % Quiver Initial Conditions
    X=linspace(max([L11(i) abs(xinter)])*0.02,max([L11(i) abs(xinter)])*1.05,25); 
    Y=linspace(max([L22(i) yinter])*0.02,max([L22(i) yinter])*1.05,25); 
    
    % Plotting Phase Arrows (See bottom of script)
    pp_plot(X,Y,mu1(i),mu2(i),L11(i),L22(i),L12(i),L21(i),i)
    
    % Quiver Initial Conditions - Lower Right Quadrant
    X=linspace(max([L11(i) xinter])*0.02,max([L11(i) xinter])*1.05,25);
    Y=-linspace(max([L22(i) yinter])*0.02,max([L22(i) yinter])*1.05,25); 
    
    % Plotting Phase Arrows (See bottom of script)
    pp_plot(X,Y,mu1(i),mu2(i),L11(i),L22(i),L12(i),L21(i),i)
    
    % Quiver Initial Conditions - Upper Left Quadrant
    X=-linspace(max([L11(i) abs(xinter)])*0.02,max([L11(i) abs(xinter)])*1.05,25);
    Y=linspace(max([L22(i) yinter])*0.02,max([L22(i) yinter])*1.05,25); 
    
    % Plotting Phase Arrows (See bottom of script)
    pp_plot(X,Y,mu1(i),mu2(i),L11(i),L22(i),L12(i),L21(i),i)
    
    % Quiver Initial Conditions - Upper Right Quadrant
    X=-linspace(max([L11(i) xinter])*0.02,max([L11(i) xinter])*1.05,25);
    Y=-linspace(max([L22(i) yinter])*0.02,max([L22(i) yinter])*1.05,25); 
    
    % Plotting Phase Arrows (See bottom of script)
    pp_plot(X,Y,mu1(i),mu2(i),L11(i),L22(i),L12(i),L21(i),i)
    
    % Plotting Phase Diagram
    plot(xts(:,1),xts(:,2),'r--','linewidth',2);    % Trajectory
    plot(co_sa_mean/co_sa_mean(1),co_pa_mean/co_pa_mean(1),'xk', ...
        'MarkerSize',10,'LineWidth',1.5)            % Observed Data
    plot(0,0,'color','#1E88E5','linewidth',2)       % Dummy Nullcline for Legend
    plot(0,0,'color','#FFC107','linewidth',2)       % Dummy Nullcline for Legend
    
    % Plotting nullclines
    x = linspace(-1e5,1e6,10000);   % Dummy array

    % Nullclines
    plot(L11(i)*(1-x/L12(i)),x,'k','LineWidth',2,'color','#1E88E5');
    plot(x,L22(i)*(1-x/L21(i)),'k','LineWidth',2,'color','#FFC107');

    % Equilibrium points
    plot(0,L22(i),'o','color',[0.7 0 0.7],'MarkerFaceColor',[0.7 0 0.7],'MarkerSize',10)
    plot(L11(i),0,'o','color',[0.9 0 0],'MarkerFaceColor',[0.9 0 0],'MarkerSize',10)
    plot(0,0,'o','color','#FFD1DC','MarkerFaceColor','#FFD1DC','MarkerSize',10)
    
    % Non-zero equilibria
    plot((L12(i)-L22(i))*L11(i)*L21(i)/(L12(i)*L21(i)-L11(i)*L22(i)), ...
        (L21(i)-L11(i))*L22(i)*L12(i)/(L12(i)*L21(i)-L11(i)*L22(i)),'o', ...
        'color',[0.9 0 0],'MarkerFaceColor',[0.9 0 0],'MarkerSize',10)

    % Zero-axis nullclines
    xh  = xline(0,'color','#1E88E5','linewidth',2,'DisplayName','Nullcline');
    xh2 = yline(0,'color','#FFC107','linewidth',2); 
    
    % Formatting of Nullclines
    Sxh = struct(xh);
    Sxh.Edge.Layer = 'back';
    Sxh = struct(xh2);
    Sxh.Edge.Layer = 'back';
    
    % Axis Aesthetics
    Xlimit = xlim; Ylimit = ylim;
    xlim([min([-0.5*10^4 xinter])*1.05 max([L11(i) xinter])*1.05])
    ylim([min([L22(i)*-0.05 yinter])*1.05 max([L22(i) yinter])*1.05])
    
    ax=gca; ax.XAxis.Exponent = 4; ax.YAxis.Exponent = 4;
    
    % Legend
    leg = legend('','','','','',' Phase Arrows',' Trajectory',' Observed Data',' y_{1} Nullcline',' y_{2} Nullcline', ...
        '','',' Stable Node',' Saddle Point',' Unstable Node','location','best');
    leg.ItemTokenSize = [20,18];
    
    % Plot aesthetics
    % Title
    LimitsX = xlim; LimitsY = ylim;
    cases = {'Case 1a', 'Case 1b', 'Case 1c', 'Case 1d'};
    title(cases(i),'HorizontAlalignment', 'right', 'position', [LimitsX(2), LimitsY(2)])

    box on; xlabel('y_{1}'); ylabel('y_{2}')
    set(gca,'FontSize',14)
    set(gcf, 'units', 'centimeters', 'position', [1 1 12.1*1.25 9*1.25]);
    % export_fig('Mean_Avg_Param_Est_Phase_Plane', gcf, '-transparent', '-png', '-r125')

    fig5 = figure(5); hold on
    
    % Subplot 1: S.aureus
    subplot(1,2,1); hold on; box on
    % Plotting trajectory
    if i == 1
        pp1 = plot(t,xts(:,1),'-.','LineWidth',2,'color','#648FFF');
    elseif i == 2
        pp1 = plot(t,xts(:,1),'--','LineWidth',2,'color','#FE6100');
    elseif i == 3
        pp1 = plot(t,xts(:,1),':','LineWidth',2,'color','#FFB000');
    else
        pp1 = plot(t,xts(:,1),'-','LineWidth',2,'color','#E72BF5');
    end
    % Plotting data for comparison/reference
    plot(0:1:13,co_sa_mean/co_sa_mean(1),'kx','MarkerSize',8,'LineWidth',2)
    % Axis
    xlim([0 15]); ylim([0 0.05*10^4])
    xlabel('Time (hr)'); ylabel('y_{1}')
    % Legend
    leg = legend(' Case 1a','',' Case 1b','',' Case 1c','',' Case 1d',' Data', ...
        'location','NorthWest');
    leg.ItemTokenSize = [20,18];
    % Plot Aesthetics
    set(gca,'FontSize',14)
    ax=gca; ax.YAxis.Exponent = 4;
    
    % Subplot 2: P.aeruginosa
    subplot(1,2,2); hold on; box on
    % Plotting trajectory
    if i == 1
        pp2 = plot(t,xts(:,2),'-.','LineWidth',2,'color','#648FFF');
    elseif i == 2
        pp2 = plot(t,xts(:,2),'--','LineWidth',2,'color','#FE6100');
    elseif i == 3
        pp2 = plot(t,xts(:,2),':','LineWidth',2,'color','#FFB000');
    else
        pp2 = plot(t,xts(:,2),'-','LineWidth',2,'color','#E72BF5');
    end
    % Plotting data for comparison/reference
    plot(0:1:13,co_pa_mean/co_pa_mean(1),'kx','MarkerSize',8,'LineWidth',2)
    % Axis
    xlim([0 15]); ylim([0 60*10^4])
    xlabel('Time (hr)'); ylabel('y_{2}')
    % Legend
    leg = legend(' Trajectory',' Data','location','northwest');
    leg.ItemTokenSize = [20,18];
    % Plot Aesthetics
    set(gca,'FontSize',14)
    ax=gca; ax.YAxis.Exponent = 4;
    % Legend
    leg = legend(' Case 1a','',' Case 1b','',' Case 1c','',' Case 1d',' Data', ...
        'location','NorthWest');
    leg.ItemTokenSize = [20,18];

    % Figure Aesthetics
    set(gcf, 'units', 'centimeters', 'position', [1 1 12.1*2 9]);

end

% Solving gLV equations in time with parameter estimated coefficients
% 1x2 Subplot showcasing S.aureus and P.aeruginosa trajectories against
% observed data

% Initialising new figure
% export_fig('Mean_Avg_Param_Est_Time_Series', gcf, '-transparent', '-png', '-r125')

%%

% Function to plot phase plane
function q = pp_plot(X,Y,mu1,mu2,L11,L22,L12,L21,i)

    % ODEs to solve
    Fm = @(t,y) [mu1.*y(1)-mu1.*y(1).*y(1)./L11-mu2.*y(1).*y(2)./L12;
                 mu2.*y(2)-mu2.*y(2).*y(2)./L22-mu1.*y(1).*y(2)./L21];
    
    % Solving the ODEs in time
    [t,~] = ode45(Fm, [0 13], [1 1]);

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
    
    if i ~= 4
        % Generating Quiver Plot --- Grid of points & Grid of solutions
        q = quiver(m,n,u,v,'color',[255/255,140/255,0,0.6],'AutoScaleFactor',1.2);
    elseif i == 4
        q = quiver(m,n,u,v,'color',[255/255,140/255,0,0.6],'AutoScaleFactor',0.1, ...
            'ShowArrowHead','on','AutoScale','on');
    end
end