% Connor Tynan - Visualising Parameter Estimates

% Miscellaneous
clc; clf; close all; clear;

% Importing mono-culture data
[mono_sa, mono_pa] = importrawdata;

% Normalising the data
for i = 1:6
    mono_sa(:,i) = mono_sa(:,i)/mono_sa(1,i);
    mono_pa(:,i) = mono_pa(:,i)/mono_pa(1,i);
end

% Dummy vector for plotting
z = linspace(-1e20,1e20,1001);

% Time parameters
t0 = 0; t1 = 13;

% Arrays for parameter estimates (median values)

mu1 = [0.91, 0.86, 0.99, 1.17, 1.01, 1.08];
mu2 = [0.94, 0.96, 0.85, 0.85, 0.76, 0.74];
L11 = [28800.88, 33287.12, 16604.13, 15629.94, 12193.25, 8416.49];
L22 = [17663.69, 12626.88, 9098.26, 3310.00, 152428.87, 163773.34];
L12 = zeros(1,6);
L21 = zeros(1,6);

figure(1)
% Initialising 2x3 subplot
set(gca,'Visible','off')
ha = tight_subplot(2,3,[.125 .065],[.08 .08],[.07 .075]);

for i = 1:6

    % S.aureus mono-culture
    axes(ha(i));
    
    % ODEs to solve
    F = @(t,y) [mu1(i)*y(1)-mu1(i)*y(1)^2/L11(i);
    mu2(i)*y(2)-mu2(i)*y(2)^2/L22(i)];

    [t,xts] = ode45(F, [t0 t1], [1 1]);
    
    % Quiver Initial Conditions
    X=linspace(L11(i)*0.03,L11(i)*0.97,25); 
    Y=linspace(L22(i)*0.03,L22(i)*0.97,25); 
    
    % Initialising Grid for Quiver
    [m,n] = meshgrid(X,Y);
    % Initialising Quiver Solution Vectors
    u = zeros(size(m));
    v = zeros(size(n));
    t = 0;                  % Dummy time value
    
    % Quiver Computation
    
    for j = 1:numel(m)
    
        % Solving System of Equations
        Yprime = F(t,[m(j); n(j)]);
        
        % Storing the solutions
        u(j) = Yprime(1);
        v(j) = Yprime(2);
    end
    
    % Generating Quiver Plot --- Grid of points & Grid of solutions
    scale = 1; q = quiver(m,n,u,v,'color',[255/255,140/255,0,0.6]);

    hold on

    % Plotting Phase Diagram
    plot(xts(:,1),xts(:,2),'r--','linewidth',2);    % Trajectory
    plot(mono_sa(:,i),mono_pa(:,i),'xk','MarkerSize',10,'LineWidth',1.5) % Observed Data
    plot(0,0,'color','#1E88E5','linewidth',2)       % Dummy Nullcline for Legend
    plot(0,0,'color','#FFC107','linewidth',2)       % Dummy Nullcline for Legend
    
    % Equilibrium points
    plot(0,L22(i),'o','color',[0.9 0 0],'MarkerFaceColor',[0.9 0 0],'MarkerSize',10)
    plot(L11(i),0,'o','color',[0.9 0 0],'MarkerFaceColor',[0.9 0 0],'MarkerSize',10)
    plot(L11(i),L22(i),'o','color',[0.7 0 0.7],'MarkerFaceColor',[0.7 0 0.7],'MarkerSize',10)
    plot(0,0,'o','color','#FFD1DC','MarkerFaceColor','#FFD1DC','MarkerSize',10)

    % y1 and y2 Nullclines
    xh = xline(0,'color','#1E88E5','linewidth',2,'DisplayName','Nullcline');
    xh2 = xline(L11(i),'color','#1E88E5','linewidth',2);
    xh3 = yline(0,'color','#FFC107','linewidth',2); 
    xh4 = yline(L22(i),'color','#FFC107','linewidth',2);
    % Formatting of Nullclines
    Sxh = struct(xh);
    Sxh.Edge.Layer = 'back';
    Sxh = struct(xh2);
    Sxh.Edge.Layer = 'back';
    Sxh = struct(xh3);
    Sxh.Edge.Layer = 'back';
    Sxh = struct(xh4);
    Sxh.Edge.Layer = 'back';
    
    Xlimit = xlim; Ylimit = ylim;
    xlim([-L11(i)/10 L11(i)*1.2])
    ylim([-Ylimit(2)/10 L22(i)*1.2])

    ax=gca; ax.XAxis.Exponent = 4; ax.YAxis.Exponent = 4;
    

    % Plot aesthetics
    LimitsX = xlim; LimitsY = ylim;
    title(sprintf('Replica %d',i),'HorizontAlalignment', 'right', 'position', [LimitsX(2), LimitsY(2)])

%     leg = legend(' Phase Arrows',' Trajectory',' Observed Data',' y_{1} Nullcline',' y_{2} Nullcline', ...
%         ' Saddle Point','',' Stable Node',' Unstable Node','location','east');

%     leg.ItemTokenSize = [20,18];

    box on; xlabel('y_{1}'); ylabel('y_{2}')
    set(gca,'FontSize',12)
    figure(2)
    subplot(1,2,1)
    plot(t,xts(:,1))
    hold on
    plot(t0:1:t1,mono_sa(:,i),'xk')
    subplot(1,2,2)
    plot(t,xts(:,2))
    hold on
    plot(t0:1:t1,mono_pa(:,i),'xk')
end

figure(1)
set(gcf, 'units', 'centimeters', 'position', [1 1 12.1*2.5 9*2]);

