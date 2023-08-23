% Connor Tynan

% Housekeeping
clc; clf; close all

% Importing data
data = readtable('C:\Users\Tynan\Desktop\Maths\Year 4\Dissertation\Final Year Project\Parameter Estimation\Bayesian Inferencing\theta1data.csv');
data1 = readtable('C:\Users\Tynan\Desktop\Maths\Year 4\Dissertation\Final Year Project\Parameter Estimation\Bayesian Inferencing\theta2data.csv');

% Converting from table to matrix format
df = table2array(data);
df1 = table2array(data1);

% Initialising Solution Matrix
hist = zeros(20,15000);

% Colours for MCMC density plots
Colors = {'#ffb3ba', '#ffdfba', '#ffffba', '#baffc9', '#bae1ff', ...
    '#9a6324', '#fffac8', '#800000', '#3cb44b', '#ffe119', '#4363d8', ...
    '#f58231', '#911eb4', '#46f0f0', '#f032e6', '#aaffc3', '#808000', ...
    '#ffd8b1', '#000075', '#808080', '#ffffff', '#000000'};

% Initialising Figure
fig = figure(1); hold on;

% Main Computation
for j = 1:2

    % Counters for updating
    k = 1;
    x = 1;

    % Inner loop - Iterations
    for i = 1:300000
        
        % theta[1] Case
        if j == 1
            if mod(i,15000) == 0 && i > 1
                hist(k,x) = df(i,5);
                k = k + 1;
                x = 0;
            elseif df(i,3) == k
                hist(k,x) = df(i,5);
            end
        % theta[2] Case
        elseif j == 2
            if mod(i,15000) == 0 && i > 1
                hist(k,x) = df1(i,5);
                k = k + 1;
                x = 0;
            elseif df1(i,3) == k
                hist(k,x) = df1(i,5);
            end
        end
        
        % Updating Counter
        x = x + 1;
    
    end
    
    % Plotting Results
    subplot(1,2,j)
    
    % Plotting Density Plot for 5 Chains
    for i = 1:5
        hold on
        % Plotting histogram
        p(i) = histogram(hist(i,:),300,'FaceAlpha',0.5,'EdgeColor','none', ...
            'Normalization', 'probability','FaceColor',Colors{i});
        % Plotting distribution line
        [values, edges] = histcounts(hist(i,:), 300, 'Normalization', 'probability');
        centers = (edges(1:end-1)+edges(2:end))/2;
        b = plot(centers, values, '-','Color',Colors{i},'LineWidth',1.8);
    end

end

% Plot Aesthetics

% theta[1] Plot

sp1 = subplot(1,2,1); hold on; grid minor; box on;

% Legend and Axis
leg1 = legend(p,'1','2','3','4','5'); title(leg1,'Chain');
ylabel('Probability','FontSize',12)
title(sp1,'theta[1]')
xlim([0.65 1.45]); ylim([0 inf])

% theta[2] Plot

sp2 = subplot(1,2,2); hold on; grid minor; box on;

% Legend and Axis
leg2 = legend(p,'1','2','3','4','5'); title(leg2,'Chain');
title(sp2,'theta[2]')
xlim([0.6 3.3*10^4]); ylim([0 inf])
ax=gca; ax.XAxis.Exponent = 0;

% Figure Aesthetics

spax=axes(fig,'visible','off'); 
spax.XLabel.Visible='on'; spax.YLabel.Visible='on';
xlabel(spax,'Parameter Value'); % ylabel(spax,'Probability');
sgtitle('MCMC Posterior Distributions')

% Export Settings
set(gcf, 'units', 'centimeters', 'position', [1 1 12.1 9]);
% export_fig('MCMC_density_plots_theta', gcf, '-transparent', '-png', '-r500')