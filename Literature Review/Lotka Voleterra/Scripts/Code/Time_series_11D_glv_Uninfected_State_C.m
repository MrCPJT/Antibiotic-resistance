% 11D gLV model - Fitting 11D gLV model to mouse gut microbiome data
% Connor Tynan - 17/10/22 - Final year project

% Datasets from Stein. et al. (2013) 

% Data includes quantified (inferred) growth and species-species interactions 
% for 11 individual microbial species.

% Experimental measures of population size with time are also included, 
% initial conditions for our model can be taken from here by taking
% population sizes of individuals at time t = 0. In total there are 9
% initial conditions to use from the data.

% Full dataset can be found here: https://doi.org/10.1371/journal.pcbi.1003388.s001
% Similarly data can be found here: https://github.com/erijones/ssr_module
% under (stein_parameters.csv) & (stein_ic.csv)

clc; clf; clear; close all;

% Importing and parsing data

df = readtable('C:\Users\Tynan\Desktop\gLV 11D\stein_parameters.csv'); % Importing table
ro = table2array(df(1:11,13));       % Per-organism growth rate
sus = table2array(df(1:11,14));      % Antibiotic susceptibility
K = table2array(df(1:11,2:12));      % Interactions between microbes - K[i][j] represents effect of pop. j on pop. i

% Importing and parsing initial conditions
% Initial conditions from 9 experiments

df2 = readtable('C:\Users\Tynan\Desktop\gLV 11D\stein_ic.csv'); % Importing table
init_ind = [2 13 17 27 38 42 53 64 68]; % Column index at which t=0 (ics)
ics = zeros(length(init_ind),11);                               % Initialise array for ics
for i = 1:length(init_ind)      
    ics(i,:) = table2array(df2(6:16,init_ind(i))); % Saving the ics as a matrix 
end

% Modification of initial condition based on outcome

ics2(1,:) = ics(2,:); ics2(2,:) = ics(5,:); ics2(3,:) = ics(8,:);

% Parameters

t0 = 0;     % Initial time
t1 = 200;   % Final time
dt = 0.01;  % Time step size

% Variables

te = t0:dt:t1;               % Time vector for Euler method
y = zeros(11,length(te));   % Initialising solution vector
s = zeros(size(K));         % Vector to hold K_{ij}*y(j) values
sigma = zeros(11,length(te));        % Vector to hold summation values
titles = ["Barnesiella", "und. Lachnospiraceae", "uncl. Lachnospiraceae", "Other"... 
    "Blautia" "und. Mollicutes" "Akkermansia" "Coprobacillus" "Clostridium Difficile"...
    "Enterococcus" "und. Enterobacteriaceae"];

fig = figure(1);

% Main Computation

% Looping over different initial conditions
for k = 1:3

    % Initialising Initial conditions
    y(:,1) = ics2(k,:);
    
    % Time series Computation - Forward Euler
    for t = 1:(numel(te)-1)
        
        % Looping over different species
        for i = 1:11        
            
            % Summation notation
            for j = 1:11    
                s(i,j) = K(i,j).* y(j,t);
            end
    
            % Calculating Sigma
            sigma(i,t) = sum(s(i,:)); 
    
            % Forward Euler
            y(i,t+1) = y(i,t) + dt*( y(i,t)*( ro(i) + sigma(i,t) ) );
            
        end
    
    end
    
    % Plotting the results
    for i = 1:11
        if i < 9 
            sp = subplot(3,4,i); hold on; 
            plot(te, y(i,:),'LineWidth',2)
            ylim([-0.5 15]) 
            title(titles(i));
            if k == 3
                grid minor; box on
            end
        elseif i >= 9
            sp = subplot(3,3,i-2); hold on;
            plot(te, y(i,:),'LineWidth',2)
            ylim([-0.5 15])
            title(titles(i));
            if k == 3
                grid minor; box on
            end
        end
    end
end

% Plot Aesthetics / Miscellaneous

spax=axes(fig,'visible','off'); 
spax.XLabel.Visible='on'; spax.YLabel.Visible='on';
xlabel(spax,'Time [hr]','FontSize',14); ylabel(spax,'Population Growth Rate','FontSize',14);
sgtitle('Species Time Series - Uninfected C Case')

set(gcf, 'units', 'centimeters', 'position', [1 1 24 16]);