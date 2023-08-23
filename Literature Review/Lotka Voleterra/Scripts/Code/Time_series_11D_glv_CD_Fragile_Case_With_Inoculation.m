% 11D gLV model - Fitting 11D gLV model to mouse gut microbiome data
% Connor Tynan - 17/10/22 - Final year project

clc; clear; clf; close;

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


% Importing and parsing growth/interaction data

df = readtable('C:\Users\Tynan\Desktop\gLV 11D\stein_parameters.csv'); clc; % Importing data
ro = table2array(df(1:11,13));       % Per-organism growth rate
sus = table2array(df(1:11,14));      % Antibiotic susceptibility
K = table2array(df(1:11,2:12));      % Interactions between microbes
                                     % K[i][j] represents effect of pop. j on pop. i

% Importing and parsing initial conditions 

df2 = readtable('C:\Users\Tynan\Desktop\gLV 11D\stein_ic.csv'); clc; % Importing experimental data
init_ind = [2 13 17 27 38 42 53 64 68]; % Column index at which t=0 (initial population size)
ics = zeros(length(init_ind),11);       % Initialise array for initial conditions

for i = 1:length(init_ind)      
    ics(i,:) = table2array(df2(6:16,init_ind(i))); % Saving initial conditions in a matrix 
end


% Parameters

t1 = 400; % Final time
dt = 0.01; % Time step size


% Variables

te = 0:dt:t1;                 % Time vector for Euler method

% Regular gLV system
y = zeros(11,length(te));     % Initialising solution vector
s = zeros(size(K));           % Vector to hold K_{ij}*y(j) values
sigma = zeros(11,length(te)); % Vector to hold summation values

% Perturbed System
yp = zeros(11,length(te));     % Initialising solution vector
sp = zeros(size(K));           % Vector to hold K_{ij}*y(j) values
sigmap = zeros(11,length(te)); % Vector to hold summation values

% Initial conditions

y(:,1) = ics(5,:);       % Initialising solution vector with initial conditions
y_acum(:,1) = ics(5,:);  % Initialising solution vector with initial conditions
yp(:,1) = ics(5,:);      % Initialising solution vector for perturbed system with initial conditions

u = zeros((numel(te)-1),1); % Initialising antibiotic dosage with initial conditions
u(5000) = 30;                   % Dosage at t = 1

% Main computation

% Time series
for t = 1:(numel(te)-1)
    
    % Looping over different species
    for i = 1:11        
        
        % Summation notation
        for j = 1:11    
            s(i,j) = K(i,j).* y(j,t);
            sp(i,j) = K(i,j).* yp(j,t);
        end

        % Calculating Sigma
        sigma(i,t) = sum(s(i,:)); 
        sigmap(i,t) = sum(sp(i,:)); 

        % Forward Euler
        y(i,t+1) = y(i,t) + dt*( y(i,t)*( ro(i) + sigma(i,t) ));
        yp(i,t+1) = yp(i,t) + dt*( yp(i,t)*( ro(i) + sigmap(i,t) + sus(i)*u(t)) );
        
        if t == 20000
            y(9,t) = y(9,t) + 1e-10;
            yp(9,t) = y(9,t) + 1e-10;           
        end

    end

end


% Plotting the results

fig = figure(1); hold on; grid minor; box on

% CD inoculation
plot([200 200], [0 26],':','color',[0.8 0 0],'LineWidth',2)

% Time series
plot(te,sum(y),'LineWidth',2,'color',[0.4 0.3 0.9]);
plot(te,sum(yp),'--k','LineWidth',2)

% Patching CD inoculation
patch([400 200 200 400], [0 0, 26 26], [0.8 0.8 0.8], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0.6 0 0]);
patch([200 0 0 200], [0 0, 26 26], [0.8 0.8 0.8], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', [0 0.5 0]);

% Plot aesthetics

% Labels
ylabel('Total Population Growth Rate');
xlabel('Time [hr]');

% Title
title('11D gLV Time Series - CD-fragile','FontSize',14)

% Limits
ylim([0 26])

% Legend
leg = legend('', 'gLV', 'gLV + Perturbation','Location','east','FontSize',12);
leg.ItemTokenSize = [21,9];
title(leg,'Model')

% Annotation
text(300, 8.5, 'Infected D','color', 'k','FontSize',13)
text(80, 5, 'Uninfected E','color', 'k','FontSize',13)
text(280, 23.5, 'Uninfected C','color', 'k','FontSize',13)
text(205, 12, 'CD Inoculation','color', [0.4 0.4 0.4],'FontSize',13)

ax = gca;
ax.XAxis.FontSize = 13;
ax.YAxis.FontSize = 13;
set(gcf, 'units', 'centimeters', 'position', [1 1 12.1 9]);
