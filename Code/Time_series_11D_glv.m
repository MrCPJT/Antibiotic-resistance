% gLV model reproduction - Jones and Carlson - Data from Stein. et al.
% (2013)
% Connor Tynan - 17/10/22 - Dissertation

clc; clf; clear; close all

% Importing and parsing data

df = readtable('C:\Users\Tynan\Desktop\gLV 11D\stein_parameters.csv'); clc; % Importing table
ro = table2array(df(1:11,13));       % Per-organism growth rate
sus = table2array(df(1:11,14));      % Antibiotic susceptibility
K = table2array(df(1:11,2:12));      % Interactions between microbes - K[i][j] represents effect of pop. j on pop. i

% Importing and parsing initial conditions
% Initial conditions from 9 experiments

df2 = readtable('C:\Users\Tynan\Desktop\gLV 11D\stein_ic.csv'); clc; % Importing table
init_ind = [2 13 17 27 38 42 53 64 68]; % Column index at which t=0 (ics)
ics = zeros(length(init_ind),11);                               % Initialise array for ics
for i = 1:length(init_ind)      
    ics(i,:) = table2array(df2(6:16,init_ind(i))); % Saving the ics as a matrix 
end

% Parameters

t1 = 300; % Final time
dt = 0.1; % Time step size

% Variables

te = 0:dt:t1;               % Time vector for Euler method
y = zeros(11,length(te));   % Initialising solution vector
s = zeros(size(K));         % Vector to hold K_{ij}*y(j) values
sigma = zeros(11,length(te));        % Vector to hold summation values
ys = zeros(11,9,1);

% Main computation

fig = figure(1);

for k = 1:9

    % Initial conditions

    y(:,1) = ics(k,:);
    
    % Main computation
    
    % Time series
    for t = 1:(numel(te)-1)
        
        % Different species
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
    
    % Plotting each initial condition
    sp = subplot(3,3,k); hold on; grid minor; box on
    plot(te, y, 'LineWidth', 2)
    title(sp,sprintf('Initial Condition %d',k))

end

spax=axes(fig,'visible','off'); 
spax.XLabel.Visible='on'; spax.YLabel.Visible='on';
xlabel(spax,'Time [hr]','FontSize',14); ylabel(spax,'Population Growth Rate','FontSize',14);
sgtitle('Time Series of 11D-gLV Model')

set(gcf, 'units', 'centimeters', 'position', [1 1 24 16]);
