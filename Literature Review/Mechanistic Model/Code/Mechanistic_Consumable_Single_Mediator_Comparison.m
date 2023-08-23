% Single (chemical - consumable) mediator interaction - Momeni. et al. (2018) - Figure 1.
% Connor Tynan - 25/10/22 - Final year project

clc; clear; clf; close all;

% Background

% https://doi.org/10.7554/eLife.25051

% S1 -> C1 -> S2
% S1 produces C1, C1 is then consumed by S2, promoting S2 growth
% C1 production rate proportional to S1 population density
% As S1 -> 0, expect C1 to begin decreasing.

% Parameters

% Time
t0 = 0;     % Initial time
t1 = 3000;  % Final time
dt = 0.001;   % Time step size

% Fitness (Evolution-controlling variables)
r0 = [-0.01 0.001]; % Basal growth rate (Per-species growth rate - constant)

% Interaction terms - please see literature for details
Kc1s2 = 1e5;  % S2 consuming C1 at 1/2 max rate - Carrying capacity
Ks2c1 = 1e5;  % C1 affecting S2 at 1/2 max rate - Carrying capacity

rs2c1 = 0.1;  % C1 affecting S2 at max rate

% Misc
a = 1; % Consumption rate
b = 0.05; % Production rate

% Variables

time = t0:dt:t1; % Time vector (x-axis)

% Initialising solution vectors
S1 = zeros(1,length(time)); % Population density of species S1
C1 = zeros(1,length(time)); % Concentration of chemical C1
S2 = zeros(1,length(time)); % Population density of species S2

S1LV = zeros(1,length(time)); % Population density of species S1 (LV)
S2LV = zeros(1,length(time)); % Population density of species S2 (LV)

% Initial Conditions

S1(1) = 20;
C1(1) = 0;
S2(1) = 10;

S1LV(1) = 20;
S2LV(1) = 10;

% Main computation

for t=1:(numel(time)-1)
    
    % Forward solving mechanistic equations
    S1(t+1) = S1(t) + dt * S1(t) * r0(1);
    C1(t+1) = C1(t) + dt * (b*S1(t) - S2(t)*a*C1(t)/(C1(t)+Kc1s2));
    S2(t+1) = S2(t) + dt * S2(t) * ( r0(2) + rs2c1*C1(t)/(C1(t)+Ks2c1) );
    
    % Forward solving LV equations
    S1LV(t+1) = S1LV(t) + dt * r0(1) * S1LV(t);
    S2LV(t+1) = S2LV(t) + dt * S2LV(t) * (r0(2) + rs2c1 * ( S1LV(t))/(S1LV(t) + Ks2c1) );

end

% Plots

grid minor; hold on; grid; box on;

% Plotting mechanistic species densities
plot(time,S1,'-','LineWidth',2,'color',[0.2 0.2 0.9])
plot(time,S2,'-','LineWidth',2,'color',[0.9 0.2 0.2])

% Plotting chemical concentration
plot(time,C1,'-.','LineWidth',2,'color',[1 0.7 0.2])

% Plotting LV species densities
plot(time,S1LV, '--','LineWidth',2,'color',[0.2 0.2 0.9]); 
plot(time,S2LV, '--','LineWidth',2,'color',[0.9 0.2 0.2])

% Plot aesthetics

legend('S1', 'S2', 'C1','S1 LV','S2 LV','Location','northwest')

title('Single Metabolite Modelling - Mechanistic vs. gLV')
xlabel('Time [hr]','FontSize',12)
ylabel('Magnitude')

% ax = gca;
% ax.XAxis.FontSize = 13;
% ax.YAxis.FontSize = 13;
set(gcf, 'units', 'centimeters', 'position', [1 1 12.1 9]);
