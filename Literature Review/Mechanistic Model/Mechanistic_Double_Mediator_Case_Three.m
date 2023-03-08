% Dual (chemical - consumable) mediator interaction - Momeni. et al. (2018) - Figure 6.
% Connor Tynan - 25/10/22 - Final year project

% Background

% https://doi.org/10.7554/eLife.25051

% Increasing the number of mediators increases the complexity of the system
% exponentially. More-involved systems such as this is where the gLV model
% typically falls down as it fails to account for cases of increased
% complexity. In addition, the gLV model fails to account for indirect
% influences, suggesting use of mechanistic models for higher
% order/complexity systems.

clc; clf; close all; clear

% Parameters

t0 = 0;     % Initial time
t1 = 3000;   % Final time
dt = 0.01;   % Time step size

r0 = [-0.02, -0.01]; % Basal growth rate (Per-species growth rate - constant)

% Interaction terms - please see literature for details

% Kij - Mediator Cj affecting species Si at 50% of max fitness effect - Kij ~ r0.5ij 
Ks1c1 = 1.1 * 1e5;
Ks1c2 = 0;
Ks2c2 = 1e5;
Ks2c1 = 2 * 1e5;

%rij - Mediator Cj affecting species Si at maximum fitness effect
rs1c1 = 15;
rs1c2 = 0;
rs2c1 = 15;
rs2c2 = 25;  
  
% Consumption rate of chemical mediators
a = 1; 
ac1s1 = a;
ac1s2 = a;
ac2s2 = a;

% Production rate of chemical mediators 
b = 0.1;
b0 = 10*b;
bc2s1 = b*10;

% Variables

time = t0:dt:t1; % Time vector (x-axis)

% Initialising solution vectors
S1 = zeros(1,length(time)); % Population density of species S1
S2 = zeros(1,length(time)); % Population density of species S2
C1 = zeros(1,length(time)); % Concentration of chemical C1
C2 = zeros(1,length(time)); % Concentration of chemical C2

% Initial Conditions

% Initial population densities
S1(1) = 5;
S2(1) = 5;

% Initial mediator concentrations
C1(1) = 0.00;
C2(1) = 0.00001;

% Main computation - Equation (47) in the literature
% Forward solving mechanistic equations

for t=1:(numel(time)-1)
    
    % Species densities

    S1(t+1) = S1(t) + dt * S1(t) * ( r0(1) + rs1c1*C1(t)/(C1(t)+Ks1c1));

    S2(t+1) = S2(t) + dt * S2(t) * ( r0(2) + rs2c2 * ...
        (((C1(t)*C2(t))/(Ks2c1*Ks2c2))/(C1(t)/Ks2c1 + C2(t)/Ks2c2)) * ...
        ((1/(C1(t)/Ks2c1 + 1)) + (1/(C2(t)/Ks2c2 + 1))));

    % Chemical concentrations

    C1(t+1) = C1(t) + dt * (b0 - ac1s1 * rs1c1 * C1(t)/(C1(t)+Ks1c1) - ac1s2 * ...
        rs2c2* (((C1(t)*C2(t))/(Ks2c1*Ks2c2))/(C1(t)/Ks2c1 + C2(t)/Ks2c2)) * ...
        ((1/(C1(t)/Ks2c1 + 1)) + (1/(C2(t)/Ks2c2 + 1))) * S2(t) );

    C2(t+1) = C2(t) + dt * (bc2s1 * S1(t) - ac2s2 * rs2c2 * ...
        (((C1(t)*C2(t))/(Ks2c1*Ks2c2))/(C1(t)/Ks2c1 + C2(t)/Ks2c2)) * ...
        ((1/(C1(t)/Ks2c1 + 1)) + (1/(C2(t)/Ks2c2 + 1))) * S2(t) ) ;
    
end

% Plotting results

% Plotting species densities

fig = figure(1); hold on; box on

plot(time,S1./(S1+S2),'-b','LineWidth',2); plot(time,S2./(S1+S2),'Color',[1, 0.5, 0.1],'LineWidth',2)

ylabel('Fractions of \color{blue}S1\color{black} and\color{orange} S2', ...
    'interpreter', 'tex')
ylim([10^-4 10^0+0.1])
yticks([10^-4 10^-2 10^0])

xlabel('Time (hrs)')
xlim([0 3000])
xticks([0 1000 2000 3000])

title('Mechanistic Model - Case III')

set(gca,'YScale','log')
set(gca,'XMinorTick','off','YMinorTick','off')

set(gcf, 'units', 'centimeters', 'position', [1 1 12.1 9]);

% Plotting chemical concentrations

fig2 = figure(2); hold on; box on

p = plot(time,C1,'-k','LineWidth',2); plot(time,C2,'Color',[1, 0.0, 0.0],'LineWidth',2)

ylabel('[C1] and \color{red}{[C2]} \color{black} concentrations (µM)', ...
    'interpreter', 'tex')
ylim([0 10^6+1000000])
yticks([0 10^0 10^2 10^4])

xlabel('Time (hrs)')
xlim([0 3000])
xticks([0 1000 2000 3000])

title('Mechanistic Model - Case III')

set(gca,'YScale','log')
set(gca,'XMinorTick','off','YMinorTick','off')

set(gcf, 'units', 'centimeters', 'position', [1 1 12.1 9]);