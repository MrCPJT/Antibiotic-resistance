% Comparison of two models, gLV and pairwise, in modelling 2D interactions.
% Connor Tynan - 03/11/22 - Final year project

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

% Current script:
% 
% Want to compare the evolution of species growing together. So, consider
% systems of 2 species and solve simultaneously. Particularly interested in
% modelling the evolution of species in the presence of C.diff (pathogenic
% bacterium). In addition, choose to model only 'identified'
% species, unidentified species are less useful in research.

% Forward Euler Numerical Scheme

% Current script compares cases of the pairwise model from Momeni. et al.
% and the gLV model

clc; clf; close all; clear

% Importing and parsing growth/interaction data

df = readtable('C:\Users\Tynan\Desktop\gLV 11D\stein_parameters.csv'); clc; % Importing data
ro = table2array(df(1:11,13));       % Per-organism growth rate
K = table2array(df(1:11,2:12));      % Interactions between microbes
                                     % K[i][j] represents effect of pop. j on pop. i

% Importing and parsing experimental results / initial conditions 

df2 = readtable('C:\Users\Tynan\Desktop\gLV 11D\stein_ic.csv'); clc; % Importing the data
init_ind = [2 13 17 27 38 42 53 64 68]; % Index at which t = 0 (initial population size)

ics = zeros(length(init_ind),11);       % Initialise array for initial conditions

% Storing initial conditions in a matrix
for i = 1:length(init_ind)      
    ics(i,:) = table2array(df2(6:16,init_ind(i))); 
end

% Creating matrices for 'identified' microbes - justfies a comparison to C.diff

% Filtering out unidentified interactivity parameters

unknown_index = [2,3,4,6,11];   % Index of unknown species relative to dataframe
K_k = K;                        % K_k = K_known
ro_k = ro;                      % ro_k = ro_known
k = 0;                          % Counter

% Removing unknown species, updating index values per iteration
for i = unknown_index
    i = i - k;      % Updating index of unknown parameters
    K_k(i,:) = [];  % Removing unknown row entries
    K_k(:,i) = [];  % Removing unknown column entries
    ro_k(i,:) = []; % Removing organism-specific growth rates of unknown species

    k = k + 1;      % Updating counter
end

% Initial conditions for known species only
ics_k = ics; 
k = 0;  % Counter
for i = unknown_index
    i = i - k;
    ics_k(:,i) = [];
    k = k + 1;
end

ncm = [1,2,3,4,6]; % Non-CD indices for 'known matrices'

% Names of identified microbes
m_names = [" Barnesiella" " Blautia" " Akkermansia" " Coprobacillus" " Clostridium Difficile"...
    " Enterococcus"]; 

% Parameters

t0 = 0;     % Initial time
t1 = 600;  % Final time
dt = 0.001; % Time step size

K21 = 3;            % Carrying capacity of the pairwise system

r0 = [0.065 0.045]; % Pairwise model per-species growth rates

% Variables

te = t0:dt:t1; % Time vector for Euler method

% Comparing different species to C.diff: Initialise one vector for species
% solutions and another to contain the relevant C.diff evolutions

% 2D gLV
yk = zeros(5,length(te));   % Initialising solution vector for different species   
ykcd = zeros(5,length(te)); % Initialising solution vector for C.diff

% Pairwise Model
ykp = zeros(5,length(te));   % Initialising solution vector for different species   
ykcdp = zeros(5,length(te)); % Initialising solution vector for C.diff

% Initial conditions

% Majority of initial conditions from Stein. et al. contain no C.diff (exception being IC8)
% and require inoculation. To remedy this, initialise C.diff densities randomly
% with arbitrarily small values comparable to the one initial condition containing C.diff 

for i = 1:5
    % 2D gLV model
    yk(i,1) = ics_k(1,ncm(i)); % Initialising solution vector with initial conditions
    ykcd(i,1) = 1e-5 + (1e-4 - 1e-5) * rand(1,1);
    
    % Pairwise model
    ykp(i,1) = ics_k(1,ncm(i)); % Initialising solution vector with initial conditions
    ykcdp(i,1) = 1e-5 + (1e-4 - 1e-5) * rand(1,1);
end

% Main computation

% Solving relevant 2D gLV equations to model the time-series evolution of
% two comparing species. Each pair of species has unique interactions.

for i = 1:5 % Updating microbial species to be compared and instance of C.diff evolution
    for t = 1:(numel(te)-1)
    
        % 2D gLV model
        % Stepping forward non-C.diff species 
        yk(i,t+1) = yk(i,t) + dt * ( yk(i,t)*( ro_k(ncm(i)) + ykcd(i,t) * K_k(ncm(i),5) + ...
            yk(i,t) * K_k(ncm(i),ncm(i))));  

        % Stepping forward C.diff species
        ykcd(i,t+1) = ykcd(i,t) + dt * ( ykcd(i,t)*( ro_k(5) + yk(i,t) * K_k(5,ncm(i)) + ...
            ykcd(i,t) * K_k(5,5)));         
        
        % Pairwise model
        % Stepping forward non-C.diff species
        ykp(i,t+1) = ykp(i,t) + dt * ykp(i,t) * ( r0(1) + ykp(i,t) * K_k(ncm(i),ncm(i))/(ykp(i,t) + K21) + ...
            ykcdp(i,t) * K_k(ncm(i),5)/(ykcdp(i,t) + K21));

        % Stepping forward C.diff species
        ykcdp(i,t+1) = ykcdp(i,t) + dt * ykcdp(i,t) * ( r0(2) + ykp(i,t) * ...
            K_k(5,ncm(i))/(ykp(i,t) + K21) + ykcdp(i,t) * K_k(5,5)/(ykcdp(i,t) + K21));

    end
end

% Plotting Results

colours = {[0 0.6 0.2], [0.7 0.5 0.2], [0.3 0.2 0.9], [0.8 0 0.8], [0.2 0.6 0.6]};

fig = figure(1); 

for i = 1:6
    % Plotting results
    if i < 6
        subplot(2,3,i); hold on; box on; grid minor
        p{i} = plot(te,yk(i,:),'Color',colours{i},'LineWidth',2,'DisplayName',m_names(ncm(i)));
        plot(te,ykp(i,:),'--','Color',colours{i},'LineWidth',2)
        cdp = plot(te,ykcd(i,:),'-r','LineWidth',2,'DisplayName',' C. Difficile');
        plot(te,ykcdp(i,:),'--r','LineWidth',2);
        dashed = plot([0 0],[0 0],'--k','LineWidth',2,'DisplayName',' Pairwise');
        solid = plot([0 0],[0 0],'-k','LineWidth',2,'DisplayName',' gLV');
    % Plotting Legends
    elseif i == 6
        fake_ax = subplot(2,3,i);
        ax_pos = fake_ax.Position;
        area([NaN NaN], NaN(2, 4));
        
        hl = legend([p{1:5},cdp]); % hl.Layout.Tile = "east";
        
        leg_pos = hl.Position;
        delete(fake_ax);
        ax = axes('Position', [ax_pos(1:2) 0 0]);
        leg = legend([p{1:5},cdp],'FontSize',12); % hl.Layout.Tile = "east";
        leg.Location = 'none';
        leg.Position = [(0.8 - leg_pos(3) / 2) (0.25 - leg_pos(4) / 2) leg_pos(3:4)];
        leg.ItemTokenSize = [12,9];
        title(leg,'Known Species')
        ax.Visible = false;
        
        fake_ax2 = subplot(2,3,i);
        ax_pos2 = fake_ax2.Position;
        area([NaN NaN], NaN(2, 4));
        
        hl2 = legend([dashed,solid]);
        
        leg_pos2 = hl2.Position;
        delete(fake_ax2);
        ax = axes('Position', [ax_pos2(1:2) 0 0]);
        leg2 = legend([dashed,solid],'FontSize',12); % hl.Layout.Tile = "east";
        leg2.Location = 'none';
        leg2.Position = [(0.8 - leg_pos2(3) / 2) (0.25 - leg_pos2(4) / 2) leg_pos2(3:4)];
        leg2.ItemTokenSize = [21,9];
        title(leg2,'Model')
        ax.Visible = false;

    end
end

spax=axes(fig,'visible','off'); 
spax.XLabel.Visible='on'; spax.YLabel.Visible='on';
xlabel(spax,'Time [hr]','FontSize',12); ylabel(spax,'Population Growth Rate','FontSize',12);
sgtitle('Modelling 2D interactions: Pairwise vs. gLV')

set(gcf, 'units', 'centimeters', 'position', [1 1 24 12]);
% export_fig('Modelling_2_species_glv_vs_pairwise', gcf, '-transparent', '-png', '-r500')