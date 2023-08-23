% Connor Tynan - Visualising S.aureus and P.aerugniosa mono and co-culture
% data.

% NOTE: Requires tight_subplot.m and barwitherr.m functions from the public
% mathworks file exchange

% Miscellaneous
clc; clf; close all; clear

% Importing the data
[mean_mono_sa, mean_mono_pa, mean_co_sa, mean_co_pa] = importmeandata;
[std_mono_sa, std_mono_pa, std_co_sa, std_co_pa] = importmeanstddata;
[mono_sa, mono_pa, co_sa, co_pa] = importrawdata;

% Figure 1 - Mean averaged with error bars
fig = figure(1);

% Initialising 2x2 subplot
set(gca,'Visible','off')
ha = tight_subplot(2,2,[.125 .065],[.075 .075],[.05 .05]);

% S.aureus mono-culture
axes(ha(1));
% Bar chart with error bars
barwitherr(cat(3,zeros(14,1),std_mono_sa),0:1:13,mean_mono_sa,'linewidth',2,'FaceColor','#1E88E5');
% Plot aesthetics
xlabel('Time (hr)'); ylabel('Bacteria Count (ml^{-1})')
LimitsX = xlim; LimitsY = ylim;
title('S.aureus Mono-Culture','HorizontAlalignment', 'right', 'position', [LimitsX(2), LimitsY(2)])
set(gca,'FontSize',12)

% P.aeruginosa mono-culture
axes(ha(2));
% Bar chart with error bars
barwitherr(cat(3,zeros(14,1),std_mono_pa),0:1:13,mean_mono_pa,'linewidth',2,'FaceColor','#FFC107');
% Plot aesthetics
xlabel('Time (hr)'); ylabel('Bacteria Count (ml^{-1})')
LimitsX = xlim; LimitsY = ylim;
title('P.aeruginosa Mono-Culture','HorizontAlalignment', 'right', 'position', [LimitsX(2), LimitsY(2)])
set(gca,'FontSize',12)

% S.aureus co-culture
axes(ha(3));
% Bar chart with error bars
barwitherr(cat(3,zeros(14,1),std_co_sa),0:1:13,mean_co_sa,'linewidth',2,'FaceColor','#004D40');
% Plot aesthetics
xlabel('Time (hr)'); ylabel('Bacteria Count (ml^{-1})')
LimitsX = xlim; LimitsY = ylim;
title('S.aureus Co-Culture','HorizontAlalignment', 'right', 'position', [LimitsX(2), LimitsY(2)])
set(gca,'FontSize',12)

% P.aeruginosa co-culture
axes(ha(4));
% Bar chart with error bars
barwitherr(cat(3,zeros(14,1),std_co_pa),0:1:13,mean_co_pa,'linewidth',2,'FaceColor','#D81B60');
% Plot aesthetics
xlabel('Time (hr)'); ylabel('Bacteria Count (ml^{-1})')
LimitsX = xlim; LimitsY = ylim;
title('P.aeruginosa Co-Culture','HorizontAlalignment', 'right', 'position', [LimitsX(2), LimitsY(2)])
set(gca,'FontSize',12)

% Export Settings
set(gcf, 'units', 'centimeters', 'position', [1 1 12.1*2 9*2]);


% Figure 2-5 - Visualising the raw data

% S.aureus mono-culture
figure(2)

% Initialising Subplot
set(gca,'Visible','off')
ha = tight_subplot(2,3,[.085 .065],[.07 .07],[.04 .04]);

% Plotting raw data
for i = 1:6
    axes(ha(i));
    % Bar chart with error bars
    bar(0:1:13,mono_sa(:,i),'LineWidth',2,'FaceColor','#1E88E5');
    axis tight
    % Plot aesthetics
    LimitsX = xlim; LimitsY = ylim;
    title(sprintf('Replica %d',i),'HorizontAlalignment', 'right', 'position', [LimitsX(2), LimitsY(2)])
    xticks([])
    set(gca,'FontSize',14)
end

% P.aeruginosa mono-culture
figure(3)

% Initialising Subplot
set(gca,'Visible','off')
ha = tight_subplot(2,3,[.085 .065],[.07 .07],[.04 .04]);

% Plotting raw data
for i = 1:6
    axes(ha(i));
    % Bar chart with error bars
    bar(0:1:13,mono_pa(:,i),'LineWidth',2,'FaceColor','#FFC107');
    axis tight
    % Plot aesthetics
    LimitsX = xlim; LimitsY = ylim;
    title(sprintf('Replica %d',i),'HorizontAlalignment', 'right', 'position', [LimitsX(2), LimitsY(2)])
    xticks([])
    set(gca,'FontSize',14)
end

% S.aureus co-culture
figure(4)

% Initialising Subplot
set(gca,'Visible','off')
ha = tight_subplot(2,3,[.085 .065],[.07 .07],[.04 .04]);

for i = 1:6
    axes(ha(i));
    % Bar chart with error bars
    bar(0:1:13,co_sa(:,i),'LineWidth',2,'FaceColor','#004D40');
    axis tight
    % Plot aesthetics
    LimitsX = xlim; LimitsY = ylim;
    title(sprintf('Replica %d',i),'HorizontAlalignment', 'right', 'position', [LimitsX(2), LimitsY(2)])
    xticks([])
    set(gca,'FontSize',14)
end

% P.aeruginosa co-culture
figure(5)

% Initialising Subplot
set(gca,'Visible','off')
ha = tight_subplot(2,3,[.085 .065],[.07 .07],[.04 .04]);

% Plotting raw data
for i = 1:6
    axes(ha(i));
    % Bar chart with error bars
    bar(0:1:13,co_pa(:,i),'LineWidth',2,'FaceColor','#D81B60');
    axis tight
    % Plot aesthetics
    LimitsX = xlim; LimitsY = ylim;
    title(sprintf('Replica %d',i),'HorizontAlalignment', 'right', 'position', [LimitsX(2), LimitsY(2)])
    xticks([])
    set(gca,'FontSize',14)
end
