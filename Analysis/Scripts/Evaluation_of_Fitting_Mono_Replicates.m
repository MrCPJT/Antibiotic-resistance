% Connor Tynan - Fitting estimates to gLV model

clc; clear; clf; close all; format long

% Importing the data - See importfile.m function file for details

sheet1 = importfile("C:\Users\Tynan\Desktop\Maths\Year 4\Dissertation\Parameter Fitting\Data\221108-Staph-Pseudomonas interaction data.xlsx","S. aureus mono-culture","A3:G16");
sheet2 = importfile("C:\Users\Tynan\Desktop\Maths\Year 4\Dissertation\Parameter Fitting\Data\221108-Staph-Pseudomonas interaction data.xlsx","P. aeruginosa mono-culture","A3:G16");
sheet3 = importfile("C:\Users\Tynan\Desktop\Maths\Year 4\Dissertation\Parameter Fitting\Data\221108-Staph-Pseudomonas interaction data.xlsx","S. aureus co-culture","A3:G16");
sheet4 = importfile("C:\Users\Tynan\Desktop\Maths\Year 4\Dissertation\Parameter Fitting\Data\221108-Staph-Pseudomonas interaction data.xlsx","P. aeruginosa co-culture","A3:G16");

% Creating matrices to work with from data

sh1 = table2array(sheet1);
sh2 = table2array(sheet2);
sh3 = table2array(sheet3);
sh4 = table2array(sheet4);

% Housekeeping

clear sheet1 sheet2 sheet3 sheet4

% Averaging %

% 1D growth averages - Averaging data for mono-culture species growth

meansh1 = zeros(14,1); meansh2 = zeros(14,1); 
meansh3 = zeros(14,1); meansh4 = zeros(14,1); 

% Cleaning our new averages - Omitting N/A & NANs

for i = 1:14
    meansh1(i) = mean(rmmissing(sh1(i,2:7)));
    meansh2(i) = mean(rmmissing(sh2(i,2:7)));
    meansh3(i) = mean(rmmissing(sh3(i,2:7)));
    meansh4(i) = mean(rmmissing(sh4(i,2:7)));
end

% Parameters

% Basic time and step parameters

t0 = 0; t1 = 15; dt = 0.001; te = t0:dt:t1;

% Variables

mu1 = [1 1.053333333 0.91 0.84 1.02 1.3 1.07 1.18]; 
M11 = [17600 20626.5 28900 41900 14210 19000 11509 8240];
mu2 = [0.82 0.895714286 0.92 0.96 0.84 0.85 1.1 0.78]; 
M22 = [19800 12314.28571 19200 12700 9200 3300 10000 12000];

for k = 1:length(mu1)

    % Initial conditions
    
    y = zeros(t1/dt,1); y(1) = (1);  % S. Aureus population
    z = zeros(t1/dt,1); z(1) = (1);  % P. Aeruginosa population
    
    for i = 1:(numel(te)-1)
    
        y(i+1) = y(i) + dt * ( y(i)*mu1(k) - 1/(M11(k)) * y(i).^2 );
        z(i+1) = z(i) + dt * ( z(i)*mu2(k) - 1/(M22(k)) * z(i).^2 );
        
    end

    if k == 1

        figure(1); hold on
        plot(t0:1:13,meansh1/meansh1(1),'ok','MarkerSize',7)  % Experimental data
        plot(te,y,'LineWidth',1.2,'Color',[0.8 0.5 0.5]) % Parameter estimation 

        figure(2); hold on
        plot(t0:1:13,meansh2/meansh2(1),'ok','MarkerSize',7)  % Experimental data
        plot(te,z,'LineWidth',1.2,'Color',[0.6 0.5 0.7]) % Parameter estimation

    end

    if k == 2

        figure(1); hold on; grid minor; grid on; box on;
        plot(t0:1:13,meansh1/meansh1(1),'ok','MarkerSize',7)  % Experimental data
        plot(te,y,'--','LineWidth',1.2,'Color',[0.8 0.5 0.5]) % Parameter estimation 
        
        figure(2); hold on; grid minor; grid on; box on;
        plot(t0:1:13,meansh2/meansh2(1),'ok','MarkerSize',7)  % Experimental data
        plot(te,z,'--','LineWidth',1.2,'Color',[0.6 0.5 0.7]) % Parameter estimation 
    end

    if k > 2

        figure(3)
        subplot(2,3,k-2); hold on; grid minor; grid on; box on;
        plot(t0:1:13,sh1(:,k-1)/sh1(1,k-1),'ok')  % Experimental data
        plot(te,y,'LineWidth',1.2,'Color',[0.8 0.5 0.5])                  % Parameter estimation 

        figure(4)
        subplot(2,3,k-2); hold on; grid minor; grid on; box on;
        plot(t0:1:13,sh2(:,k-1)/sh2(1,k-1),'ok')  % Experimental data
        plot(te,z,'LineWidth',1.2,'Color',[0.6 0.5 0.7])                  % Parameter estimation 
    end

end

% Plots / Plotting

% Replicate plots

figure(1)
ylim([0 max(meansh1/meansh1(1))*1.1])

ylabel('Normalised Population Size','interpreter','tex','FontWeight','bold','FontSize',14);
xlabel('Time [s]','interpreter','tex','FontWeight','bold','FontSize',14);

legend('','Mean-averaged','','Replicate-averaged','FontSize',15,'Location','NorthWest')

title('S. Aureus - Parameter estimates','FontSize',16)

figure(2)
ylim([0 max(meansh2/meansh2(1))*1.1])

ylabel('Normalised Population Size','interpreter','tex','FontWeight','bold','FontSize',14);
xlabel('Time [s]','interpreter','tex','FontWeight','bold','FontSize',14);

legend('','Mean-averaged','','Replicate-averaged','FontSize',15,'Location','NorthWest')

title('P. Aeruginosa - Parameter estimates','FontSize',16)

figure(3); 

ylim([0 35000])

subfig = axes(gcf,'visible','off');
subfig.Title.Visible='on';
subfig.XLabel.Visible='on';
subfig.YLabel.Visible='on';

ylabel(subfig,'Normalised Population Size','interpreter','tex','FontWeight','bold','FontSize',20,'Position',[-0.05 0.5]);
xlabel(subfig,'Time [s]','interpreter','tex','FontWeight','bold','FontSize',20);

sgtitle('S. Aureus - Parameter Estimates vs. Observation Data','FontSize',24)

figure(4); 

ylim([0 25000])

subfig = axes(gcf,'visible','off');
subfig.Title.Visible='on';
subfig.XLabel.Visible='on';
subfig.YLabel.Visible='on';

ylabel(subfig,'Normalised Population Size','interpreter','tex','FontWeight','bold','FontSize',20,'Position',[-0.05 0.5]);
xlabel(subfig,'Time [s]','interpreter','tex','FontWeight','bold','FontSize',20);

sgtitle('P. Aeruginosa - Parameter Estimates vs. Observation Data','FontSize',24)