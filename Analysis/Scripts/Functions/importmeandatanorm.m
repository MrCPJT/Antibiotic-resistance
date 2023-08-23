function [meansh1norm,meansh2norm,meansh3norm,meansh4norm] = importmeandatanorm

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

meansh1norm = meansh1/meansh1(1);
meansh2norm = meansh2/meansh2(1);
meansh3norm = meansh3/meansh3(1);
meansh4norm = meansh4/meansh4(1);

end