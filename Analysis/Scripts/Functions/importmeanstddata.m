function [meansh1std,meansh2std,meansh3std,meansh4std] = importmeanstddata

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

meansh1std = zeros(14,1); meansh2std = zeros(14,1); 
meansh3std = zeros(14,1); meansh4std = zeros(14,1); 

% Cleaning our new averages - Omitting N/A & NANs

for i = 1:14
    meansh1std(i) = std(rmmissing(sh1(i,2:7)));
    meansh2std(i) = std(rmmissing(sh2(i,2:7)));
    meansh3std(i) = std(rmmissing(sh3(i,2:7)));
    meansh4std(i) = std(rmmissing(sh4(i,2:7)));
end

end