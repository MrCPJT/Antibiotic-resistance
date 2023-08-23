function [sh1,sh2,sh3,sh4] = importrawdata

sheet1 = importfile("C:\Users\Tynan\Desktop\Maths\Year 4\Dissertation\Parameter Fitting\Data\221108-Staph-Pseudomonas interaction data.xlsx","S. aureus mono-culture","A3:G16");
sheet2 = importfile("C:\Users\Tynan\Desktop\Maths\Year 4\Dissertation\Parameter Fitting\Data\221108-Staph-Pseudomonas interaction data.xlsx","P. aeruginosa mono-culture","A3:G16");
sheet3 = importfile("C:\Users\Tynan\Desktop\Maths\Year 4\Dissertation\Parameter Fitting\Data\221108-Staph-Pseudomonas interaction data.xlsx","S. aureus co-culture","A3:G16");
sheet4 = importfile("C:\Users\Tynan\Desktop\Maths\Year 4\Dissertation\Parameter Fitting\Data\221108-Staph-Pseudomonas interaction data.xlsx","P. aeruginosa co-culture","A3:G16");

% Creating matrices to work with from data

sh1 = table2array(sheet1);
sh2 = table2array(sheet2);
sh3 = table2array(sheet3);
sh4 = table2array(sheet4);

sh1 = sh1(:,2:7);
sh2 = sh2(:,2:7);
sh3 = sh3(:,2:7);
sh4 = sh4(:,2:7);

end