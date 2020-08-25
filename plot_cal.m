# plot_cal: plots pre-/post- calibration data
#
# (c) 2020 David Goncalves
# MIT License, see LICENSE.txt 

function [A_1, b, t, std_err] = plot_cal(s)
# Generates magnetic calibration parameters, 
# a calibrated data set, and pre-/post- plots

#Set limits and major labels
figure;
xlim([-1,1]);
ylim([-1,1]);
zlim([-1,1]);
xlabel('X');
ylabel('Y');
zlabel('Z');
hold all;

#plot unit sphere and set perspective
[x y z] = sphere;
usph=surf(x, y, z);
set(usph,'facealpha',0.5);
set(usph,'handlevisibility','off');
axis equal;
daspect([1 1 1])
view(30,10)

#plot uncalibrated data
scatter3(s(:,1), s(:,2), s(:,3), [], 'r', 'filled')

#get calibrated data
[A_1, b, t, std_err] = mag_cal(s);

#plot calibrated data
scatter3(t(:,1), t(:,2), t(:,3), [], 'b','filled')
legend('Raw Data', 'Calibrated')

#pretty title with standard error
title(cstrcat("Magnetic Ellipsoidal Fit Calibration\nStd. Err: ", num2str(std_err)));
