clc;
clear;
close all;

%% Adding paths
addpath('Exercise1');
addpath('Exercise2');
addpath('Exercise3');
%% Generation of Helix
tic;
torous_radius = 5;
ncurls = 17;
helix = generate_helix(torous_radius,ncurls,0);
draw_fig(helix,500);
hold on;
%% Generate Circle
circle_radius = 5;
circle = generate_circle(circle_radius,0);
draw_fig(circle,500);
toc;
title(['Radius of Circle = ',num2str(circle_radius) ' ; Radius of torous = ',...
    num2str(torous_radius) ' ; ncurl = ',num2str(ncurls)]);
%% Integration
val_integral = integrate_hw01(circle,helix);
fprintf('Value of the integral is %f\n',val_integral);

