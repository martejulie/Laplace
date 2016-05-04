% This code is part of the Master's degree thesis "Estimation of oxygen
% consumption from partial pressure gradients in cortex", written by 
% Marte Julie Sætra, May 2016.
% *****************************
% An example on how to generate example data, explained in
% detailt in Chapter 6.

clear all
close all

% set parameters
params = [80, 0.0065]; % po2 at the vessel wall, oxygen consumption M
rves = 6.0; % vessel radius
rt = 100.0; % tissue radius

% load dataset 5
load('/home/martejulie/master_project_data/dataset5/po2_data.mat');
load('/home/martejulie/master_project_data/dataset5/dataStruct_delay15.mat');
load('/home/martejulie/master_project_data/dataset5/dist_2_vessel.mat');

% extract data
xvec = po2_data(:,1); 
yvec = po2_data(:,2);
po2vec = po2_data(:,3);
evec = dataStruct.pO2_std_err;

% extract spatial grid from dataset 5
[hx dummy xi] = unique(xvec);       
[hy dummy yi] = unique(yvec);
% calculate spatial spacing
dx = (max(hx(:))- min(hx(:))) / (length(hx)-1);
dy = (max(hy(:))- min(hy(:))) / (length(hy)-1);
d = mean([dx, dy]);
% create grid with uniform spatial spacing
Hx = min(hx(:)):d:max(hx(:));       
Hy = min(hy(:)):d:max(hy(:));
[X, Y] = meshgrid(Hx, Hy);
% calculate distance to vessel
r = sqrt( (X - vessel_coords(1)).^2 + (Y - vessel_coords(2)).^2 );
r(r < rves) = rves;
% calculate pO2 values from the Standard Krogh Model
P_ideal = params(1) + 0.25 * params(2) * (r.^2 - rves) - 0.5 * params(2) * rt.^2 * log(r ./ rves);

% add noise (exponential model)
exptbl = table(po2vec, log10(evec));
expmdl = fitlm(exptbl);
mystd = 10.^(expmdl.Coefficients.Estimate(1)+expmdl.Coefficients.Estimate(2).*P_ideal);
P_noisy = normrnd(P_ideal, mystd);

% smooth with CSAPS
xx = {Hy, Hx};
p = 0.1;
P_smoothed = csaps(xx, P_noisy, p, xx);

% calculate M
del2P = 4*del2(P_smoothed, d);


% ******************************
% plot
P_ideal = flipud(P_ideal);
P_noisy = flipud(P_noisy);
P_smoothed = flipud(P_smoothed);
del2P = flipud(del2P);

map = makeColorMap([1,1,1], [1,0,0], 100);

figure(1);
imagesc(P_noisy, [0, max(P_noisy(:))]);
colormap(map);
axis xy;
colorbar;
title('\textbf{$\mathrm{pO_2}$ for Noisy Example Data}', 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$ ', 'Interpreter', 'latex');
set(gca,'fontsize',16)

figure(2); 
fim = imagesc(P_smoothed, [0, max(P_smoothed(:))]);
colormap(map);
colorbar;
axis xy;
title(['\textbf{Smoothed $\mathrm{pO_2},\, p = ', num2str(p), '$}'], 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
set(gca,'fontsize',16)

figure(3); 
mim = imagesc(del2P); 
colormap(NegativeEnhancingColormap(100, [min(del2P(:)) max(del2P(:))], [0 0 1], [1 0 0], 1));
colorbar;
axis xy;
title(['\textbf{$\nabla^2 \mathrm{pO_2} = M,\, p = ', num2str(p), '$}'], 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
set(gca,'fontsize',16)

figure(4); 
plot(P_smoothed(:), del2P(:),'.b');
title(['\textbf{$M$ vs. $\mathrm{pO_2},\, p = ', num2str(p), '$}'], 'Interpreter', 'latex');
xlabel('$\mathrm{pO_2}$ [mmHg]', 'Interpreter', 'latex');
ylabel('$M$ [mmHg/$\mu m$]', 'Interpreter', 'latex');
set(gca,'fontsize',16)

