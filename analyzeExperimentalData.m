% This code is part of the Master's degree thesis "Estimation of oxygen
% consumption from partial pressure gradients in cortex", written by 
% Marte Julie Sætra, May 2016.
% *****************************
% An example on how to extract the oxygen consumption rate
% from measurements of oxygen partial pressue pO2, by using the Laplace
% method.
% The method is described in Chapter 2. Smoothing is described in 
% Chapter 3. Further detailt are covered by Chapter 7.

clear all
close all

% choose data set
dataset = 5;
load(['/home/martejulie/master_project_data/dataset', num2str(dataset),'/po2_data.mat']);
load(['/home/martejulie/master_project_data/dataset', num2str(dataset),'/dist_2_vessel.mat']);

% create P matrix
[P_original, Hx, Hy, d] = extractPO2matrix(po2_data);
xx = {Hy, Hx};
% create R matrix
r = extractRmatrix(Hx, Hy, d, vessel_coords);

% Estimate M as a function of p_cutoff, for different smoothing factors 
p_cutoff = 10:0.1:70;
p_smooth = [0.9, 0.5, 0.1, 0.01, 0.005];
M_vectors = {};
for j = 1:length(p_smooth)
    % smooth P
    P_smoothed = csaps(xx, P_original, p_smooth(j), xx);
    % estimate M
    del2P = 4*del2(P_smoothed, d);
    M_vec = [];
    for i = 1:length(p_cutoff)
        M_matrix = del2P(P_smoothed <= p_cutoff(i));
        M_vec(i) = mean(M_matrix);        
    end
    M_vectors{j} = M_vec;
end

% Estimate M for different regions parted with respect to r
p = 0.1;
P_smoothed = csaps(xx, P_original, p, xx);
del2P = 4*del2(P_smoothed, d);
r_low = 10:10:40;
r_high = 30:10:60;
MofR = estimateMforDifferentRegions(r, r_low, r_high, del2P);

% Estimate M for different regions parted with respect to pO2
p_low = 0:10:10;
p_high = 10:10:40;
MofPO2 = estimateMforDifferentRegions(P_smoothed, p_low, p_high, del2P);

% ******************************
% plot
P_original = flipud(P_original);
P_smoothed = flipud(P_smoothed);
del2P = flipud(del2P);

map = makeColorMap([1,1,1], [1,0,0], 100);

figure(1);
imagesc(P_original, [0, max(P_original(:))]);
colormap(map);
axis xy;
colorbar;
title(['\textbf{$\mathrm{pO_2}$ for dataset $', num2str(dataset),'$}'], 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16)

figure(2); 
fim = imagesc(P_smoothed, [0, max(P_smoothed(:))]);
colormap(map);
colorbar;
axis xy;
title(['\textbf{$\mathrm{pO_2}$ for dataset $', num2str(dataset), ', \,p = ', num2str(p), '$}'], 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16)

figure(3); 
mim = imagesc(del2P); 
colormap(NegativeEnhancingColormap(100, [min(del2P(:)) max(del2P(:))], [0 0 1], [1 0 0], 1));
colorbar;
axis xy;
title(['\textbf{$\nabla^2 \mathrm{pO_2} = M, \,p = ', num2str(p), '$}'], 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$y$', 'Interpreter', 'latex');
set(gca, 'fontsize', 16)

figure(4); 
plot(P_smoothed(:), del2P(:),'.b');
title(['\textbf{$M$ vs. $\mathrm{pO_2},\, p = ', num2str(p), '$}'], 'Interpreter', 'latex');
xlabel('$\mathrm{pO_2}$ [mmHg]', 'Interpreter', 'latex');
ylabel('Local $M$ [mmHg/$\mu m$]', 'Interpreter', 'latex');
set(gca, 'fontsize', 16)

figure(5)
hold on
legendInfo = {};
colors = {'b', 'r', 'k', [0,153./255,76./255], [204./255, 204./255 0]};
for i = 1:length(M_vectors)
    plot(p_cutoff, M_vectors{i}, 'Color', colors{i})
    legendInfo{i} = ['p = ', num2str(p_smooth(i))];
end
legend(legendInfo)
grid minor
title(['\textbf{$M\mathrm{_{est}}$ vs. $p\mathrm{_{cutoff}}$ for dataset $', num2str(dataset),'$}'], 'Interpreter', 'latex');
xlabel('$\mathrm{pO_2}$ [mmHg]', 'Interpreter', 'latex');
ylabel('$M$ [mmHg/$\mu m$]', 'Interpreter', 'latex');
set(gca, 'fontsize', 16)
%axis([10 60 0 0.1])
hold off

figure(6);
fig = imagesc(r_low, r_high, MofR, [0.0, 0.05]);
set(gca, 'fontsize', 16, 'YTick', [30:10:60])
colormap(map);
colorbar;
axis xy;
title(['\textbf{$M, \, p = ', num2str(p), '$}'], 'Interpreter', 'latex');
xlabel('$r\mathrm{_{low}}$', 'Interpreter', 'latex');
ylabel('$r\mathrm{_{high}}$', 'Interpreter', 'latex');

figure(7);
fig = imagesc(p_low, p_high, MofPO2, [0.0 0.05]);
set(gca, 'fontsize', 16, 'XTick', [0,10], 'YTick', [10,20,30,40]);
colormap(map);
colorbar;
axis xy;
title(['\textbf{$M, \, p = ', num2str(p), '$}'], 'Interpreter', 'latex');
xlabel('$p\mathrm{_{low}}$', 'Interpreter', 'latex');
ylabel('$p\mathrm{_{high}}$', 'Interpreter', 'latex');