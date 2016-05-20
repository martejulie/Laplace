function [ P, Hx, Hy, d ] = generatePmatrix( po2_data, pixels_per_100u )
% Generates a pO2 grid from experimental data
% Arg:
%   po2_data (matrix): x values, y values, pO2 values
%
% Returns:
%   P (matrix): pO2 grid
%   Hx (vector): x vector
%   Hy (vector): y vector
%   d (double): spatial spacing

% extract data
xvec = po2_data(:,1)*100/pixels_per_100u; 
yvec = po2_data(:,2)*100/pixels_per_100u; 
pO2vec = po2_data(:,3);

% extract spatial grid
[hx dummy xi] = unique(xvec); 
[hy dummy yi] = unique(yvec); 
% calculate spatial spacing
dx = (max(hx(:))- min(hx(:))) / (length(hx)-1);
dy = (max(hy(:))- min(hy(:))) / (length(hy)-1);
d = mean([dx, dy]);
% create grid with uniform spatial spacing
Hx = min(hx(:)):d:ceil(max(hx(:))/d)*d;
Hy = min(hy(:)):d:ceil(max(hy(:))/d)*d;
% create pO2 grid
nx = length(Hx);
ny = length(Hy);
P = NaN(ny,nx);
for i=1:length(pO2vec)
    P(yi(i),xi(i)) = pO2vec(i);
end
end

