function [ r ] = generateRmatrix( Hx, Hy, d, vessel_coords, pixels_per_100u )
% Generates a matrix with values of the distance to vessel
% Arg:
%   Hx (vector): x vector
%   Hy (vector): y vector
%   d (double): spatial spacing
%   vessel_coords (vector): vessel possition
%
% Returns:
%   r (matrix): distance to vessel
%
[X, Y] = meshgrid(Hx, Hy);
r = sqrt( (X - vessel_coords(1)*100/pixels_per_100u).^2 + (Y - vessel_coords(2)*100/pixels_per_100u).^2 );
end

