function [ r ] = generateRmatrix( Hx, Hy, d, vessel_coords )
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
r = sqrt( (X - floor(vessel_coords(1)/d)*d).^2 + (Y - floor(vessel_coords(2)/d)*d).^2 );
end

