function [x,dx,y,dy] = grid2d(a,b,nx,c,d,ny)
%% Description:

% INPUTS
% [a,b]: x range, nx: the desided number of cells.
% [c,d]: y range, ny: the desided number of Cells.

% OUTPUTS
% x: vector x's coordinates, dx: cell size.
% y: vector y's coordinates, dy: cell size.

%% Compute cell size
dx = (b-a)/nx; dy = (d-c)/ny;

%% Compute vector arrays
[x,y] = meshgrid(a+dx/2:dx:b,c+dy/2:dy:d);