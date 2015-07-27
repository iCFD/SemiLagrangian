function [dF] = Upwindflux(a,uo,ind)
% Vectorized version of Two sided upwinding flux
%
% coded by Manuel Diaz, NTU, 2012.12.20
%
% Domain cells (I{i}) reference:
%
%                |           |   u(i)    |           |
%                |  u(i-1)   |___________|           |
%                |___________|           |   u(i+1)  |
%                |           |           |___________|
%             ...|-----0-----|-----0-----|-----0-----|...
%                |    i-1    |     i     |    i+1    |
%                |-         +|-         +|-         +|
%              i-3/2       i-1/2       i+1/2       i+3/2
%
% Vector variables, (periodic BC behavior included)
um =circshift(uo,ind);  %u(i-1)
up =circshift(uo,-ind); %u(i+1)

% for flux splitting
a_p = max(a,0); % a{+}
a_m = min(a,0); % a{-}

% Double sided upwind fluxes
fl = a_p.*um + a_m.*uo; % left
fr = a_p.*uo + a_m.*up; % right

% df values
dF = fr-fl;