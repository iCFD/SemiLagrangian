function u0 = IC2d(x,y,option)
%% Initial Condition (IC)
% This subroutine creates an special IC for testing purposes.
% The present IC models a rectangular domain with four equally sized
% regions with diferrent initial state value, u, an adimensional value.
%
% by Manuel Diaz, manuel.ade'at'gmail.com 
% Institute of Applied Mechanics, 2012.09.06

%% Select your IC

switch option
    case {1} % 4 quadrants
        %  Region 1 = 1.0 [-]
        %  Region 2 = 0.0 [-]
        %  Region 3 = 0.5 [-]
        %  Region 4 = 1.5 [-]

        % Parameters
        nx = size(x,2);
        ny = size(y,1);

        % Preallocate u0,
        u0 = zeros(ny,nx);

        % Parameters of regions dimensions
        x_middle = ceil(nx/2);
        y_middle = ceil(ny/2);
        l_1 = 1:x_middle; l_2 = x_middle+1:nx;
        h_1 = 1:y_middle; h_2 = y_middle+1:ny;

        % Initial Condition for our 2D domain
        u0(h_1,l_1) = 0.70; % region 1
        u0(h_1,l_2) = 0.10; % region 2
        u0(h_2,l_1) = 0.90; % region 3
        u0(h_2,l_2) = 0.50; % region 4
        
    case {2} % Square Jump
        % parameters
        nx = size(x,2);
        ny = size(y,1);
        
        % Preallocate u0
        u0 = zeros(ny,nx);
        
        % Parameters of region
        x_1 = ceil(2*nx/5);
        x_2 = ceil(3*nx/5);
        y_1 = ceil(3*ny/5);
        y_2 = ceil(4*ny/5);
        l_1 = x_1:x_2;
        h_1 = y_1:y_2;
        
        u0(h_1,l_1) = 2; %Jump
        
    case {3} % non-centered Gaussian
        % assuming a domain [0,1]x[0,1]
        x0=0.5; y0=0.7; amp=2.; sigma=0.1;
        u0=amp*exp(-((x-x0).^2 + (y-y0).^2)/sigma^2);
end
