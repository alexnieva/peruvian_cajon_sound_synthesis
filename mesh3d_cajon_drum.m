% mesh3d_cajon_drum.m
%
% A Matlab script implementing a rectilinear waveguide mesh.  This
% scheme is based on the STK Mesh2D class.  The calculations make use
% of two sets of wave variable matrices that are alternated to save
% state.  Note that the mesh boundaries do not appear to be perfectly
% fixed in the plots.  This occurs because the plotted values do not
% include the outer boundary junctions.
%
% by Gary Scavone
% MUMT614, McGill University, 2004.
% Modified by Alex Nieva from mesh2d.m
%    
%  z_ _ _ _ _ _ x
%  |
%  |
%  |
%  |
%  |
%  |
%  y

% Length of signal to calculate.
N = 5000;

% Sample rate for playback.
fs = 22050;

% Number of (x,y,z) grid junctions (assuming rectangular shape).
% Resolution deltaS = sqrt(3)*c/fs, number of bins L/deltaS

NH = 17; % y axis
NW = 11; % x axis
NP = 11; % z axis

% Initialize calculation matrices.
vxp = zeros(NH,NW,NP);
vxm = zeros(NH,NW,NP);
vyp = zeros(NH,NW,NP);
vym = zeros(NH,NW,NP);
vzp = zeros(NH,NW,NP);
vzm = zeros(NH,NW,NP);

vxp1 = zeros(NH,NW,NP);
vxm1 = zeros(NH,NW,NP);
vyp1 = zeros(NH,NW,NP);
vym1 = zeros(NH,NW,NP);
vzp1 = zeros(NH,NW,NP);
vzm1 = zeros(NH,NW,NP);

v = zeros(NH-1,NW-1,NP-1);
y2 = zeros(1,N);

% Initialization of velocity on the front plate.
valm = zeros(NH-1,NW-1);
valm(7:11,4:7) = 0.9;

vxp1(1:NH-1,1:NW-1,1) = valm;
vyp1(1:NH-1,1:NW-1,1) = valm;
vzp1(1:NH-1,1:NW-1,1) = valm;
vxm1(1:NH-1,2:NW,1) = valm;
vym1(2:NH,1:NW-1,1) = valm;
vzm1(1:NH-1,1:NW-1,2) = valm;

% Initial velocity
v = (1/3) * (vxp1(1:NH-1,1:NW-1,1:NP-1) + vxm1(1:NH-1,2:NW,1:NP-1) + ...
           vyp1(1:NH-1,1:NW-1,1:NP-1) + vym1(2:NH,1:NW-1,1:NP-1) + ...
           vzp1(1:NH-1,1:NW-1,1:NP-1) + vzm1(1:NH-1,1:NW-1,2:NP));

for n = 1:N,

  if ( mod(n,2) == 0 ), % tick0

v = (1/3)* (vxp(1:NH-1,1:NW-1,1:NP-1) + vxm(1:NH-1,2:NW,1:NP-1) + ...
           vyp(1:NH-1,1:NW-1,1:NP-1) + vym(2:NH,1:NW-1,1:NP-1) + ...
           vzp(1:NH-1,1:NW-1,1:NP-1) + vzm(1:NH-1,1:NW-1,2:NP));

    % Update outgoing junction wave components.
    vxp1(1:NH-1,2:NW,1:NP-1) = v - vxm(1:NH-1,2:NW,1:NP-1);
    vyp1(2:NH,1:NW-1,1:NP-1) = v - vym(2:NH,1:NW-1,1:NP-1);
    vzp1(1:NH-1,1:NW-1,2:NP) = v - vzm(1:NH-1,1:NW-1,2:NP);
    vxm1(1:NH-1,1:NW-1,1:NP-1) = v - vxp(1:NH-1,1:NW-1,1:NP-1);
    vym1(1:NH-1,1:NW-1,1:NP-1) = v - vyp(1:NH-1,1:NW-1,1:NP-1);
    vzm1(1:NH-1,1:NW-1,1:NP-1) = v - vzp(1:NH-1,1:NW-1,1:NP-1);

    % Do boundary reflections.
    vxp1(1:NH-1,1,1:NP-1) = -0.9999*vxm(1:NH-1,1,1:NP-1); % Left
    vxm1(1:NH-1,NW,1:NP-1) = -0.9999*vxp(1:NH-1,NW,1:NP-1); % Right
    vyp1(1,1:NW-1,1:NP-1) = -0.9999*vym(1,1:NW-1,1:NP-1); % Top
    vym1(NH,1:NW-1,1:NP-1) = -0.9999*vyp(NH,1:NW-1,1:NP-1); % Bottom
    vzp1(1:NH-1,1:NW-1,1) = -0.1*vzp(1:NH-1,1:NW-1,1); % Front
    vzm1(1:NH-1,1:NW-1,NP) = -0.9999*vzm(1:NH-1,1:NW-1,NP); % Back
  
  else % tick1

    v = (1/3) * (vxp1(1:NH-1,1:NW-1,1:NP-1) + vxm1(1:NH-1,2:NW,1:NP-1) + ...
               vyp1(1:NH-1,1:NW-1,1:NP-1) + vym1(2:NH,1:NW-1,1:NP-1) + ...
               vzp1(1:NH-1,1:NW-1,1:NP-1) + vzm1(1:NH-1,1:NW-1,2:NP));

    % Update outgoing junction wave components.
    vxp(1:NH-1,2:NW,1:NP-1) = v - vxm1(1:NH-1,2:NW,1:NP-1);
    vyp(2:NH,1:NW-1,1:NP-1) = v - vym1(2:NH,1:NW-1,1:NP-1);
    vzp(1:NH-1,1:NW-1,2:NP) = v - vzm1(1:NH-1,1:NW-1,2:NP);
    vxm(1:NH-1,1:NW-1,1:NP-1) = v - vxp1(1:NH-1,1:NW-1,1:NP-1);
    vym(1:NH-1,1:NW-1,1:NP-1) = v - vyp1(1:NH-1,1:NW-1,1:NP-1);
    vzm(1:NH-1,1:NW-1,1:NP-1) = v - vzp1(1:NH-1,1:NW-1,1:NP-1);
    
    % Do boundary reflections.
    vxp(1:NH-1,1,1:NP-1) = -0.9999*vxm1(1:NH-1,1,1:NP-1); % Left
    vxm(1:NH-1,NW,1:NP-1) = -0.9999*vxp1(1:NH-1,NW,1:NP-1); % Right
    vyp(1,1:NW-1,1:NP-1) = -0.9999*vym1(1,1:NW-1,1:NP-1); % Top
    vym(NH,1:NW-1,1:NP-1) = -0.9999*vyp1(NH,1:NW-1,1:NP-1); % Bottom
    vzp(1:NH-1,1:NW-1,1) = -0.1*vzp1(1:NH-1,1:NW-1,1); % Front
    vzm(1:NH-1,1:NW-1,NP) = -0.9999*vzm1(1:NH-1,1:NW-1,NP); % Back
  end

  % Do output calculation.
  y2(n) = v(NH-8,NW-5,NP-1);
%keyboard;
end

plot(y2);
sound(0.99*y2/max(abs(y2)), fs);
