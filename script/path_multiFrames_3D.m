close all; clear; clc;
% Illustration of multi-frame interpolation
% Author: Sipu Ruan, ruansp@jhu.edu, 2017

%% Define original frames
figure; hold on; axis equal; grid on;
drawFrame(eye(4),10);

% Define frames
axisF = [[0;0;1] [1;0;0] [1;2;0] [1;0;0] [3;2;6]];
thF = [pi/3 pi/6 0 -pi/4 pi/5];
dF = [[0;0;1] [0;0;1.1] [0;1;0] [0;1;1] [0.5;0.5;2]];

% Convert to Homo transformation matrix
Nframes = size(axisF,2);
for i = 1:Nframes
    R = rot_axis(axisF(:,i), thF(:,i));
    dF(:,i) = 30*dF(:,i);
    
    H(:,:,i) = [R dF(:,i); 0 0 0 1];
    drawFrame(H(:,:,i), 5);
end
plot3(dF(1,:),dF(2,:),dF(3,:),'k');
plot3(dF(1,:),dF(2,:),dF(3,:),'k*');

%% Interpolation
% Choice of parameterization of original time steps
t0 = timeParam(Nframes, 'even', [], []);
t0R = timeParam(Nframes, 'dist', H, 'R');
t0PCG = timeParam(Nframes, 'dist', H, 'PCG');
t0SE = timeParam(Nframes, 'dist', H, 'SE');

% Time steps for interpolation
sc = 10;
dt = 1/(sc*Nframes-1);
t = 0:dt:1;

% Interpolation
[HMultiR, muR] = interpMultiPt( t0R, H, t, 'R' );
[HMultiSE, muSE] = interpMultiPt( t0SE, H, t, 'SE' );
[HMultiPCG, muPCG] = interpMultiPt( t0PCG, H, t, 'PCG' );
HSE = interpX( t0, H, t, 'SE' );
HPCG = interpX( t0, H, t, 'PCG' );
for i = 1:size(t,2)
    xMultiR(:,i) = HMultiR(1:3,4,i);
    xMultiSE(:,i) = HMultiSE(1:3,4,i);
    xMultiPCG(:,i) = HMultiPCG(1:3,4,i);
    xSE(:,i) = HSE(1:3,4,i);
    xPCG(:,i) = HPCG(1:3,4,i);
end

plot3(xMultiR(1,:),xMultiR(2,:),xMultiR(3,:), 'k.');
plot3(xMultiSE(1,:),xMultiSE(2,:),xMultiSE(3,:), 'b--');
plot3(xMultiPCG(1,:),xMultiPCG(2,:),xMultiPCG(3,:), 'r');

%% Function for 3D rotation
% rotating along an axis defined by u, with angle theta
function R = rot_axis(u, theta)
u = u./norm(u,2);
S = [    0  u(3) -u(2);
      -u(3)   0   u(1);
      u(2) -u(1)   0  ];
R = eye(3) + sin(theta)*S + (1-cos(theta))*S^2;
end