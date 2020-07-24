close all; clear; clc;
% Animation of multi-frame interpolation
% Author: Sipu Ruan, ruansp@jhu.edu, 2017

%% Define original frames
vidObj = VideoWriter('uav_multi_frames.avi');
open(vidObj);

figure; hold on; axis equal; axis off;
view([1,1,0.5])

drawFrame(eye(4),10);

% Define frames
axisF = [[0;0;1] [1;0;0] [1;2;0] [1;0;0] [3;2;6]];
thF = [pi/3 pi/6 0 -pi/4 pi/5];
dF = [[0;0;1] [0;0;1.1] [0;1;0] [0;1;0.8] [0.5;0.5;1]];

% Convert to Homo transformation matrix
Nframes = size(axisF,2);
for i = 1:Nframes
    R = rot_axis(axisF(:,i), thF(:,i));
    dF(:,i) = 30*dF(:,i);
    
    H(:,:,i) = [R dF(:,i); 0 0 0 1];
    drawFrame(H(:,:,i), 5);
end
plot3(dF(1,:),dF(2,:),dF(3,:),'k*');

%% Interpolation
% Choice of parameterization of original time steps
t0 = timeParam(Nframes, 'dist', H, 'PCG');

% Time steps for interpolation
sc = 20;
dt = 1/(sc*Nframes-1);
t = 0:dt:1;

% Interpolation
[HMulti, mu] = interpMultiPt( t0, H, t, 'PCG' );
for i = 1:size(t,2)
    xMulti(:,i) = HMulti(1:3,4,i);
end

plot3(xMulti(1,:),xMulti(2,:),xMulti(3,:), 'b--');

%% Animation
for i = 1:size(t,2)   
    plot3(dF(1,:),dF(2,:),dF(3,:),'k*');
    
    hold on; axis off; axis equal;
    view([1,0,0.2])
    
    plot3(xMulti(1,1:i),xMulti(2,1:i),xMulti(3,1:i), 'b-');
    
    polyEllipsoidQuadrotor(HMulti(:,:,i))
    
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
    
    hold off;
end

close(vidObj);