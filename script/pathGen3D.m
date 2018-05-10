close all; clear; clc
% Path interpolation in 2D
% Author: Sipu Ruan, ruansp@jhu.edu, 2017

%% frame O
O = [0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 1 1];
plot3([O(1,1), O(1,2)], [O(2,1), O(2,2)], [O(3,1), O(3,2)],...
    'r-', 'LineWidth', 2)
hold on;
plot3([O(1,1), O(1,3)], [O(2,1), O(2,3)], [O(3,1), O(3,3)],...
    'g-', 'LineWidth', 2)
plot3([O(1,1), O(1,4)], [O(2,1), O(2,4)], [O(3,1), O(3,4)],...
    'b-', 'LineWidth', 2)

%% frame A, B
% Transform Ho_oa, Ha_ab, Ho_ob
% Rotation
sigma = 0.01;
Ro_oa = eulerzyz(2*pi/3,0,0);
Ro_ob = eulerzyz(-pi/3+sigma,-pi+sigma,pi/6);
% Translation
do_oa = 3*[-1; 2; 2]; do_ob = 3*[1; .5; 1];
% Homogeneous Matrix
Ho_oa = [Ro_oa do_oa; 0 0 0 1];
Ho_ob = [Ro_ob do_ob; 0 0 0 1];
Ha_ab = inv(Ho_oa) * Ho_ob;

% frames A and B
A = Ho_oa * O;
B = Ho_ob * O;

% plot the frames
plot3([A(1,1), A(1,2)], [A(2,1), A(2,2)], [A(3,1), A(3,2)],...
    'r-', 'LineWidth', 2)
plot3([A(1,1), A(1,3)], [A(2,1), A(2,3)], [A(3,1), A(3,3)],...
    'g-', 'LineWidth', 2)
plot3([A(1,1), A(1,4)], [A(2,1), A(2,4)], [A(3,1), A(3,4)],...
    'b-', 'LineWidth', 2)

plot3([B(1,1), B(1,2)], [B(2,1), B(2,2)], [B(3,1), B(3,2)],...
    'r-', 'LineWidth', 2)
plot3([B(1,1), B(1,3)], [B(2,1), B(2,3)], [B(3,1), B(3,3)],...
    'g-', 'LineWidth', 2)
plot3([B(1,1), B(1,4)], [B(2,1), B(2,4)], [B(3,1), B(3,4)],...
    'b-', 'LineWidth', 2)

%% Path Generation: SE(3)
% Ho_ab, Ro_ab, do_ab
Ho_ab = Ho_oa * Ha_ab * inv(Ho_oa);
Ro_ab = Ho_ab(1:3,1:3);
do_ab = Ho_ab(1:3,4);

% path generation: a helical path
t = 0:0.01:1;
for i=1:length(t)
    Hse(:,:,i) = expm(t(i)*logm(Ho_ab)) * Ho_oa;
    Ase(:,:,i) = Hse(:,:,i) * O;
    
    if mod(i,10) == 1
        plot3([Ase(1,1,i), Ase(1,2,i)], [Ase(2,1,i), Ase(2,2,i)],...
            [Ase(3,1,i), Ase(3,2,i)], 'r--', 'LineWidth', 1)
        plot3([Ase(1,1,i), Ase(1,3,i)], [Ase(2,1,i), Ase(2,3,i)],...
            [Ase(3,1,i), Ase(3,3,i)], 'g--', 'LineWidth', 1)
        plot3([Ase(1,1,i), Ase(1,4,i)], [Ase(2,1,i), Ase(2,4,i)],...
            [Ase(3,1,i), Ase(3,4,i)], 'b--', 'LineWidth', 1)
    end
    
    Ase_origin(:,i) = [Ase(1,1,i); Ase(2,1,i); Ase(3,1,i)];
end
plot3(Ase_origin(1,:), Ase_origin(2,:), Ase_origin(3,:),...
    'k--', 'LineWidth', 1.2)
%% Path Generation: PCG(3)
% to_ab
to_oa = Ho_oa(1:3,4); Ro_oa = Ho_oa(1:3,1:3);
to_ab = do_ab - (eye(3)-Ro_ab)*to_oa;

% path generation
for i=1:length(t)
    Rpcg(:,:,i) = expm(t(i)*logm(Ro_ab)) * Ro_oa;
    tpcg(:,:,i) = t(i)*to_ab + to_oa;
    Apcg(:,:,i) = Rpcg(:,:,i)*O(1:3,:) + tpcg(:,:,i);
    
    if mod(i,10) == 1
        plot3([Apcg(1,1,i), Apcg(1,2,i)], [Apcg(2,1,i), Apcg(2,2,i)],...
            [Apcg(3,1,i), Apcg(3,2,i)], 'r-', 'LineWidth', 1)
        plot3([Apcg(1,1,i), Apcg(1,3,i)], [Apcg(2,1,i), Apcg(2,3,i)],...
            [Apcg(3,1,i), Apcg(3,3,i)], 'g-', 'LineWidth', 1)
        plot3([Apcg(1,1,i), Apcg(1,4,i)], [Apcg(2,1,i), Apcg(2,4,i)],...
            [Apcg(3,1,i), Apcg(3,4,i)], 'b-', 'LineWidth', 1)
    end
    
    Apcg_origin(:,i) = [Apcg(1,1,i); Apcg(2,1,i); Apcg(3,1,i)];
end
% plot the traj of origin of the frames
plot3(Apcg_origin(1,:), Apcg_origin(2,:), Apcg_origin(3,:),...
    'k-', 'LineWidth', 1.2)
axis equal;
grid on;

%% Function for 3D rotation
% rotating along an axis defined by u, with angle theta
function R = rot_axis(u, theta)
u = u./norm(u,2);
S = [    0  u(3) -u(2);
      -u(3)   0   u(1);
      u(2) -u(1)   0  ];
R = eye(3) + sin(theta)*S + (1-cos(theta))*S^2;
end

% zyz Euler Angles
function R = eulerzyz(alpha, beta, gamma)
Rz1 = [cos(alpha) -sin(alpha) 0;
       sin(alpha)  cos(alpha) 0;
                0           0 1];
Ry2 = [ cos(beta) 0 sin(beta);
                0 1         0;
       -sin(beta) 0 cos(beta)];
Rz3 = [cos(gamma) -sin(gamma) 0;
       sin(gamma)  cos(gamma) 0;
                0           0 1];
R = Rz1*Ry2*Rz3;
end