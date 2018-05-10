close all; clear; clc
% Path interpolation in 2D
% Author: Sipu Ruan, ruansp@jhu.edu, 2017

%% frame O
O = [0 1 0; 0 0 1; 1 1 1];
plot(O(1,1:2), O(2,1:2), 'r-', 'LineWidth', 2)
hold on;
plot([O(1,1), O(1,3)], [O(2,1), O(2,3)], 'b-', 'LineWidth', 2)

%% frame A, B
% Transform Ho_oa, Ha_ab, Ho_ob
thA = 2*pi/3; do_oa = [-5; 4];
thB = -pi/3; do_ob = [4; 2];
Ro_oa = rot2d(thA);
Ro_ob = rot2d(thB);

Ho_oa = [Ro_oa do_oa; 0 0 1];
Ho_ob = [Ro_ob do_ob; 0 0 1];
Ha_ab = inv(Ho_oa) * Ho_ob;

% frames A and B
A = Ho_oa * O;
B = Ho_ob * O;

% plot the frames
plot(A(1,1:2), A(2,1:2), 'r-', 'LineWidth', 2)
plot([A(1,1), A(1,3)], [A(2,1), A(2,3)], 'b-', 'LineWidth', 2)
plot(B(1,1:2), B(2,1:2), 'r-', 'LineWidth', 2)
plot([B(1,1), B(1,3)], [B(2,1), B(2,3)], 'b-', 'LineWidth', 2)

%% Path Generation: SE(2)
% Ho_ab, Ro_ab, do_ab
Ho_ab = Ho_oa * Ha_ab * inv(Ho_oa);
Ro_ab = Ho_ab(1:2,1:2);
do_ab = Ho_ab(1:2,3);

% path generation: a helical path
t = 0:0.01:1;
for i=1:length(t)
    % transformation of SE(2)
    Hse(:,:,i) = expm(t(i)*logm(Ho_ab)) * Ho_oa;
    % moving frame
    Ase(:,:,i) = Hse(:,:,i) * O;
    
    % plot the moving frame
    if mod(i,10) == 1
        plot(Ase(1,1:2,i), Ase(2,1:2,i), 'r--', 'LineWidth', 1)
        plot([Ase(1,1,i), Ase(1,3,i)], [Ase(2,1,i), Ase(2,3,i)],...
            'b--', 'LineWidth', 1)
    end
    
    % track the origin of the moving frame
    Ase_origin(:,i) = [Ase(1,1,i); Ase(2,1,i)];
end
plot(Ase_origin(1,:), Ase_origin(2,:), 'k--', 'LineWidth', 1.2)

%% Path Generation: PCG(2)
% to_ab
to_oa = Ho_oa(1:2,3); Ro_oa = Ho_oa(1:2,1:2);
to_ab = do_ab - (eye(2)-Ro_ab)*to_oa;

% path generation
for i=1:length(t)
    Rpcg(:,:,i) = expm(t(i)*logm(Ro_ab)) * Ro_oa;
    tpcg(:,:,i) = t(i)*to_ab + to_oa;
    Apcg(:,:,i) = Rpcg(:,:,i)*O(1:2,:) + tpcg(:,:,i);
    
    if mod(i,10) == 1
    plot(Apcg(1,1:2,i), Apcg(2,1:2,i), 'r-', 'LineWidth', 1)
    plot([Apcg(1,1,i), Apcg(1,3,i)], [Apcg(2,1,i), Apcg(2,3,i)],...
        'b-', 'LineWidth', 1)
    end
    
    Apcg_origin(:,i) = [Apcg(1,1,i); Apcg(2,1,i)];
end
plot(Apcg_origin(1,:), Apcg_origin(2,:), 'k', 'LineWidth', 1.2)

axis equal;
grid on;

%% Function for 2D rotation
function R = rot2d(th)
R = [cos(th) -sin(th); sin(th) cos(th)];
end