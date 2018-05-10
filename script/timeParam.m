function t = timeParam(N, method, H, m)
% TIMEPARAM implements the initial time parameterization for a sequence
% Inputs:
%    N      ----- Number of time steps
%    method ----- Method of parameterization: evenly distributed ("even"),
%                   distance-based ("dist")
%    H      ----- A sequence of frames
%    m      ----- Method for computing the distance metric between two
%                   frames: Euclidean distance ("R"), SE(3) metric ("SE"),
%                   PCG(3) metric ("PCG")
%
% Output:
%    t      ----- Parameterized time scale
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2018

if strcmp(method, 'even')
    dt = 1/(N-1);
    t = 0:dt:1;
elseif strcmp(method, 'dist')
    t = zeros(1,N);
    for i = 1:N-1
        d(i) = metric(H(:,:,i), H(:,:,i+1), m);
    end
    D = sum(d);
    for i = 2:N
        t(i) = t(i-1)+d(i-1)/D;
    end
end

function d = metric(H1, H2, method)
if strcmp(method, 'R')
    X1 = H1(1:3,4);
    X2 = H2(1:3,4);
    d = norm(X1-X2);
elseif strcmp(method, 'SE')
    d = norm(logm(H1\H2), 'fro');
elseif strcmp(method, 'PCG')
    R1 = H1(1:3,1:3); R2 = H2(1:3,1:3);
    x1 = H1(1:3,4); x2 = H2(1:3,4);
    
    d = norm(x1-x2)+norm(logm(R1\R2),'fro');
    d = norm(x1-x2);
end
