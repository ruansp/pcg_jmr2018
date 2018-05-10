function [H,Mu] = interpMultiPt(t0, H0, t, method)
% INTERPMULTIPT implements the interpolation with multiple points
% Inputs:
%   t0     ----- initial time scale
%   H0     ----- initial trajectory of frames
%   t      ----- interpolation time scale
%   method ----- interpolation methods: Euclidean space ("R"), Special 
%                 Euclidean group ("SE"), Pose Change Group ("PCG")
%
% Outputs:
%   H      ----- Interpolated framed trajectory
%   Mu     ----- Average frame
%
% Author: Sipu Ruan, ruansp@jhu.edu, September 2017

% Degree of the polynomial
n = size(H0,3);
deg = size(H0,3)-1;
m = 1;
for i = deg:-1:0
    T(m,:) = t.^i;
    m = m+1;
end

if strcmp(method, 'R')
    % Extract translation part
    for i = 1:n
        X0(:,i) = H0(1:3,4,i);
    end
    % Mean
    Mu = mean(X0,2);
    % Fit polynomial
    P = vPolyfit(t0, X0);
    X = P*T;
    
    for i = 1:size(X,2)
        H(:,:,i) = [eye(3) X(:,i); 0 0 0 1];
    end
elseif strcmp(method, 'SE') 
    % Transfer to se(3)
    for i = 1:n
        h0Hat = logm(H0(:,:,i));
        h0(1:3,i) = vex(h0Hat(1:3,1:3));
        h0(4:6,i) = h0Hat(1:3,4);
    end
    
    % Mean of frames
    mu = 1/n*sum(h0,2);
    muHat = [skew(mu(1:3)) mu(4:6); zeros(1,4)];
    Mu = expm(muHat);
    % Shift the original frames
    for i = 1:n
        HShift = Mu \ H0(:,:,i);
        hShiftHat = logm(HShift);
        hShift(1:3,i) = vex(hShiftHat(1:3,1:3));
        hShift(4:6,i) = hShiftHat(1:3,4);
    end
    
    % Fit a polynomial curve in se(3)
    P = vPolyfit(t0, hShift);   
    h = P * T;
    
    % Transfer back to SE(3)
    for i = 1:size(h,2)
        hHat = [skew(h(1:3,i)) h(4:6,i); zeros(1,4)];
        H(:,:,i) = Mu * expm(hHat);
    end
elseif strcmp(method, 'PCG')
    for i = 1:deg+1
        r0Hat = logm(H0(1:3,1:3,i));
        r0(:,i) = vex(r0Hat);
        x0(:,i) = H0(1:3,4,i);
    end
    
    % Mean of frames
    mur = 1/n*sum(r0,2);
    mux = 1/n*sum(x0,2);
    Mu = [expm(skew(mur)) mux; 0 0 0 1];
    % Shift the original frames
    for i = 1:n
        HShift = Mu \ H0(:,:,i);
        rShiftHat = logm(HShift(1:3,1:3));
        rShift(:,i) = vex(rShiftHat);
        xShift(:,i) = HShift(1:3,4);
    end
    
    % Fit a polynomial curve separately
    Pr = vPolyfit(t0, rShift);
    Px = vPolyfit(t0, xShift);
    
    r = Pr * T;
    x = Px * T;
    
    % Transfer back to PCG(3)
    for i = 1:size(r,2)
        rHat = skew(r(:,i));
        H(:,:,i) = Mu * [expm(rHat) x(:,i); 0 0 0 1];
    end
end
end

function P = vPolyfit(t, X)
% VPOLYFIT implements polynomial curve fitting in vector space
P = vander(t) \ X';
P = P';
end