function polyEllipsoidQuadrotor(H)
% Generate and plot a quadrotor constructed by a union of ellipsoids
% Input:
%   H ----- Homogeneous transformation matrix for the base of quadrotor
%
% Author: Sipu Ruan, ruansp@jhu.edu, 2017

% center of each part
base = [0;0;0];
wing = base+[1;0;0];
wing2 = wing+[1;0;0.25];
propell = wing2+[0;0;0.25];

% base
[Xb, Yb, Zb] = ellipsoid(base(1),base(2),base(3),1,1,0.3);
[Xb, Yb, Zb] = ellipsoidXform(Xb,Yb,Zb,H);
surf(Xb,Yb,Zb, 'FaceColor', 'r', 'EdgeColor', 'none');

% wings and propellers
[Xw, Yw, Zw] = ellipsoid(wing(1),wing(2),wing(3),1,0.2,0.2);
[Xw2,Yw2,Zw2] = ellipsoid(wing2(1),wing2(2),wing2(3),0.1,0.1,0.25);
[Xp,Yp,Zp] = ellipsoid(propell(1),propell(2),propell(3),0.75,0.75,0.02);

xWing(1,:) = reshape(Xw, 1, size(Xw,1)*size(Xw,1));
xWing(2,:) = reshape(Yw, 1, size(Xw,1)*size(Xw,1));
xWing(3,:) = reshape(Zw, 1, size(Xw,1)*size(Xw,1));

xWing2(1,:) = reshape(Xw2, 1, size(Xw2,1)*size(Xw2,1));
xWing2(2,:) = reshape(Yw2, 1, size(Xw2,1)*size(Xw2,1));
xWing2(3,:) = reshape(Zw2, 1, size(Xw2,1)*size(Xw2,1));

xPro(1,:) = reshape(Xp, 1, size(Xp,1)*size(Xp,1));
xPro(2,:) = reshape(Yp, 1, size(Xp,1)*size(Xp,1));
xPro(3,:) = reshape(Zp, 1, size(Xp,1)*size(Xp,1));

[Xw, Yw, Zw] = ellipsoidXform(Xw,Yw,Zw,H);
[Xw2, Yw2, Zw2] = ellipsoidXform(Xw2,Yw2,Zw2,H);
[Xp, Yp, Zp] = ellipsoidXform(Xp,Yp,Zp,H);
surf(Xw,Yw,Zw,'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.9); 
surf(Xw2,Yw2,Zw2,'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.9); 
surf(Xp,Yp,Zp,'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.9);

R = H(1:3,1:3);
t = H(1:3,4);
for i = 2:4
    xWing(:,:,i) = R*rot_axis([0,0,1], pi/2*(i-1))*xWing(:,:,1)+t;
    xWing2(:,:,i) = R*rot_axis([0,0,1], pi/2*(i-1))*xWing2(:,:,1)+t;
    xPro(:,:,i) = R*rot_axis([0,0,1], pi/2*(i-1))*xPro(:,:,1)+t;
    
    Xw(:,:,i) = reshape(xWing(1,:,i), size(Xw,1), size(Xw,1));
    Yw(:,:,i) = reshape(xWing(2,:,i), size(Xw,1), size(Xw,1));
    Zw(:,:,i) = reshape(xWing(3,:,i), size(Xw,1), size(Xw,1));
    
    Xw2(:,:,i) = reshape(xWing2(1,:,i), size(Xw2,1), size(Xw2,1));
    Yw2(:,:,i) = reshape(xWing2(2,:,i), size(Xw2,1), size(Xw2,1));
    Zw2(:,:,i) = reshape(xWing2(3,:,i), size(Xw2,1), size(Xw2,1));
    
    Xp(:,:,i) = reshape(xPro(1,:,i), size(Xp,1), size(Xp,1));
    Yp(:,:,i) = reshape(xPro(2,:,i), size(Xp,1), size(Xp,1));
    Zp(:,:,i) = reshape(xPro(3,:,i), size(Xp,1), size(Xp,1));
    
    surf(Xw(:,:,i),Yw(:,:,i),Zw(:,:,i),'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    surf(Xw2(:,:,i),Yw2(:,:,i),Zw2(:,:,i),'FaceColor', 'b', 'EdgeColor', 'none', 'FaceAlpha', 0.9);
    surf(Xp(:,:,i),Yp(:,:,i),Zp(:,:,i),'FaceColor', 'g', 'EdgeColor', 'none', 'FaceAlpha', 0.9);
end
end

function [XX,YY,ZZ] = ellipsoidXform(X,Y,Z,H)
x(1,:) = reshape(X, 1, size(X,1)*size(X,1));
x(2,:) = reshape(Y, 1, size(Y,1)*size(Y,1));
x(3,:) = reshape(Z, 1, size(Z,1)*size(Z,1));

x2 = H*[x;ones(1,size(x,2))];

XX = reshape(x2(1,:), size(X,1), size(X,1));
YY = reshape(x2(2,:), size(Y,1), size(Y,1));
ZZ = reshape(x2(3,:), size(Z,1), size(Z,1));
end
