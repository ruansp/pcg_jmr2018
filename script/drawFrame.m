function drawFrame(H,scale)
% DRAWFRAME.m draws coordinate frame according to transformation H
%  Inputs:
%    H     ----- Homogeneous transform
%    scale ----- Scaling factor
%
% Author: Sipu Ruan, ruansp@jhu.edu, August 2017

if size(H,2) == 3
    X = scale * [0 1 0 ; 0 0 1];
    X = [X; 1 1 1 ];
    Xtrans = H * X;
    
    hold on;
    plot([Xtrans(1,1) Xtrans(1,2)], [Xtrans(2,1) Xtrans(2,2)], 'r');
    plot([Xtrans(1,1) Xtrans(1,3)], [Xtrans(2,1) Xtrans(2,3)], 'b');

elseif size(H,2) == 4
    X = scale * [0 1 0 0; 0 0 1 0; 0 0 0 1];
    X = [X; 1 1 1 1];
    Xtrans = H * X;
    
    hold on;
    plot3([Xtrans(1,1) Xtrans(1,2)], [Xtrans(2,1) Xtrans(2,2)], [Xtrans(3,1) Xtrans(3,2)], 'r');
    plot3([Xtrans(1,1) Xtrans(1,3)], [Xtrans(2,1) Xtrans(2,3)], [Xtrans(3,1) Xtrans(3,3)], 'g');
    plot3([Xtrans(1,1) Xtrans(1,4)], [Xtrans(2,1) Xtrans(2,4)], [Xtrans(3,1) Xtrans(3,4)], 'b');
end

end