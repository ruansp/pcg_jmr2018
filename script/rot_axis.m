% rotating along an axis defined by u, with angle theta
function R = rot_axis(u, theta)
u = u./norm(u,2);
S = [    0  u(3) -u(2);
      -u(3)   0   u(1);
      u(2) -u(1)   0  ];
R = eye(3) + sin(theta)*S + (1-cos(theta))*S^2;
end