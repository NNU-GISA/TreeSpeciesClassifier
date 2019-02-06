X =[2 4 3;5 9 4; 7 9 6;12 19 2; 2 4 6; 5 8 7; 6 8 3 ; 8 1 3; 9 8 3;9 8 7; 7 6  5;];

plot3(X(:,1),X(:,2),X(:,3),'*')
hold on;
[x0n, an, phin, rn, d, sigmah, conv, Vx0n,...
    Van, uphin, urn, GNlog, a, R0, R] = lscone(X, [5 2 6]', [4 1 4]', 0.5, 4, 0.01, 0.01, 1);


r = rn;
h = 15;
m = h/r;
[R,A] = meshgrid(linspace(0,r,11),linspace(0,2*pi,41));
X = R .* cos(A);
Y = R .* sin(A);
Z = m*R;
% Cone around the z-axis, point at the origin
mesh(X,Y,Z)

xlabel('x')
ylabel('y')
zlabel('z')



% hold on
% axis equal
% axis([-3 3 -3 3 0 3])

% phi = -pi/3;
% X1 = X*cos(phi) - Z*sin(phi);
% Y1 = Y;
% Z1 = X*sin(phi) + Z*cos(phi);
% % Previous cone, rotated by angle phi about the y-axis
% mesh(X1,Y1,Z1)
% 
% theta = pi/4;
% X2 = X1*cos(theta) - Y1*sin(theta);
% Y2 = X1*sin(theta) + Y1*cos(theta);
% Z2 = Z1;
% Second cone rotated by angle theta about the z-axis
%mesh(X2,Y2,Z2)