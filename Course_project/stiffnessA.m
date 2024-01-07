function [Ke, Pe] = stiffnessA(x_a,y_b,t,E,v,q0)
syms x y;
D0 = E*t^3/12/(1 - v^2);

D = D0*[1 v 0;v 1 0; 0 0 (1 - v)/2];
P_bar = [1 x y x^2 x*y y^2 x^3 x^2*y x*y^2 y^3 x^3*y x*y^3];

C_bar = [P_bar; diff(P_bar,x); diff(P_bar,y)];

nodes = [0 0; x_a 0; x_a y_b; 0 y_b];

C = zeros(12,12); 

% Populating the matrix C by evaluating the C_bar at each ode

for i =1:4
    C(3*i-2:3*i,:) = subs(C_bar,[x y],nodes(i,:));
end

% Obtaining the B matrix 

B = [-diff(diff(P_bar,x),x); -diff(diff(P_bar,y),y); -2*diff(diff(P_bar,x),y)]/C;

% Obtaining the element level stiffness matrix

Ke = int(int(transpose(B)*D*B,x,0,x_a),y,0,y_b);

% Obtaining the Shape-functions

N = P_bar/C;

% Obtaining the Load Vector

Pe = int(int(transpose(N)*q0,x,0,x_a),y,0,y_b);
