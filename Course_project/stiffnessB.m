function [Ke, Pe] = stiffnessB(x_a, y_b, t, E, v, q0, Int)

syms zeta_i eta_i zeta eta;
x = zeta*(x_a/2);
y = eta*(y_b/2);
G = E/2/(1+v);
D0 = E*t^3/12/(1 - v^2);
alpha = 5/6;
Ds = alpha*G*t*[1 0;0 1];
Db = D0*[1 v 0;v 1 0; 0 0 (1 - v)/2];

% Assigning the Gauss quadrature weights and points

wts_1 = 2;
wts_2 = [1 1; 1 1; 1 1; 1 1];

pts_1 = [0 0]; 
pts_2 = [-1/sqrt(3) -1/sqrt(3); 1/sqrt(3) -1/sqrt(3); 1/sqrt(3) 1/sqrt(3); -1/sqrt(3) 1/sqrt(3)];

% Defining the nodes

nodes = [-1 -1; 1 -1; 1 1; -1 1];

Ni = (1/4)*(1 + eta_i*eta)*(1 + zeta_i*zeta);

% Defining the Jacobian matrix 

J = [diff(x,zeta) diff(x,eta); diff(y,zeta) diff(y,eta)];

Ni_derivatives = J\[diff(Ni,zeta); diff(Ni,eta)];
Ni_x = Ni_derivatives(1);
Ni_y = Ni_derivatives(2);

var = [0 -Ni_x 0;
       0 0 -Ni_y;
       0 -Ni_y -Ni_x];
var2 = [Ni_x -Ni 0;
        Ni_y 0 -Ni];

% Populating the B matrix for both bending and shear part

Bb = sym('Bb', [3,12]);
Bs = sym('Bs',[2,12]);
for i=1:4
   Bb(:,3*i-2:3*i) = subs(var,[zeta_i eta_i],nodes(i,:));
end

for i=1:4
   Bs(:,3*i-2:3*i) = subs(var2,[zeta_i eta_i],nodes(i,:));
end

% Obtaining the stiffness matrix for both shear and bending

Kb_before_integral = transpose(Bb)*Db*Bb;
Kb_after_integral = Kb_before_integral;

Ks_before_integral = transpose(Bs)*Ds*Bs;
Ks_after_integral = Ks_before_integral;

% Evalutaing the integral as per Gauss quadrature rule depending on Reduced
% order integration or complete order integration

if Int(1) == 1
    for i=1:12
      for j=1:12 
          sum = 0;
          for k=1:1
             sum = sum + wts_1(1)^2*det(J)*subs(subs(Kb_before_integral(i,j),zeta, pts_1(k,1)),eta,pts_1(k,2));
          end
          Kb_after_integral(i,j) = sum;
      end
    end
elseif Int(1) == 2
    for i=1:12
      for j=1:12 
          sum = 0;
          for k=1:4
             sum = sum + wts_2(k,1)*wts_2(k,2)*det(J)*subs(subs(Kb_before_integral(i,j),zeta, pts_2(k,1)),eta,pts_2(k,2));
          end
          Kb_after_integral(i,j) = sum;
      end
    end
end

if Int(2) == 1
    for i=1:12
      for j=1:12 
          sum = 0;
          for k=1:1
             sum = sum + wts_1(1)^2*det(J)*subs(subs(Ks_before_integral(i,j),zeta, pts_1(k,1)),eta,pts_1(k,2));
          end
          Ks_after_integral(i,j) = sum;
      end
    end
elseif Int(2) == 2
    for i=1:12
      for j=1:12 
          sum = 0;
          for k=1:4
             sum = sum + wts_2(k,1)*wts_2(k,2)*det(J)*subs(subs(Ks_before_integral(i,j),zeta, pts_2(k,1)),eta,pts_2(k,2));
          end
          Ks_after_integral(i,j) = sum;
      end
    end
end

% Obtaining the shape functions

N = sym('N',[12,1]);

for i=1:4
N(3*i-2:3*i) = subs([Ni;0;0],[zeta_i eta_i],nodes(i,:));
end

% Obtaining the load vector

Pe = int(int(N*q0*det(J),eta,-1,1),zeta,-1,1);

Ke = Kb_after_integral + Ks_after_integral;
