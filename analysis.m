%
%Element : element type (1 for Kirchoff plate, 2 for Mindlin plate)
% Int : Integration scheme (only applicable for Mindlin plate)
% x_a, y_b : element dimensions
% t : plate thickness
% E : Young's modulus
% v : Poisson's ratio
% D : flexural rigidity
% q0 : unoformly distributed load
% n_el : total element numbers
% n_eq : total unconstrained D.O.F.
% LM : element to global direction array
% Ke, Pe : element to global direction array
% K, P : global stiffness matrix and load vector
% d : solved D.O.F. array
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [d]=analysis(Element,Int,x_a,y_b,t,E,v,q0,LM,n_eq,n_el)
K=zeros(n_eq); P=zeros(n_eq,1);
% compute the element stiffness matrix (Ke) and load vector (Pe)
% (computation performed only once because of uniform mesh)
if Element==1
    [Ke, Pe] = stiffnessA(x_a,y_b,t,E,v,q0);
elseif Element==2
    [Ke, Pe] = stiffnessB(x_a, y_b, t, E, v, q0, Int);
end
%
% Construct the global stiffness matrix and load vector
j=1;
for e=1:n_el
    gl=LM(:,e);
    i=find(gl~=0);
    I=gl(i);
    K(I,I)=K(I,I) + Ke(i,i);
    P(I) = P(I) + Pe(i);
    j=j+1;
end
%
% Compute the Deformations
size(K)
size(P);
d = K\P;