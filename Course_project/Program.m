%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      PROGRAM.m                            %
% This is a finite-element program for plate analysis using %
% either 4-node Kirchoff or Mindlin plate elements.         % 
%                       developed by J. Zhang 12/1.2007     %
% Input Variables:                                          %
% el_row : element numbers along x                          %
% el_col : element numbers along y                          %
% a      : plate length                                     %
% b      : plate width                                      %
% t      : plate thickness                                  %
% E      : Young's modulus                                  %
% v      : Poisson's ratio                                  %
% D      : flexural rigidity                                %
% q0     : uniformly distributed loading                    %
% Element: Kirchoff plate (Element-1)                       %
%          Mindlin plate (element-2)                        %
% Edges  : boundary condition (0-free, 1-simply supported   %
%          2-fixed                                          %
% Output Variables:                                         %
% d_max : maxim displacement                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
clear all; close all; clc
tic
%
% Input rectangular plate mesh density
el_row = 20; el_col = 20;

% Element type (1-Kirchoff; 2-Reisner/Mindlin and integration scheme for
% R/M only)

Element=2; Int=[2 1];

%Overall structure geometry
a = 4; b = 4; t=1;

%Material properties
E = 200e9; v = 0.3;
D = E*t^3/(12*(1 - v^2));

%Loading
q0 = 30e3;

%Boundary conditions (0-free; 1-simple support; 2-fixed)
Edges = [0 2 2 1]; %[0 2 2 1] for (free; fixed; fixed; simply-supported)
%
%      y-axis
%      ^
%      |    edge 3
%      |--------------|
%      |              |
%      |              |edge 2
%edge 4|              |
%      |              | 
%      |              |
%      |--------------|---> x axis
%          edge 1

disp('---------------------FEM Solution-----------------------')
[X,Y,U,x_a,y_b,n_el,n_np,IEN,d_max] = main(el_row, el_col, Element, Int, a, b,t,E,v,q0,Edges);
format long;
resp_const = d_max*D/(q0*a^4);

% Exact solution of Mindlin Plate with all four boundary simply-supported
w_exact_M = 0; w_exact_K = 0; w_exact_M2 = 0;
alpha = 5/6;
G = E/2/(1 + v);
x0= a/2; y0=b/2;
for i=1:2:31
    for j=1:2:31
        qmn = 16*q0/pi^2/i/j;
        w_k = qmn/D/pi^4*sin(i*pi*x0/a)*sin(j*pi*y0/b)/(i^2/a^2+j^2/b^2)^2;
        w_s = qmn/(alpha*G*t)/pi^2*sin(i*pi*x0/a)*sin(j*pi*y0/b)/(i^2/a^2 + j^2/b^2);
        w_exact_M = w_exact_M + w_k + w_s;
        w_exact_K = w_exact_K + w_k;
        w_exact_M2 = w_exact_M2 + w_k + 1/alpha/G/t/pi^2*qmn/(i^2/a^2 + j^2/b^2)*sin(i*pi*x0/a)*sin(j*pi*y0/b);
    end
end

disp('---Exact Kirchoff sol for Simply-supported plate')
w_exact_K
factor_K = w_exact_K*D/q0/a^4;
disp('----Mindlin Sol based on Conjugate Plate analogy---')
w_exact_M
factor_M = w_exact_M*D/q0/a^4;
w_exact_M2;
toc