function [X,Y,U,x_a,y_b,n_el,n_np,IEN,d_max] = main(el_row,el_col,Element,Int,a,b,t,E,v,q0,Edges)
%
% Direction array, node properties
[LM,ID,IEN,x_a,y_b,n_eq,n_el,n_np] = direction(el_row,el_col,Edges,a,b);
%
% Stiffness matrix, load vector, and displacements
[d] = analysis(Element,Int,x_a,y_b,t,E,v,q0,LM,n_eq,n_el);
%
% Results
[X,Y,U,d_max] = postprocess2(d,x_a,y_b,ID,el_row,el_col,n_np);