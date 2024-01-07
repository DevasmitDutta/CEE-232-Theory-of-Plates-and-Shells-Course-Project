% Subroutine to post process the results and plot
% 
function [X,Y,U,d_max]=postprocess2(d,x_a,y_b,ID,el_row,el_col,n_np);
%
% Transfer [d] into nodal displacement array [U] 
j=1;
for A=1:n_np
    for i=1:3
        if ID(i,A)==0
            U(j)=0;
            j=j+1;
        else
            U(j)=d(ID(i,A),1);
            j=j+1;
        end
    end
end
%
% Vertical displacements @ each node
% n_pr-nodes per row; n_pc-nodes per column
n_pr = el_col + 1; n_pc = el_row+1;
x=0;y=0;
j=1;
for nr=1:n_pc
    for nc=1:n_pr
        X(nr,nc)=x;
        Y(nr,nc)=y;
        U_Z(nr,nc)=U(j);
        j=j+3;
        x=x+x_a;
    end
    x=0;
    y=y+y_b;
end
%
% PLOT DEFORMED SHAPE
figure(1)
mesh(X,Y,-U_Z)
title('Deformed Shape of Plate')
xlabel('X')
ylabel('Y')
zlabel('Deflection')
figure(2)
surf(X,Y,-U_Z)
colormap hsv
colorbar

shading interp
title('Deformed shape of Plate')
xlabel('X')
ylabel('Y')
zlabel('Deflection')
% Find max deflection 
format long
d_max = max(max(U_Z(:,:)))