function [LM,ID,IEN,x_a,y_b,n_eq,n_el,n_np] = direction(el_row,el_col,Edges,a,b)
%
% el_row : element number along x                                           
% el_col : element number along y
% Edges  : Boundary conditions for four edges
% n_el : total node numbers
% n_ep : nodes per element
% n_en : nodes per element
% x_a, y_b : element dimensions
% NND : Nodal constraint array
% ID : global D.O.F. array
% n_eq : total unconstrained D.O.F.
% IEN : element connectivity array
% LN : element to global direction array
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
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% Given the mesh Density calculate:
n_el = el_row*el_col;
n_np = (el_col + 1)*(el_row + 1);
n_en  =4; n_dof = 3;
x_a = a/el_col; y_b = b/el_row;
%
% Implement the specified boundary conditions into a Node/D.O.F ID Array
% For simply supported edge, only w=0 is implemented for equal
% applicability to Kirchoff and Mindlin Plates
%
% Specify displacements constraints (A's are node numbers)
for A=1:n_np
    NND(A,1)=A;
    NND(A,2:4)=1;
end

% Enter 0's for constrained d.o.f.s
if Edges(1)==1         % Edge 1 : simply supported
    for A=1:el_col+1
        NND(A,2)=0;
        NND(A,3)=0;
    end
elseif Edges(1)==2      % Edge 1 fixed
    for A=1:el_col+1
        NND(A,2:4)=0;
    end
end

if Edges(3)==1          % Edge 3: simply supported
    for A=n_np-el_col:n_np
        NND(A,2)=0;
        NND(A,3)=0;
    end 
elseif Edges(3)==2      % Edge 3: fixed
    for A=n_np-el_col:n_np
        NND(A,2:4)=0;
    end
end

if Edges(4)==1           % Edge 4: simply supported
    for A=1:el_col+1:n_np-el_col
        NND(A,2)=0;
        NND(A,4)=0;
    end
elseif Edges(4)==2        % Edge 4: fixed
    for A=1:EL_COL+1:N_np-el_col
        NND(A,2:4)=0;
    end
end

if Edges(2)==1             % Edge 2: simply supported
    for A=el_col+1:el_col+1:n_np
        NND(A,2)=0;
        NND(A,4)=0;
    end
elseif Edges(2)==2         % Edge 2: fixed
    for A=el_col+1:el_col+1:n_np
        NND(A,2:4)=0;
    end
end

%
% Build ID Array
j=1; jj=1;

NN=NND(:,1);
%
for A=1:n_np
    for i=1:n_dof
        if (A==NN(jj,1) & (NND(jj,i+1)==0))
            ID(i,A)=0;
            if (i==n_dof)
                jj=jj+1;
            end
        else
            ID(i,A)=j;
            j=j+1;
            if (i==n_dof)
                jj=jj+1;
            end
        end
    end
end

n_eq = j-1;
ID;
%
% Build Element Connectivity Array
for nc=1:el_col
    IEN(1,nc)=nc;
    IEN(2,nc) = nc+1;
    IEN(3,nc) = IEN(2,nc) + (el_col + 1);
    IEN(4, nc) = IEN(1,nc) + (el_col +1);
end

for nc=el_col+1:n_el
    for nr=1:n_en
        IEN(nr,nc) = IEN(nr,nc-el_col) + (el_col + 1);
    end
end

IEN
%
% Build element stiffness array
for e=1:n_el
    for i=1:n_dof
        for a=1:n_en
            p=n_dof*(a-1)+i;
            LM(p,e)=ID(i,IEN(a,e));
        end
    end
end
LM;