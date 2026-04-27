
function [coord,elem,surf,neumann,Q]=mesh_P1_2D_n(N_h,x1,x2,x3,y1,y2)

% =========================================================================
%
%  This function creates triangular mesh for P1 elements
%
%  input data (integers):
%    N_h      - an integer defining a density of a uniform mesh
%    x1       - length of the body in front of the slope
%    x2       - length of the the slope in x-direction
%    x3       - length of the body behind the slope
%    y1       - hight of the body below the slope
%    y2       - height of the slope
%
%  output data:
%    coord     - coordinates of the nodes, size(coord)=(2,n_n) where n_n
%                is a number of nodes
%    elem      - array containing numbers of nodes defining each element,
%                size(elem)=(3,n_e), n_e = number of elements
%    surf      - array containing numbers of nodes defining each surface
%                element, size(surf)=(2,n_s), n_s = number of surface elements
%    neumann   - array containing numbers of nodes defining each
%                surface element, size(neuman)=(2,n_e_s). The surface 
%                is the following side of the body: (0,size_xy) x size_xy, 
%                where the nonhomogeneous Neumann boundary condition 
%                is considered.
%    Q         - logical array indicating the nodes where the Dirichlet
%                boundary condition is considered, size(Q)=(2,n_n)
%
% ======================================================================
%

%
% numbers of segments, nodes and elements
%

  N_x = (x1+x2+x3)*N_h; % number of segments in x direction
  N_y = (y1+y2)*N_h;    % number of segments in y direction
  N1_x=(x1+x2)*N_h;     % number of segments in x direction in front of the slope    
  N2_x=N_x-N1_x;        % number of segments in x direction right from the slope
  N1_y=y1*N_h;          % number of segments in y direction below the slope
  N2_y=N_y-N1_y;        % number of segments in y direction in the slope
  % 
  n_n = (N_x+1)*(N1_y+1)+(N1_x+1)*N2_y;     % number of nodes 
  n_e = 2*(N_x*N1_y+N1_x*N2_y);             % number of elements
 
%
% C - 2D auxilliary array that contains node numbers and that is important 
% for the mesh construction. Since the body is a union of two rectangles
% the array C also consists of two auxilliary 2D arrays, C1 and C2.
%
  C=zeros(N_x+1,N_y+1);
  C1=reshape(1:(N_x+1)*(N1_y+1),N_x+1,N1_y+1);
  C2=reshape(((N_x+1)*(N1_y+1)+1):n_n,N1_x+1,N2_y);
  C(1:(N_x+1),1:(N1_y+1))=C1;
  C(1:(N1_x+1),(N1_y+2):(N_y+1))=C2;
  
%
% coordinates of nodes
%
  % coordinates in directions x and y
  coord_x=linspace(0,x1+x2+x3,N_x+1);
  coord_y=linspace(0,y1+y2,N_y+1);
  %
  aux=((y1+y2-coord_y((N1_y+2):(N_y+1)))*x2+y2*x1)/(y2*(x1+x2));
  cy_x=repmat(coord_x(1:(N1_x+1)),N2_y,1).*repmat(aux',1,N1_x+1);
  cy_x=cy_x';

  % long 1D arrays containing coordinates of all nodes in x,y directions
  c_x=[repmat(coord_x,1,N1_y+1), cy_x(:)'];     
  c_y=[repmat(kron(coord_y(1:(N1_y+1)),ones(1,N_x+1)),1),...
           repmat(kron(coord_y((N1_y+2):(N_y+1)),ones(1,N1_x+1)),1)];  
  
  % the required array of coordinates, size(coord)=(2,n_n)
  coord=[c_x; c_y] ;       
 

% 
% construction of the array elem
%
  % ordering of the nodes creating the unit cube:
  %  V1 -> [0 0], V2 -> [1 0], V3 -> [1 1], V4 -> [0 1]
  %  V1,...,V4 are logical 2D arrays which enable to select appropriate
  %  nodes from the array C.

  V1=false(N_x+1,N_y+1);
  V1(1:N_x,1:N1_y)=1;
  V1(1:N1_x,(N1_y+1):N_y) =1; 
  %
  V2=false(N_x+1,N_y+1);
  V2(2:(N_x+1),1:N1_y)=1;
  V2(2:(N1_x+1),(N1_y+1):N_y) =1;
  %
  V3=false(N_x+1,N_y+1);
  V3(2:(N_x+1),2:(N1_y+1))=1;
  V3(2:(N1_x+1),(N1_y+2):(N_y+1)) =1;
  %
  V4=false(N_x+1,N_y+1);
  V4(1:N_x,2:(N1_y+1))=1;
  V4(1:N1_x,(N1_y+2):(N_y+1)) =1; 

  % used division of a rectangle into 2 triangles:   
  %   V1 V2 V4
  %   V2 V3 V4
  % size(aux_elem)=(2*3,n_e/2)
  aux_elem=[C(V1)'; C(V2)'; C(V4)';
            C(V2)'; C(V3)'; C(V4)' ];
        
  % the array elem, size(elem)=(3,n_e)          
  elem=reshape(aux_elem,3,n_e);     
 
%
% Surface of the body - the array "surf"
%
  
  % For each egde of the body, we define the restriction C_s of the array C 
  % and logical 2D arrays V1_s, V2_s which enable to select appropriate
  % vertices from the array C_s. We keep anticlockwise ordering of the 
  % vertices in order to detect the outward direction to the domain.
  
  % Edge 1: y=0 (the bottom of the body)
  C_s=zeros(N_x+1,1);
  C_s(:)=C(:,1);  
  V1_s=false(N_x+1,1);  V1_s(1:N_x,1)=1;
  V2_s=false(N_x+1,1);  V2_s(2:(N_x+1),1)=1;  
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf1=reshape(aux_surf,2,N_x);       
  
  % Edge 2: x=x1+x2+x3 (the right hand side of the body)
  C_s=zeros(N_y+1,1);
  C_s(:)=C(end,:);  
  V1_s=false(N_y+1,1);  V1_s(1:N1_y,1)=1;
  V2_s=false(N_y+1,1);  V2_s(2:(N1_y+1),1)=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf2=reshape(aux_surf,2,N1_y);       
  
  % Edge 3: y=y1+y2 (the top of the body)
  C_s=zeros(N_x+1,1);
  C_s(:)=C(:,end);  
  V1_s=false(N_x+1,1);  V1_s(2:(N1_x+1),1)=1;  
  V2_s=false(N_x+1,1);  V2_s(1:N1_x    ,1)=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf3=reshape(aux_surf,2,N1_x);       
  
  % Edge 4: x=0 (the left hand side of the body)
  C_s=zeros(N_y+1,1);
  C_s(:)=C(1,:);  
  V1_s=false(N_y+1,1);  V1_s(2:(N_y+1),1)=1;  
  V2_s=false(N_y+1,1);  V2_s(1:N_y    ,1)=1; 
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf4=reshape(aux_surf,2,N_y);     
  
  % Edge 5: y=y_1 (the top of the foundation)
  C_s=zeros(N_x+1,1);
  C_s(:)=C(:,N1_y+1);  
  V1_s=false(N_x+1,1);  V1_s((N1_x+2):(N_x+1),1)=1;
  V2_s=false(N_x+1,1);  V2_s((N1_x+1):(N_x),1)=1; 
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf5=reshape(aux_surf,2,N2_x);       
  
  % Edge 6: x=x1+x2 (the slope segment)
  C_s=zeros(N_y+1,1);
  C_s(:)=C(N1_x+1,:);  
  V1_s=false(N_y+1,1);  V1_s((N1_y+1):(N_y),1)=1;
  V2_s=false(N_y+1,1);  V2_s((N1_y+2):(N_y+1),1)=1;
  aux_surf=[C_s(V1_s)'; C_s(V2_s)' ];
  surf6=reshape(aux_surf,2,N2_y);  
  
  % the array "surf"
  surf = [surf1 surf2 surf3 surf4 surf5 surf6] ;
  
%
% Boundary conditions
%
  
  % nonhomogeneous Neumann boundary conditions on Face 3
  neumann=surf3 ; 

  % logical array indicating the nodes with the Dirichlet boundary cond.
  Q = coord>0 ;
  Q(1,(coord(1,:) >= x1+x2+x3-1e-10)) = 0;
  Q(2,(coord(2,:) == 0)) = 0;
end
