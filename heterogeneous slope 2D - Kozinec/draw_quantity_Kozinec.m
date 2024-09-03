function draw_quantity_Kozinec(coord,elem,U,Q_node)

% =========================================================================
%
%  This function depicts prescribed nodal quantity
%
%  input data:
%    coord - coordinates of the nodes, size(coord)=(2,n_n) where n_n is a
%            number of nodes
%    elem - n_p x n_e array containing numbers of nodes defining each
%           element, n_e = number of elements
%    U - nodal displacements, size(U)=(2,n_n) to catch deformed shape
%        if the deformed shape is not required then set 0*U
%    Q_node - prescribed nodal quantity, size(Q_node)=(1,n_n)
%
% ======================================================================
%

  figure;
  hold on;
  
  patch('Faces',elem(1:3,:)','Vertices',coord'+U',...
        'FaceVertexCData',Q_node','FaceColor','interp','EdgeColor','none'); 
    

  % na i-tom riadku je index a suradnice i-teho vrcholu
  C = load('kozinec2_Coord.txt'); C=C(:,2:3);

  % indexy dublovanych vrcholov:
  % I(i,1)  je nahradeny indexom  I(i,2)
  I = load('kozinec2_Dupl.txt');

  % indexy orientovanych hranic oblasti; oblasti su separovane hodnotou 0 
  D = load('kozinec2_Dom0.txt'); m=size(D,2);

  % z vektoru D nachadzame indexy vrcholov jednotlivych oblasti:
  % J(i,1:2) obsahuje pociatocny a koncovy index hranice i-tej oblasti
  iB=1; iE=iB; nD=0;
  while iE<m
    nD=nD+1; J(nD,1)=iB;
    while D(iE)>0, iE=iE+1; end
    J(nD,2)=iE-1;
    iB=iE+1; iE=iB;
  end
  %

 % figure(1); hold on; axis off;
  for i=1:nD
    ID=D(J(i,1):J(i,2)); ID=[ID ID(1)];
    plot(C(ID,1),C(ID,2),'black');
  end
  %
  n=size(C,1);
  II=[1:n;1:n]; I=I'; II(2,I(1,:))=I(2,:);

%   for i=1:n
%     j=II(2,i); text(C(j,1),C(j,2),num2str(j));
%   end  
  
  box on
  axis equal;
  hold off;
  axis off;
end