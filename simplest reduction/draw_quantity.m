function draw_quantity(coord,elem,U,Q_node,x1,x2,x3,y1,y2)

% =========================================================================
%
%  This function depicts prescribed nodal quantity
%
%  input data:
%    coord    - coordinates of the nodes, size(coord)=(2,n_n) where n_n is
%               a number of nodes
%    elem     - array containing numbers of nodes defining each element, 
%               size(elem)=(n_p,n_e), n_e = number of elements (triangles)
%    U - nodal displacements, size(U)=(2,n_n) to catch deformed shape
%        if the deformed shape is not required then set 0*U
%    Q_node - prescribed nodal quantity, size(Q_node)=(1,n_n)
%    x1 - length of the body in front of the slope
%    x2 - length of the slope in x-direction
%    x3 - length of the body behind the slope
%    y1 - hight of the body below the slope
%    y2 - height of the slope
%
% ======================================================================
%

  figure;
  hold on;
  
  % visualization of the quantity
  patch('Faces',elem(1:3,:)','Vertices',coord'+U',...
        'FaceVertexCData',Q_node','FaceColor','interp','EdgeColor','none'); 
 
%  alpha(s,.5);
%  colorbar;
  
  % undeformed shape of the body
  plot([0,x1+x2+x3],[0,0])
  plot([x1+x2+x3,x1+x2+x3],[0,y1])
  plot([x1+x2+x3,x1+x2],[y1,y1])
  plot([x1+x2,x1],[y1,y1+y2])
  plot([x1,0],[y1+y2,y1+y2]) 
  plot([0,0],[y1+y2,0]) 
  
  %
  box on
  view(2);
  axis equal;
  hold off;
  axis off;
end