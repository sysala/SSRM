function [HatP,DHatP1,DHatP2] = local_basis_volume_2D(elem_type, Xi)

%--------------------------------------------------------------------------
% This function evaluates local basis functions and their derivatives at
% prescribed quadrature points depending on a chosen finite elements.
%
%  input data: 
%       elem_type - the type of Lagrange finite elements
%       Xi - coordinates of the quadrature points, size(Xi)=(2,n_q)
%
%  output data:
%       HatP - values of basis functions at the quadrature points,
%              size(HatP)=(n_p,n_q)
%       DHatP1 - derivatives of basis functions at the quadrature points 
%                in the direction xi_1, size(HatP)=(n_p,n_q)
%       DHatP2 - derivatives of basis functions at the quadrature points 
%                in the direction xi_2, size(DHatP2)=(n_p,n_q)
%       n_p - number of basis functions
%       n_q - number of integration points within one element
%--------------------------------------------------------------------------

xi_1 = Xi(1,:); xi_2 = Xi(2,:);

switch(elem_type)
  case 'P1'
    % - the reference triangle with coordinates:
    %               [0,0], [1,0], [0,1]
    % - n_p=3, n_q=length(xi_1)
    HatP = [1-xi_1-xi_2; xi_1; xi_2] ;
    DHatP1 = [ -1; 1; 0];
    DHatP2 = [ -1; 0; 1];
    
  case 'P2'
    % - the reference triangle with coordinates:
    %               [0,0], [1,0], [0,1]
    % - coordinates of midpoints:
    %             [1/2,1/2], [0,1/2], [1/2,0], 
    % - n_p=6, n_q=length(xi_1)
    % - barycentric coordinates: xi_0, xi_1, xi_2, where
    xi_0 = 1-xi_1-xi_2;
    n_q = length(xi_1);
    HatP = [ xi_0.*(2*xi_0-1); xi_1.*(2*xi_1-1); xi_2.*(2*xi_2-1);
                 4*xi_1.*xi_2;     4*xi_0.*xi_2;     4*xi_0.*xi_1 ];
    DHatP1 = [ -4*xi_0+1; 4*xi_1-1; zeros(1,n_q); 
                  4*xi_2;  -4*xi_2;  4*(xi_0-xi_1) ];
    DHatP2 = [ -4*xi_0+1; zeros(1,n_q); 4*xi_2-1; 
                  4*xi_1;  4*(xi_0-xi_2);  -4*xi_1 ];      
              
    case 'P4'
    % - the reference triangle with coordinates:
    %               [0,0], [1,0], [0,1]
    % - coordinates of midpoints:
    %   [1/2,0], [1/2,1/2], [0,1/2], [1/4,0], [3/4,0], [3/4,1/4]
    %   [1/4,3/4], [0,3/4], [0,1/4], [1/4,1/4], [1/2,1/4], [1/4,1/2]
    % - n_p=15, n_q=length(xi_1)
    % - barycentric coordinates: xi_0, xi_1, xi_2, where
    xi_0 = 1-xi_1-xi_2;
    n_q = length(xi_1);
    HatP = [ xi_0.*(4*xi_0-1).*(4*xi_0-2).*(4*xi_0-3)/6
             xi_1.*(4*xi_1-1).*(4*xi_1-2).*(4*xi_1-3)/6
             xi_2.*(4*xi_2-1).*(4*xi_2-2).*(4*xi_2-3)/6
             4*xi_0.*xi_1.*(4*xi_0-1).*(4*xi_1-1)
             4*xi_1.*xi_2.*(4*xi_1-1).*(4*xi_2-1)
             4*xi_0.*xi_2.*(4*xi_0-1).*(4*xi_2-1)
             8*xi_0.*xi_1.*(4*xi_0-1).*(4*xi_0-2)/3
             8*xi_0.*xi_1.*(4*xi_1-1).*(4*xi_1-2)/3
             8*xi_1.*xi_2.*(4*xi_1-1).*(4*xi_1-2)/3
             8*xi_1.*xi_2.*(4*xi_2-1).*(4*xi_2-2)/3
             8*xi_0.*xi_2.*(4*xi_2-1).*(4*xi_2-2)/3
             8*xi_0.*xi_2.*(4*xi_0-1).*(4*xi_0-2)/3
             32*xi_0.*xi_1.*xi_2.*(4*xi_0-1)
             32*xi_0.*xi_1.*xi_2.*(4*xi_1-1)
             32*xi_0.*xi_1.*xi_2.*(4*xi_2-1) ];
    DHatP1 = [ -((4*xi_0-1).*(4*xi_0-2).*(4*xi_0-3)+4*xi_0.*(4*xi_0-2).*(4*xi_0-3)+4*xi_0.*(4*xi_0-1).*(4*xi_0-3)+4*xi_0.*(4*xi_0-1).*(4*xi_0-2))/6
                ((4*xi_1-1).*(4*xi_1-2).*(4*xi_1-3)+4*xi_1.*(4*xi_1-2).*(4*xi_1-3)+4*xi_1.*(4*xi_1-1).*(4*xi_1-3)+4*xi_1.*(4*xi_1-1).*(4*xi_1-2))/6
                zeros(1,n_q)
                4*(-xi_1.*(4*xi_0-1).*(4*xi_1-1)+xi_0.*(4*xi_0-1).*(4*xi_1-1)-4*xi_0.*xi_1.*(4*xi_1-1)+4*xi_0.*xi_1.*(4*xi_0-1))
                4*(xi_2.*(4*xi_1-1).*(4*xi_2-1)+4*xi_1.*xi_2.*(4*xi_2-1))
                4*(-xi_2.*(4*xi_0-1).*(4*xi_2-1)-4*xi_0.*xi_2.*(4*xi_2-1))
                8*(-xi_1.*(4*xi_0-1).*(4*xi_0-2)+xi_0.*(4*xi_0-1).*(4*xi_0-2)-4*xi_0.*xi_1.*(4*xi_0-2)-4*xi_0.*xi_1.*(4*xi_0-1))/3
                8*(-xi_1.*(4*xi_1-1).*(4*xi_1-2)+xi_0.*(4*xi_1-1).*(4*xi_1-2)+4*xi_0.*xi_1.*(4*xi_1-2)+4*xi_0.*xi_1.*(4*xi_1-1))/3
                8*(xi_2.*(4*xi_1-1).*(4*xi_1-2)+4*xi_1.*xi_2.*(4*xi_1-2)+4*xi_1.*xi_2.*(4*xi_1-1))/3
                8*xi_2.*(4*xi_2-1).*(4*xi_2-2)/3
               -8*xi_2.*(4*xi_2-1).*(4*xi_2-2)/3
                8*(-xi_2.*(4*xi_0-1).*(4*xi_0-2)-4*xi_0.*xi_2.*(4*xi_0-2)-4*xi_0.*xi_2.*(4*xi_0-1))/3
                32*(-xi_1.*xi_2.*(4*xi_0-1)+xi_0.*xi_2.*(4*xi_0-1)-4*xi_0.*xi_1.*xi_2)
                32*(-xi_1.*xi_2.*(4*xi_1-1)+xi_0.*xi_2.*(4*xi_1-1)+4*xi_0.*xi_1.*xi_2)
                32*(-xi_1.*xi_2.*(4*xi_2-1)+xi_0.*xi_2.*(4*xi_2-1)) ];
    DHatP2 = [ -((4*xi_0-1).*(4*xi_0-2).*(4*xi_0-3)+4*xi_0.*(4*xi_0-2).*(4*xi_0-3)+4*xi_0.*(4*xi_0-1).*(4*xi_0-3)+4*xi_0.*(4*xi_0-1).*(4*xi_0-2))/6 
                zeros(1,n_q) 
                ((4*xi_2-1).*(4*xi_2-2).*(4*xi_2-3)+4*xi_2.*(4*xi_2-2).*(4*xi_2-3)+4*xi_2.*(4*xi_2-1).*(4*xi_2-3)+4*xi_2.*(4*xi_2-1).*(4*xi_2-2))/6
                4*(-xi_1.*(4*xi_0-1).*(4*xi_1-1)-4*xi_0.*xi_1.*(4*xi_1-1))
                4*(xi_1.*(4*xi_1-1).*(4*xi_2-1)+4*xi_1.*xi_2.*(4*xi_1-1))
                4*(-xi_2.*(4*xi_0-1).*(4*xi_2-1)+xi_0.*(4*xi_0-1).*(4*xi_2-1)-4*xi_0.*xi_2.*(4*xi_2-1)+4*xi_0.*xi_2.*(4*xi_0-1))
                8*(-xi_1.*(4*xi_0-1).*(4*xi_0-2)-4*xi_0.*xi_1.*(4*xi_0-2)-4*xi_0.*xi_1.*(4*xi_0-1))/3
               -8*xi_1.*(4*xi_1-1).*(4*xi_1-2)/3
                8*xi_1.*(4*xi_1-1).*(4*xi_1-2)/3
                8*(xi_1.*(4*xi_2-1).*(4*xi_2-2)+4*xi_1.*xi_2.*(4*xi_2-2)+4*xi_1.*xi_2.*(4*xi_2-1))/3
                8*(-xi_2.*(4*xi_2-1).*(4*xi_2-2)+xi_0.*(4*xi_2-1).*(4*xi_2-2)+4*xi_0.*xi_2.*(4*xi_2-2)+4*xi_0.*xi_2.*(4*xi_2-1))/3
                8*(-xi_2.*(4*xi_0-1).*(4*xi_0-2)+xi_0.*(4*xi_0-1).*(4*xi_0-2)-4*xi_0.*xi_2.*(4*xi_0-2)-4*xi_0.*xi_2.*(4*xi_0-1))/3
                32*(-xi_1.*xi_2.*(4*xi_0-1)+xi_0.*xi_1.*(4*xi_0-1)-4*xi_0.*xi_1.*xi_2)
                32*(-xi_1.*xi_2.*(4*xi_1-1)+xi_0.*xi_1.*(4*xi_1-1))
                32*(-xi_1.*xi_2.*(4*xi_2-1)+xi_0.*xi_1.*(4*xi_2-1)+4*xi_0.*xi_1.*xi_2) ];   
 
  otherwise; disp('Bad choise of element type');
end
