function [U_it, lagr_it, it, criterion]=...
                 newton(tol1,it_newt_max,it_damp_max,U_ini,alpha,lambda,lagr,...
                        WEIGHT,B,K_elast,Q,f,c_bar,sin_phi,shear,bulk,lame)
                               
%--------------------------------------------------------------------------
% The Newton method with damping for solution of the system 
%              find U:   F_{alpha}(U)=f
%
% Input data:
%   U_ini - initial choice of U
%   alpha - prescribed value of the parameter alpha
%   WEIGHT - weight coefficients at integration points, size(WEIGHT)=(1,n_int)
%   B - the strain-displacement matrix, size(B)=(3*n_int,2*n_n)
%   K_elast - the elastic stiffness matrix, size(K_elast)=(2*n_n, 2*n_n)
%   Q - logical 2xn_n array restricting nodes with Dirichlet boundary cond.
%   f - vector of external load
%   shear,bulk,lame - elastic material parameters (1 x n_int arrays)
%   c_bar,sin_phi - 1 x n_int arrays with heterogeneous parameters
%
% Output data:
%   U_it - the solution of the nonlinear system of equations
%   lagr_it - value of the cost function
%   it - number of Newton's iteration
%   criterion - value of stopping criterion
%
%--------------------------------------------------------------------------

%
% Initialization     
%
                       
  n_n=size(U_ini,2);     % number of nodes
  n_int=length(WEIGHT);  % number of integration points
  dU = zeros(2,n_n);     % Newton's increment 
  F = zeros(2,n_n) ;     % vector of internal forces
  E = zeros(3,n_int);    % strain tensors at integration points               
  U_it=U_ini;            % initial displacements

  % constitutive operator and its derivative
  E(:) = B*U_it(:) ;   % strain at integration points
  [S,DS]=constitutive_problem(E,alpha,c_bar,sin_phi,shear,bulk,lame);
  % vector of internal forces
  F(:) = B'*reshape(repmat(WEIGHT,3,1).*S(1:3,:), 3*n_int,1) ;      

%  
% Newton's solver (the semismooth Newton method)
%
      
  it=0;               % iteration number
  while true         
      
     it=it+1;       
         
     % assembling of the tangential stiffness matrix
     AUX=reshape(1:3*n_int,3,n_int);
     iD=repmat(AUX,3,1); 
     jD=kron(AUX,ones(3,1));
     vD = repmat(WEIGHT,9,1).*DS ;      
     D_p = sparse( iD(:),jD(:),vD(:), 3*n_int,3*n_int ) ;   
     K_tangent = B'*D_p*B;    
         
     % Newton's increment
     dU(Q) = K_tangent(Q,Q)\(f(Q)-F(Q)); 
     
     % damping parameter
     omega=damping(it_damp_max,U_it,dU,F,f,B,WEIGHT,alpha,c_bar,sin_phi,shear,bulk,lame);
     
     % recomputation with a regularized stiffness matrix when omega=0
     if omega==0 
         r=0.01;
         K_tangent = r*alpha*K_elast+(1-r)*K_tangent;
         dU(Q) = K_tangent(Q,Q)\(f(Q)-F(Q)); 
         omega=1;
     end
     
     % next iteration
     U_it= U_it + omega*dU ; 
     E(:) = B*U_it(:) ;   % strain at integration points    
                            
     % Lagrangian value
     Psi= potential(E,alpha,c_bar,sin_phi,shear,bulk,lame); 
     lagr_it=lambda+WEIGHT*Psi'-f(Q)'*U_it(Q);
    
     % constitutive operator and its derivative
     [S,DS]=constitutive_problem(E,alpha,c_bar,sin_phi,...
                                                   shear,bulk,lame);
     % vector of internal forces
     F(:) = B'*reshape(repmat(WEIGHT,3,1).*S(1:3,:), 3*n_int,1) ; 
 
     % stopping criterion and its test
     criterion = norm(F(Q)-f(Q))/norm(f(Q));         
     if  criterion < tol1
        fprintf('     Newton method converges:  '); 
        fprintf(' number of iteration=%d  ',it); 
        fprintf(' stopping criterion=%e  ',criterion); 
        fprintf('\n'); 
        break
     end      
     
     % test on Lagrangian value
     if  lagr_it < min(lagr,0)
        fprintf('     Lagrangian value significantly decreased: lagr_it=%e  ',lagr_it)
        fprintf('\n'); 
        break
     end  
          
     % test on number of iteration
     if  it == it_newt_max
        fprintf('     Newton solver converges slowly: stopping criterion=%e  ',criterion)
        fprintf('\n'); 
        break
     end     
     
  end % true     
  
end % function