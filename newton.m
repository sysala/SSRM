function [U_it, zeta_it, it, criterion]=...
               newton(U_ini,alpha,WEIGHT,B,iD,jD,K_elast,Q,f,...
                      c_bar,sin_phi,shear,bulk,lame)
                               
%--------------------------------------------------------------------------
% The Newton method for solution of the system 
%              find U, zeta:   F_alpha(U)=zeta*f_t, f_t'*U=1
%
% Input data:
%   U_ini - initial choice of U
%   alpha - prescribed value of the parameter alpha
%   WEIGHT - weight coefficients of integration points, size(WEIGHT)=(1,n_int)
%   B - the strain-displacement matrix, size(B)=(3*n_int,2*n_n)
%   D_elast=sparse(iD,jD,vD_elast) - elastic stress-strain matrix,
%                                    size(D_elast)=(3*n_int,3*n_int)
%   K_elast=B'*D_elast*B - the elastic stiffness matrix, size(K)=(2*n_n, 2*n_n)
%   Q - logical 2xn_n array restricting nodes with Dirichlet boundary cond.
%   f - vector of external load
%   c_bar,sin_phi,shear,bulk,lame - material parameters at integration points
%
% Output data:
%   U_it, zeta_it - approximation of the solution
%   it - number of Newton's iteration
%   criterion - value of stopping criterion
%
%--------------------------------------------------------------------------

%
% Initialization     
%
                       
  n_n=size(U_ini,2);
  n_int=length(WEIGHT);
  V = zeros(2,n_n);
  W = zeros(2,n_n);
  F = zeros(2,n_n) ;   % vector of internal forces
  E = zeros(3,n_int);  % strain tensors at integration points               
  U_it=U_ini;          % initial displacements

%  
% Newton's solver (the semismooth Newton method)
%
      
  it=0;               % iteration number
  while true         
      
     it=it+1;  
     
     % constitutive operator and its derivative
     E(:) = B*U_it(:) ;   % strain at integration points
     [S,DS]=constitutive_problem(E,alpha,c_bar,sin_phi,...
                                       shear,bulk,lame);
                          % solution of the constitutive problem
     
     % tangential stiffness matrix
     vD = repmat(WEIGHT,9,1).*DS ;      
     D_p = sparse( iD(:),jD(:),vD(:), 3*n_int,3*n_int ) ;   
     K_tangent = B'*D_p*B;   
     K_tangent=(K_tangent+K_tangent')/2;
 
     % vector of internal forces
     F(:) = B'*reshape(repmat(WEIGHT,3,1).*S(1:3,:), 3*n_int,1) ;      
     zeta_0=f(Q)'*F(Q)/norm(f(Q))^2;
     
     % Newton's increment
     W(Q) = K_tangent(Q,Q)\f(Q); 
     V(Q) = K_tangent(Q,Q)\(zeta_0*f(Q)-F(Q)); 
     zeta_it=zeta_0-f(Q)'*V(Q)/(f(Q)'*W(Q));
     dU=V+(zeta_it-zeta_0)*W;       
     
     % damping parameter
     omega=damping(U_it,dU,F,B,WEIGHT,alpha,c_bar,sin_phi,shear,bulk,lame);
     
     % recomputation with a regularized stiffness matrix when omega=0
     if omega==0 
         r=0.01;
         K_tangent = r*alpha*K_elast+(1-r)*K_tangent;
         W(Q) = K_tangent(Q,Q)\f(Q); 
         V(Q) = K_tangent(Q,Q)\(zeta_0*f(Q)-F(Q)); 
         zeta_it=zeta_0-f(Q)'*V(Q)/(f(Q)'*W(Q));
         dU=V+(zeta_it-zeta_0)*W;     
         omega=1;
     end
     
     % stopping criterion 
     U_new= U_it + dU ;
     q1 = sqrt( dU(:)'*K_elast*dU(:) ) ;
     q2 = sqrt(  U_it(:)'*K_elast*U_it(:)  ) ;
     q3 = sqrt( U_new(:)'*K_elast*U_new(:) ) ;
     criterion = q1/(q2+q3);         
  
     % test on stopping criterion
     if  criterion < 1e-12
        U_it=U_new; 
%         fprintf('     Newton method converges: plastic integration points: %d (of %d), ',n_plast,n_int); 
%         fprintf(' number of iteration=%d  ',it); 
%         fprintf(' stopping criterion=%e  ',criterion); 
%         fprintf('\n'); 
        break
     end      
     
     % next iteration
     U_it= U_it + omega*dU ;        
      
     % test on number of iteration
     if  it == 25
        fprintf('     Newton solver converges slowly: stopping criterion=%e  ',criterion)
        fprintf('\n'); 
        break
     end     
     
  end % true     
  
end % function