function [U_it, t_it, flag_N]=...
               newton_ALG4(U_ini,t_ini,omega,it_newt_max,it_damp_max,...
                           tol,r_min,r_damp,WEIGHT,B,K_elast,Q,f,...
                           c_bar,sin_phi,shear,bulk,lame)
                               
%--------------------------------------------------------------------------
% The Newton method with damping and regularization for solving the system 
%              find U, t:   F(U)=t*f, f_t'*U=omega
%
% Input data:
%   U_ini - initial choice of U
%   t_ini - initial choice of t
%   omega - prescribed value of the parameter omega
%   it_newt_max - number of Newton's iterations
%   it_damp_max - number of iterations within line search
%   tol - relative tolerance for Newton's solvers
%   r_min - basic regularization of the stiffness matrix
%   r_damp - regularization of the stiffness matrix if
%            it is impossible to find descent direction
%   WEIGHT - weight coefficients of integration points
%   B - the strain-displacement matrix
%   K_elast - the elastic stiffness matrix
%   Q - logical array restricting nodes with Dirichlet boundary cond.
%   f - vector of external load
%   c_bar,sin_phi,shear,bulk,lame - material parameters at integration points
%
% Output data:
%   U_it, t_it - approximation of the solution
%   flag_N - logical array indicating possible inconvergence of the method
%
%--------------------------------------------------------------------------

%
% Initialization     
%
                       
  n_n=size(U_ini,2);
  n_int=length(WEIGHT);
  dim=size(U_ini,1);         % dimension (dim=2 or dim=3)
  n_strain=dim*(dim+1)/2;    % number of strain (stress) components
  V = zeros(dim,n_n);
  W = zeros(dim,n_n);
  F = zeros(dim,n_n) ;       % vector of internal forces
  E = zeros(n_strain,n_int); % strain tensors at integration points            
  U_it=U_ini;                % initial displacements
  t_it=t_ini;                % initial value of t
  flag_N=0;

%  
% Newton's solver (the semismooth Newton method)
%
      
  it=0;         % iteration number
  r=r_damp;     % stronger regularization of the tangent stiffness matrix
  while true         
      
     it=it+1;  
     
     % constitutive operator and its derivative
     E(:) = B*U_it(:) ;   % strain at integration points
     if dim==2            % solution of the constitutive problem
         [S,DS]=constitutive_problem_2D(E,c_bar,sin_phi,shear,bulk,lame);
     else
         [S,DS]=constitutive_problem_3D(E,c_bar,sin_phi,shear,bulk,lame);
     end     
   
     % vector of internal forces
     F(:) = B'*reshape(repmat(WEIGHT,n_strain,1).*S(1:n_strain,:), n_strain*n_int,1) ;         

     % stopping criterion and its test
     criterion = norm(F(Q)-t_it*f(Q))/norm(f(Q));         
     if  (criterion < tol)&&(it>1)
        fprintf('     Newton method converges:  '); 
        fprintf(' number of iteration=%d  ',it); 
        fprintf(' t=%e  ',t_it); 
        fprintf('\n'); 
        break
     end          
     
     % tangential stiffness matrix and its regularization
     AUX=reshape(1:n_strain*n_int,n_strain,n_int);
     iD=repmat(AUX,n_strain,1); 
     jD=kron(AUX,ones(n_strain,1));
     vD = repmat(WEIGHT,n_strain^2,1).*DS ;      
     D_p = sparse( iD(:),jD(:),vD(:), n_strain*n_int,n_strain*n_int ) ;   
     K_tangent = B'*D_p*B;   
     K_tangent=(K_tangent+K_tangent')/2;
     K_r = r_min*K_elast+(1-r_min)*K_tangent; 
        
     % Newton's increment
     W(Q) = K_r(Q,Q)\f(Q); 
     V(Q) = K_r(Q,Q)\(t_it*f(Q)-F(Q)); 
     d_t=-f(Q)'*V(Q)/(f(Q)'*W(Q));
     dU=V+d_t*W;       
     
     % damping parameter
     alpha=damping_ALG4(it_damp_max,U_it,dU,F,B,WEIGHT,c_bar,sin_phi,shear,bulk,lame);
     
     % recomputation with a regularized stiffness matrix when omega=0
     if alpha==0
       while alpha==0 
         r=r*10;
         K_r = r*K_elast+(1-r)*K_tangent;
         W(Q) = K_r(Q,Q)\f(Q); 
         V(Q) = K_r(Q,Q)\(t_it*f(Q)-F(Q)); 
         d_t=-f(Q)'*V(Q)/(f(Q)'*W(Q));
         dU=V+d_t*W;     
         alpha=damping_ALG4(it_damp_max,U_it,dU,F,B,WEIGHT,c_bar,sin_phi,shear,bulk,lame);
         disp('   Additional regularization of the stiffness matrix was necessary.')
       end % while
       r=r/20;
     end %if
     
     % next iteration
     U_it= U_it + alpha*dU ;  
     U_it=omega*U_it/(f(Q)'*U_it(Q));
     t_it=t_it+d_t;
      
     % test on number of iteration
     if  it == it_newt_max
        fprintf('     Newton solver does not converge: stopping criterion=%e  ',criterion)
        fprintf('\n'); 
        flag_N=1;
        break
     end     
     
  end % true     
  
end % function