function [U_it, flag_N]=...
                 newton(U_ini,tol,it_newt_max,it_damp_max,r_min,r_damp,...
                       WEIGHT,B,K_elast,Q,f,c_bar,sin_phi,shear,bulk,lame)
                               
%--------------------------------------------------------------------------
% The Newton method with damping and regularization for solving the system 
%                     find U:   F(U)=f
%
% Input data:
%   U_ini - initial choice of U
%   it_newt_max - number of Newton's iterations
%   it_damp_max - number of iterations within line search
%   tol - relative tolerance for Newton's solvers
%   eps - tolerance for numerical derivatives
%   r_min - basic regularization of the stiffness matrix
%   r_damp - regularization of the stiffness matrix if
%            it is impossible to find descent direction
%   WEIGHT - weight coefficients at integration points
%   B - the strain-displacement matrix
%   K_elast - the elastic stiffness matrix
%   Q - logical 2xn_n array restricting nodes with Dirichlet boundary cond.
%   f - vector of external load
%   shear,bulk,lame - arrays with elastic parameters
%   c_bar,sin_phi - arrays with inelastic parameters
%
% Output data:
%   U_it - the solution of the nonlinear system of equations
%   flag_N - logical array indicating possible inconvergence of the method
%
%--------------------------------------------------------------------------

%
% Initialization     
%
                       
  n_n=size(U_ini,2);         % number of nodes
  n_int=length(WEIGHT);      % number of integration points
  dim=size(U_ini,1);         % dimension (dim=2 or dim=3)
  dU = zeros(dim,n_n);       % Newton's increment 
  U_it=U_ini;                % initial displacements
  F = zeros(dim,n_n) ;       % vector of internal forces
  n_strain=dim*(dim+1)/2;    % number of strain (stress) components
  E = zeros(n_strain,n_int); % strain tensors at integration points               


  % constitutive operator and its derivative
  E(:) = B*U_it(:) ;   % strain at integration points
  if dim==2            % solution of the constitutive problem
         [S,DS]=constitutive_problem_2D(E,c_bar,sin_phi,shear,bulk,lame);
     else
         [S,DS]=constitutive_problem_3D(E,c_bar,sin_phi,shear,bulk,lame);
  end     
 
  % vector of internal forces
  F(:) = B'*reshape(repmat(WEIGHT,n_strain,1).*S(1:n_strain,:), n_strain*n_int,1) ;      

%  
% Newton's solver (the semismooth Newton method)
%
      
  it=0;         % iteration number
  flag_N=0;     % if flag_N==1, the Newton method does not converge
  r=r_damp;     % stronger regularization of the tangent stiffness matrix
  while true         
      
     it=it+1;       
         
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
     dU(Q) = K_r(Q,Q)\(f(Q)-F(Q)); 
     
     % damping parameter
     alpha=damping(it_damp_max,U_it,dU,F,f,B,WEIGHT,c_bar,sin_phi,shear,bulk,lame);
     
     % recomputation with a regularized stiffness matrix when omega=0
     if alpha==0 
       while alpha==0   
         r=r*10;
         K_r = r*K_elast+(1-r)*K_tangent;
         dU(Q) = K_r(Q,Q)\(f(Q)-F(Q)); 
         alpha=damping(it_damp_max,U_it,dU,F,f,B,WEIGHT,c_bar,sin_phi,shear,bulk,lame);
         disp('   Additional regularization of the stiffness matrix was necessary.')
       end
       r=r/20;
     end
     
     % next iteration
     U_it= U_it + alpha*dU ; 
     E(:) = B*U_it(:) ;   % strain at integration points    

     if dim==2            % solution of the constitutive problem
         [S,DS]=constitutive_problem_2D(E,c_bar,sin_phi,shear,bulk,lame);
     else
         [S,DS]=constitutive_problem_3D(E,c_bar,sin_phi,shear,bulk,lame);
     end     
    
     % vector of internal forces
     F(:) = B'*reshape(repmat(WEIGHT,n_strain,1).*S(1:n_strain,:), n_strain*n_int,1) ; 
 
     % stopping criterion and its test
     criterion = norm(F(Q)-f(Q))/norm(f(Q));         
     if  criterion < tol
        fprintf('     Newton method converges:  '); 
        fprintf(' number of iteration=%d  ',it); 
        fprintf(' stopping criterion=%e  ',criterion); 
        fprintf('\n'); 
        break
     end           
      
     % test on number of iteration
     if isnan(criterion)||(it == it_newt_max)
        flag_N=1;
        fprintf('     Newton solver does not converge: stopping criterion=%e  ',criterion)
        fprintf('\n'); 
        break
     end     
     
  end % true     
  
end % function