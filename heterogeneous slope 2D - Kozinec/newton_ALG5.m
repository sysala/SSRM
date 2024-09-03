function [U_it,lambda_it,flag_N]=newton_ALG5(U_ini,omega,lambda_ini,...
                                    it_newt_max,it_damp_max,tol,...
                                    r_min,r_damp,WEIGHT,B,K_elast,Q,f,...
                                    c0,phi,psi,Davis_type,shear,bulk,lame)

%--------------------------------------------------------------------------
% The Newton method with regularization solving the system
%    find (U,lambda):  F_lambda(U)=f, f'*U=omega  for given omega
%
% Input data:
%   U_ini - initial choice of U
%   lambda_ini - initial choice of lambda
%   omega - prescribed value of the parameter omega
%   it_newt_max - number of Newton's iterations
%   tol - relative tolerance for Newton's solvers
%   r_min - basic regularization of the stiffness matrix
%   r_damp - regularization of the stiffness matrix if
%            it is impossible to find descent direction
%   WEIGHT - weight coefficients of integration points
%   B - the strain-displacement matrix
%   K_elast - the elastic stiffness matrix
%   Q - logical array restricting nodes with Dirichlet boundary cond.
%   f - vector of external load
%   c0,phi,psi,shear,bulk,lame - material parameters at integration points
%   Davis_type - a type of Davis' approximation
%
% Output data:
%   U_it, t_it - approximation of the solution
%   flag_N - logical array indicating possible inconvergence of the method
%
%--------------------------------------------------------------------------

%
% Initialization     
%

  n_n=size(U_ini,2);         % number of nodes
  n_int=length(WEIGHT);      % number of integration points
  dim=size(U_ini,1);         % dimension (dim=2 or dim=3)
  V = zeros(dim,n_n);
  W = zeros(dim,n_n);
  U_it=U_ini;                % initial displacements
  F = zeros(dim,n_n) ;       % vector of internal forces
  F_eps = zeros(dim,n_n) ;   % vector of internal forces
  n_strain=dim*(dim+1)/2;    % number of strain (stress) components
  E = zeros(n_strain,n_int); % strain tensors at integration points  
  lambda_it=lambda_ini;      % initial value of lambda
  eps=tol/1000;              % precision of numerical derivatives 
  norm_f=norm(f(Q));
  flag_N=0;
 
%  
% Newton's solver (the semismooth Newton method with regularization)
%

  it=0;         % iteration number
  r=r_damp;     % stronger regularization of the tangent stiffness matrix
  while true         
      
     it=it+1;  

     % Davis' modifications of strength parameters using current lambda 
     [c_bar,sin_phi]=reduction(c0,phi,psi,lambda_it,Davis_type);
     
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
     criterion = norm(F(Q)-f(Q));         
     if  criterion < tol*norm_f
        fprintf('     Newton method converges:  '); 
        fprintf(' number of iteration=%d  ',it); 
        fprintf(' stopping criterion=%e  ',criterion/norm_f); 
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
 
     % Numerical derivative of F w.r.t. lambda    
     lambda_eps=lambda_it+eps;
     [c_bar,sin_phi]=reduction(c0,phi,psi,lambda_eps,Davis_type);
     E(:) = B*U_it(:) ;   % strain at integration points
     if dim==2            % solution of the constitutive problem
         S=constitutive_problem_2D(E,c_bar,sin_phi,shear,bulk,lame);
     else
         S=constitutive_problem_3D(E,c_bar,sin_phi,shear,bulk,lame);
     end         
     F_eps(:) = B'*reshape(repmat(WEIGHT,n_strain,1).*S(1:n_strain,:), n_strain*n_int,1);
     G=(F_eps-F)/eps;
     
     % Newton's increment
     W(Q) = -K_r(Q,Q)\G(Q); 
     V(Q) = K_r(Q,Q)\(f(Q)-F(Q)); 
     d_l=-f(Q)'*V(Q)/(f(Q)'*W(Q));
     d_U=V+d_l*W;  

     % damping parameter
     alpha=damping_ALG5(it_damp_max,U_it,lambda_it,d_U,d_l,f,criterion,...
                        B,WEIGHT,Q,c0,phi,psi,Davis_type,shear,bulk,lame);
                           
     % recomputation with a regularized stiffness matrix when omega=0
     if alpha==0
        fprintf('     Newton solver does not converge: stopping criterion=%e  ',criterion)
        fprintf('\n'); 
        flag_N=1;
        break
%        while alpha==0 
%          r=r*10;
%          K_r = r*K_elast+(1-r)*K_tangent;
%          W(Q) = -K_r(Q,Q)\G(Q); 
%          V(Q) = K_r(Q,Q)\(f(Q)-F(Q)); 
%          d_l=-f(Q)'*V(Q)/(f(Q)'*W(Q));
%          d_U=V+d_l*W;  
%          alpha=damping_ALG5(it_damp_max,U_it,lambda_it,d_U,d_l,f,criterion,...
%                             B,WEIGHT,Q,c0,phi,psi,Davis_type,shear,bulk,lame);
%          disp('   Additional regularization of the stiffness matrix was necessary.')
%        end % while
%        r=r/20;
     end %if

     % update of U and lambda
     U_it=U_it+alpha*d_U;
     U_it=omega*U_it/(f(Q)'*U_it(Q));
     lambda_it=lambda_it+alpha*d_l;     
      
     % unsuccessful stopping criterion
     if  it == it_newt_max
        fprintf('     Newton solver does not converge: stopping criterion=%e  ',criterion)
        fprintf('\n'); 
        flag_N=1;
        break
     end   
     
  end % while      

end %function