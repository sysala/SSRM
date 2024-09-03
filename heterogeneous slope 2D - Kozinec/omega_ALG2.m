function [U,omega,flag]=omega_ALG2(lambda_ini,U_ini,eps,d_lambda,...
                          it_newt_max,it_damp_max,tol,r_min,r_damp,...
                          WEIGHT,B,K_elast,Q,f,...
                          c0,phi,psi,Davis_type,shear,bulk,lame)

%--------------------------------------------------------------------------
% Computation of U and omega related to ALG2 for given lambda.
% In particular, F_lambda(U)=f, omega=-G'(lambda) 
%
% Input data:
%   lambda - current value of the parameter lambda
%   U_ini - initial guess of the displacemnt vector
%   eps - precision of the numerical derivative
%   d_lambda - current increment of the parameter lambda
%   it_newt_max - number of Newton's iterations
%   it_damp_max - number of iterations within line search
%   tol - relative tolerance for Newton's solvers
%   r_min - basic regularization of the stiffness matrix
%   r_damp - regularization of the stiffness matrix if
%            it is impossible to find descent direction
%   WEIGHT - weight coefficients at integration points
%   B - the strain-displacement matrix
%   K_elast - the elastic stiffness matrix
%   Q - logical array restricting nodes with Dirichlet boundary cond.
%   f - vector of external load
%   shear,bulk,lame - elastic material parameters at integration points
%   c_0, phi, psi - inelastic material parameters at integration points
%
% Output data:
%   U - computed displacement field
%   omega - computed value of the parameter omega
%   flag - logical value indicating Newton's converge/inconvergence
%--------------------------------------------------------------------------

% Initialization 
  omega=0;  % unrealistic value of omega in case of flag=1
  dim=size(Q,1);             % dimension (dim=2 or dim=3)
  n_strain=dim*(dim+1)/2;    % number of strain (stress) components

% Minimizer U of the function J_{lambda} and the value J_{lambda}(U)
  [c_bar,sin_phi]=reduction(c0,phi,psi,lambda_ini,Davis_type);
  [U,flag]=newton(U_ini,tol,it_newt_max,it_damp_max,r_min,r_damp,...
                       WEIGHT,B,K_elast,Q,f,c_bar,sin_phi,shear,bulk,lame);
  if flag==1
      return
  end
  E = reshape( B*U(:),n_strain,[] ) ;  % strain at integration points      
  if dim==2
     Psi= potential_2D(E,c_bar,sin_phi,shear,bulk,lame); 
  else
     Psi= potential_3D(E,c_bar,sin_phi,shear,bulk,lame); 
  end
  J=WEIGHT*Psi'-f(Q)'*U(Q);         

% Minimizer U_eps of the function J_{lambda-eps}
  beta=min(1,eps/d_lambda);
  U_beta=beta*U_ini+(1-beta)*U;
  [c_bar,sin_phi]=reduction(c0,phi,psi,lambda_ini-eps,Davis_type);
  [U_eps,flag]=newton(U_beta,tol,it_newt_max,it_damp_max,r_min,r_damp,...
                       WEIGHT,B,K_elast,Q,f,c_bar,sin_phi,shear,bulk,lame);
  if flag==1
      return
  end

% Value of the cost function J_{lambda-eps} for given U
  E = reshape( B*U_eps(:),n_strain,[]);  % strain at integration points 
  if dim==2
     Psi= potential_2D(E,c_bar,sin_phi,shear,bulk,lame); 
  else
     Psi= potential_3D(E,c_bar,sin_phi,shear,bulk,lame); 
  end
  J_eps=WEIGHT*Psi'-f(Q)'*U_eps(Q);   

% Computed value of the parameter omega
  omega=(J_eps-J)/eps;

end % function