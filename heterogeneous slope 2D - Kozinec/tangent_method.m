function [U,lambda,flag]=tangent_method(U_ini,omega,lambda_ini,d_lambda,...
               it_newt_max,it_damp_max,tol,r_min,r_damp,...
               WEIGHT,B,K_elast,Q,f,c0,phi,psi,Davis_type,shear,bulk,lame)

%--------------------------------------------------------------------------
% Algorithm for solving the system of equations
%    find (U,lambda):  F_lambda(U)=f, f'*U=omega  for given omega.
% Its idea is to solve the following auxiliary system
%     find (U,t):  F_lambda(U)=t*f, f'*U=omega  for given omega and lambda.
% It enables us to define the function ell_omega(lambda)=t. We are looking 
% for such lambda for which ell_omega(lambda)=1. To this end, the tangent 
% method is used.
%
% Input data:
%   U_ini - initial choice of U
%   omega - given value of the parameter omega
%   lambda_ini - initial value of the parameter lambda
%   d_lambda - ititial increment of the parameter lambda
%   d_lambda_min - minimal increment of lambda
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
%   U,lambda - the solution of the investigated system
%   flag - a logical value indicating possible inconvergence
%--------------------------------------------------------------------------

%
% Initialization
%

  lambda=lambda_ini;   
  d_lambda_min=10*tol;
  fprintf('  lambda=%e  ',lambda);
  fprintf('\n'); 
  t_ini=1;
  eps=tol*10;                   % precision of numerical derivatives 
  flag=0;

  % Davis' modifications of strength parameters using initial lambda 
  [c_bar,sin_phi]=reduction(c0,phi,psi,lambda,Davis_type);
  
  % Newton's solver for the following system of nonlinear equation:
  %  given lambda,omega, find (U,t):  F_lambda(U)=tf, f'*u=omega  
  [U,t,flag_N]=newton_ALG4(U_ini,t_ini,omega,it_newt_max,it_damp_max,...
                           tol,r_min,r_damp,WEIGHT,B,K_elast,Q,f,...
                           c_bar,sin_phi,shear,bulk,lame);
  % test on convergence: 
  if flag_N==1 % the solver was not successfull       
       flag=1;
       return
  else                                 
       % the solver was successful; 
       if t>1
           lambda_min=lambda; 
           lambda_max=lambda+10*d_lambda;
       else
           lambda_min=lambda-d_lambda; 
           lambda_max=lambda;
       end
  end

%
% tangent method for solving the equation \ell_omega(lambda)=1
%
    
  while abs(d_lambda)>d_lambda_min

       %
       % Numerical derivative of \ell_omega at the current value of lambda
       %
        lambda_eps=lambda+eps;
        % Davis' modifications of strength parameters using current lambda 
        [c_bar,sin_phi]=reduction(c0,phi,psi,lambda_eps,Davis_type);
        % Newton's solver for extended system of non-linear equation
        [U_it,t_it,flag_N]=newton_ALG4(U,t,omega,it_newt_max,it_damp_max,...
                           tol,r_min,r_damp,WEIGHT,B,K_elast,Q,f,...
                           c_bar,sin_phi,shear,bulk,lame);
        % evaluation of the solver 
        if flag_N==1  % the solver was not successfull;           
          flag=1;
          return
        end
        D_ell=(t_it-t)/eps; % 

       %
       % Next Newton's iteration
       %
        if abs(D_ell)<1e-18 % almost singular Jacobian
            lambda=(lambda_min+lambda_max)/2;
            d_lambda=(lambda_max-lambda_min)/2;
        else                                 
          lambda_it=lambda+(1-t)/D_ell;
          lambda_it=min(lambda_max-d_lambda_min,max(lambda_min+d_lambda_min,lambda_it));
          d_lambda=lambda_it-lambda;
          lambda=lambda_it;
        end % if
        fprintf('  lambda=%e  ',lambda);
        fprintf('\n'); 

       %
       % Value of \ell_omega at the current value of lambda
       %
        % Davis' modifications of strength parameters using current lambda 
        [c_bar,sin_phi]=reduction(c0,phi,psi,lambda,Davis_type);
        % Newton's solver for extended system of non-linear equation
        [U_it,t_it,flag_N]=newton_ALG4(U,t,omega,it_newt_max,it_damp_max,...
                           tol,r_min,r_damp,WEIGHT,B,K_elast,Q,f,...
                           c_bar,sin_phi,shear,bulk,lame);
        % evaluation of the solver 
        if flag_N==1
          % the solver was not successfull; reduction of the increment
          flag=1;
          return
        else                                 
          % the solver was successfull; update of variables
          t=t_it;
          U=U_it;
          if t>1
              lambda_min=lambda;
          else
              lambda_max=lambda;
          end % if 
        end % if     

   end %while    

end %function