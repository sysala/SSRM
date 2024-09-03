function [U,lambda,flag]=bisection_method(U_ini,omega,lambda_ini,d_lambda,...
               it_newt_max,it_damp_max,tol,r_min,r_damp,...
               WEIGHT,B,K_elast,Q,f,c0,phi,psi,Davis_type,shear,bulk,lame)

%--------------------------------------------------------------------------
% Algorithm for solving the system of equations
%    find (U,lambda):  F_lambda(U)=f, f'*U=omega  for given omega.
% Its idea is to solve the following auxiliary system
%     find (U,t):  F_lambda(U)=t*f, f'*U=omega  for given omega and lambda.
% It enables us to define the function ell_omega(lambda)=t. We are looking 
% for such lambda for which ell_omega(lambda)=1. To this end, the bisection 
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
  fprintf('     lambda=%e  ',lambda);
  fprintf('\n'); 
  t_ini=1;
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
       % the solver was successfull; 
       if t<=1
           return % solution was found
       end
  end

%
% Finding an upper bound of lambda
%
    
    % initialization
    lambda_min=lambda; 
    d_lambda_min=10*tol;
    t_min=t; % t_min>1
    U_min=U;
    
    % while cyclus enlarging the parameter lambda 
    while t>1
        lambda_it=lambda+d_lambda/2;
        fprintf('     lambda=%e  ',lambda);
        fprintf('\n'); 
        % Davis' modifications of strength parameters using current lambda 
        [c_bar,sin_phi]=reduction(c0,phi,psi,lambda_it,Davis_type);
        % Newton's solver for extended system of non-linear equation
        [U_it,t_it,flag_N]=newton_ALG4(U,t,omega,it_newt_max,it_damp_max,...
                           tol,r_min,r_damp,WEIGHT,B,K_elast,Q,f,...
                           c_bar,sin_phi,shear,bulk,lame);
        % evaluation of the solver 
        if flag_N==1    % the solver was not successfull
          d_lambda=d_lambda/2;      
          if d_lambda<d_lambda_min
             warning('It is impossible to enlarge parameter lambda.')
             flag=1;
             return
          end
        else                                 
          % the solver was successfull; update of variables
          t=t_it;
          U=U_it;
          lambda=lambda_it;
          if t>1
              t_min=t; % t_min>1
              U_min=U;
              lambda_min=lambda;
          end % if   
          if lambda>2*lambda_ini
              warning('I seems that lambda* is not finite.')
              flag=1;
              return
          end
        end % if
    end % while
    t_max=t; % t_max<=1
    U_max=U;
    lambda_max=lambda;
    d_lambda=lambda_max-lambda_min;

%
% Finding lambda by the secant method
%
    
    while d_lambda>d_lambda_min
        % initialization
        lambda=(lambda_max+lambda_min)/2;
        fprintf('     lambda=%e  ',lambda);
        fprintf('\n'); 
        U=(U_max+U_min)/2;
        t=(t_max+t_min)/2;
        % Davis' modifications of strength parameters using current lambda 
        [c_bar,sin_phi]=reduction(c0,phi,psi,lambda,Davis_type);
        % Newton's solver for extended system of non-linear equation
        [U_it,t_it,flag_N]=newton_ALG4(U,t,omega,it_newt_max,it_damp_max,...
                           tol,r_min,r_damp,WEIGHT,B,K_elast,Q,f,...
                           c_bar,sin_phi,shear,bulk,lame);
        % evaluation of the solver 
        if flag_N==1  % the solver was not successfull;           
          flag=1;
          return
        else                                 
          % the solver was successfull; update of variables
          t=t_it;
          U=U_it;
          if t>1
              t_min=t; % t_min>1
              U_min=U;
              lambda_min=lambda;
          else
              t_max=t; % t_max<=1
              U_max=U;
              lambda_max=lambda;
          end % if 
          d_lambda=lambda_max-lambda_min;
        end % if
    end
    

end %function