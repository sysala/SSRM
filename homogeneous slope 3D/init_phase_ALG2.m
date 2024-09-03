function [U2,omega1,omega2,lambda2]=init_phase_ALG2...
                                 (lambda1,d_lambda1,d_lambda_min,...
                                 it_newt_max,it_damp_max,tol,eps,r_min,r_damp,...
                                 WEIGHT,B,K_elast,Q,f,...
                                 c0,phi,psi,Davis_type,shear,bulk,lame)

%--------------------------------------------------------------------------
% Initial phase for the continuation method given by ALG2. The aim is to
% complete (compute) the remaining input data of these algorithms.
%
% Input data:
%   lambda1 - initial value of the parameter lambda
%   d_lambda1 - initial increment of lambda
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
%   omega1  - value of the parameter omega corresponding to lambda1
%   lambda2 - the second value of the parameter lambda 
%   U2      - displacement field corresponding to lambda2
%   omega2  - value of the parameter omega corresponding to lambda2
%--------------------------------------------------------------------------

%  
% Initilization
%  
  n_n=size(Q,2);                % number of nodes
  dim=size(Q,1);                % dimension (dim=2 or dim=3)
  U_ini = zeros(dim,n_n);       % initial displacement field

%
% First step of the continuation method   
%
  
  fprintf('\n'); 
  fprintf(' Step=%d  ',1); 
  fprintf('\n'); 
  
  % Computation of U and omega for lambda1
  [U1,omega1,flag]=omega_ALG2(lambda1,U_ini,eps,1000,...
                          it_newt_max,it_damp_max,tol,r_min,r_damp,...
                          WEIGHT,B,K_elast,Q,f,...
                          c0,phi,psi,Davis_type,shear,bulk,lame);
    
  % test on convergence: 
  if flag==1   % the solver was not succesfull                      
     error('Initial choice of lambda seems to be too large.')
  end %if

  disp(['   lambda=', num2str(lambda1), ...
           ', d_lambda=', num2str(d_lambda1), ...
           ', omega=', num2str(omega1)]) ; 

%    
% Second step of the continuation method  
%

  fprintf('\n'); 
  fprintf(' Step=%d  ',2); 
  fprintf('\n'); 

  d_lambda=d_lambda1;       % initial increment of lambda
  while true                   
                   
    % update of the parameter lambda 
    lambda_it=lambda1+d_lambda;
                
    % Computation of U and omega for lambda_it
    [U_it,omega_it,flag]=omega_ALG2(lambda_it,U1,eps,d_lambda,...
                          it_newt_max,it_damp_max,tol,r_min,r_damp,...
                          WEIGHT,B,K_elast,Q,f,...
                          c0,phi,psi,Davis_type,shear,bulk,lame);
        
    % test on convergence: 
    if (flag==1)||(omega_it<=omega1)  % the solver was not succesfull                     
       % decrease of d_lambda 
       d_lambda=d_lambda/2;
    else        % the solver was successfull 
       U2=U_it;
       omega2=omega_it;
       lambda2=lambda_it;  
       break 
    end % if

    % stopping criteria
    if d_lambda<d_lambda_min
          error('It seems that FoS is equal to lambda_init.')
    end
  end % while    
  disp(['   lambda=', num2str(lambda2), ...
           ', d_lambda=', num2str(lambda2-lambda1), ...
           ', omega=', num2str(omega2),...
           ', d_omega=', num2str(omega2-omega1)]) ;

end % function