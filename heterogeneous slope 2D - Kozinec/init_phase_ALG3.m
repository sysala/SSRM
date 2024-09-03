function [U1,U2,omega1,omega2,lambda2]=init_phase_ALG3...
                                 (lambda1,d_lambda1,d_lambda_min,...
                                 it_newt_max,it_damp_max,tol,r_min,r_damp,...
                                 WEIGHT,B,K_elast,Q,f,...
                                 c0,phi,psi,Davis_type,shear,bulk,lame)

%--------------------------------------------------------------------------
% Initial phase for the continuation methods ALG2 and ALG3. The aim is to
% complete (compute) the remaining input data of these algorithms.
%
% Input data:
%   lambda1 - initial value of the parameter lambda
%   d_lambda_init - initial increment of lambda
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
%   U1 - displacement field corresponding to lambda1
%   U2 - displacement field corresponding to lambda2
%   omega1  - value of the parameter omega corresponding to lambda1
%   omega2  - value of the parameter omega corresponding to lambda2
%   lambda2 - the second value of the parameter lambda for which the Newton
%             method converges
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

  % Davis' modifications of strength parameters using initial lambda 
  [c_bar,sin_phi]=reduction(c0,phi,psi,lambda1,Davis_type);
               
  % Newton's solver for initial lambda     
  [U_it,flag_N]=newton(U_ini,tol,it_newt_max,it_damp_max,r_min,r_damp,...
                       WEIGHT,B,K_elast,Q,f,c_bar,sin_phi,shear,bulk,lame);
    
  % test on convergence: 
  if flag_N==1   % the solver was not succesfull                      
     error('Initial choice of lambda seems to be too large.')
  else           % the solver was successfull 
     U1=U_it;
     omega1=f(Q)'*U1(Q);     
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
                
    % update of strength parameters using lambda and Davis' approach
    [c_bar,sin_phi]=reduction(c0,phi,psi,lambda_it,Davis_type);
               
    % Newton's solver 
    [U_it,flag_N]=newton(U1,tol,it_newt_max,it_damp_max,r_min,r_damp,...
                       WEIGHT,B,K_elast,Q,f,c_bar,sin_phi,shear,bulk,lame);
    
    % test on convergence: 
    if flag_N==1  % the solver was not succesfull                     
       % decrease of d_lambda 
       d_lambda=d_lambda/2;
    else          % the solver was successfull 
       U2=U_it;
       omega2=f(Q)'*U2(Q);     
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
         ', omega=', num2str(omega1),...
         ', d_omega=', num2str(omega2-omega1)]) ;

end % function