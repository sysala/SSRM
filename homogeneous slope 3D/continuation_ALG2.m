function [U,lambda_hist,omega_hist]=continuation_ALG2(...
                lambda_init,d_lambda_init,d_lambda_min,step_max,...
                it_newt_max,it_damp_max,tol,r_min,r_damp,...
                WEIGHT,B,K_elast,Q,f,c0,phi,psi,Davis_type,shear,bulk,lame)

%--------------------------------------------------------------------------
% Direct continuation method estimating the safety factor lambda* and
% constructing the curve between the strength reduction parameter lambda
% and the parameter omega=-G'(lambda).
%
% Input data:
%   lambda_init - initial value of the parameter lambda
%   d_lambda_init - initial increment of lambda
%   d_lambda_min - minimal increment of lambda
%   step_max - maximal number of continuation steps
%   it_newt_max - number of Newton's iterations
%   it_damp_max - number of iterations within line search
%   tol - relative tolerance for Newton's solvers
%   eps - precision of numerical derivatives
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
%   U - displacement field corresponding to lambda*
%   lambda_hist, omega_hist - arrays for construction of the curve between 
%                             lambda and omega
%
%--------------------------------------------------------------------------

%
% Initialization 
%

  % basic initial data
  lambda_hist=zeros(1,1000);    % history of lambda-values
  omega_hist=zeros(1,1000);     % history of omega-value
  eps=tol*10;                   % precision of numerical derivatives 

  % first two steps of the continuation method
  [U,omega_old,omega,lambda]=init_phase_ALG2...
                            (lambda_init,d_lambda_init,d_lambda_min,...
                             it_newt_max,it_damp_max,tol,eps,r_min,r_damp,...
                             WEIGHT,B,K_elast,Q,f,...
                             c0,phi,psi,Davis_type,shear,bulk,lame); 

  % storage of the computed values 
  omega_hist(1)=omega_old; 
  lambda_hist(1)=lambda_init;
  omega_hist(2)=omega; 
  lambda_hist(2)=lambda; 

  % other initial data
  d_omega=omega-omega_old;     % current increment of omega
  d_lambda=lambda-lambda_init; % current increment of lambda 
   
%
% While cyclus over the parameter lambda
%    
  step=2;     % current number of continuation steps
  while true 

    fprintf('\n'); 
    fprintf(' Step=%d  ',step+1); 
    fprintf('\n');   
                   
    % update of the parameter lambda 
    lambda_it=lambda+d_lambda;

    % Computation of U and omega for given lambda_it
    [U_it,omega_it,flag]=omega_ALG2(lambda_it,U,eps,d_lambda,...
                          it_newt_max,it_damp_max,tol,r_min,r_damp,...
                          WEIGHT,B,K_elast,Q,f,...
                          c0,phi,psi,Davis_type,shear,bulk,lame);
    
    % test on convergence: 
    d_omega_test=omega_it-omega_hist(step);
    if (flag==1)||(d_omega_test<0)
       % the solver was not succesfull - decrease of d_lambda 
       d_lambda=d_lambda/2;
    elseif d_omega_test>2*d_omega
       % too large increment of omega
    %   U=(U+U_it)/2;
       d_lambda=d_lambda/2;
    else  % the solver was successful 
       U=U_it;
       omega=omega_it;
       lambda=lambda_it;        
       step=step+1;
       lambda_hist(step)=lambda;
       omega_hist(step)=omega;        
       % display of outputs
       disp(['   lambda=', num2str(lambda),...
              ', d_lambda=',num2str(d_lambda), ...
              ', omega=', num2str(omega),...
              ', d_omega=', num2str(d_omega)]) ;  
       % update of d_lambda and d_omega
       if d_omega_test>1.5*d_omega
            d_lambda=d_lambda/2;           
       end       
       d_omega=d_omega_test;
  
    end % if flag

    % stopping criteria
    if d_lambda<d_lambda_min
          disp('Minimal increment of lambda was achieved.')
          break
    end
    if step>=step_max
          disp('Maximal number of steps was achieved.')
          break
    end
 
  end % while lambda     
  
  % clipping of the output arrays
  lambda_hist=lambda_hist(1:step);
  omega_hist=omega_hist(1:step);  
    
end % function