function [U,lambda_hist,omega_hist]=continuation_ALG3(...
                lambda_init,d_lambda_init,d_lambda_min,step_max,...
                it_newt_max,it_damp_max,tol,r_min,r_damp,...
                WEIGHT,B,K_elast,Q,f,c0,phi,psi,Davis_type,shear,bulk,lame)

%--------------------------------------------------------------------------
% Indirect continuation method estimating the safety factor lambda* and
% constructing the curve between the strength reduction parameter lambda
% and the parameter omega which represents the work of external forces.
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
%   lambda_hist, omega_hist - arrays for the construction of the curve 
%                             between lambda and omega
%--------------------------------------------------------------------------

%
% Initialization 
%

  % basic initial data
  lambda_hist=zeros(1,1000);    % history of lambda-values
  omega_hist=zeros(1,1000);     % history of omega-value

  % first two steps of the continuation method
  [U_old,U,omega_old,omega,lambda]=init_phase_ALG3...
                             (lambda_init,d_lambda_init,d_lambda_min,...
                             it_newt_max,it_damp_max,tol,r_min,r_damp,...
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
% While cyclus over the parameter omega
% 
  step=2;            % number of continuation steps
  n_omega=0;         % number of reductions of omega
  n_omega_max=5;     % maximal number of reductions of omega
  while true

      fprintf('\n'); 
      fprintf(' Step=%d  ',step+1); 
      fprintf('\n'); 
      
      % update of the parameter omega
      omega_it=omega+d_omega;

      % initial estimate of the displacement field
      U_ini=d_omega*(U-U_old)/(omega-omega_old)+U; 

      % solvers for the system F_lambda(U)=f, f'*U=omega_it

      % a) Newton's method
      [U_it,lambda_it,flag]=newton_ALG5(U_ini,omega_it,lambda,...
                                    it_newt_max,it_damp_max,tol,...
                                    r_min,r_damp,WEIGHT,B,K_elast,Q,f,...
                                    c0,phi,psi,Davis_type,shear,bulk,lame);
      
%       % b) bisection method w.r.t lambda
%       [U_it,lambda_it,flag]=bisection_method(U_ini,omega_it,lambda,d_lambda,...
%                it_newt_max,it_damp_max,tol,r_min,r_damp,...
%                WEIGHT,B,K_elast,Q,f,c0,phi,psi,Davis_type,shear,bulk,lame);

%       % c) tangent method w.r.t lambda
%       [U_it,lambda_it,flag]=tangent_method(U_ini,omega_it,lambda,d_lambda,...
%                it_newt_max,it_damp_max,tol,r_min,r_damp,...
%                WEIGHT,B,K_elast,Q,f,c0,phi,psi,Davis_type,shear,bulk,lame); 


      % evaluation of the solver and update
      if flag==1  % solver was not successful
         d_omega=d_omega/2;      
         n_omega=n_omega+1;
      else        % solver was successful
         % update 
         step=step+1;
         U_old=U;
         U=U_it;
         omega_old=omega;
         omega=omega_it;
         d_lambda=lambda_it-lambda;
         lambda=lambda_it;
         n_omega=0;
         % storage
         omega_hist(step)=omega;
         lambda_hist(step)=lambda;    
         % display
         disp(['   lambda=', num2str(lambda),...
                ', d_lambda=',num2str(d_lambda), ...
                ', omega=', num2str(omega),...
                ', d_omega=', num2str(d_omega)]) ;  
         % criterion for the change of d_omega
         if (lambda_hist(step)-lambda_hist(step-1))<...
                     1.1*(lambda_hist(step-1)-lambda_hist(step-2))
             d_omega=2*d_omega;
         end
      end % if flag
            
      % stopping criteria for the indirect continuation method
      if n_omega>=n_omega_max
          disp('It is impossible to increase omega.')
          break
      end      
      if step>=step_max
          disp('Maximal number of steps was achieved.')
          break
      end
      % stopping criteria
      if d_lambda<d_lambda_min
          disp('Minimal increment of lambda was achieved.')
          break
      end
            
  end %true  
  
  % clipping of the output arrays
  lambda_hist=lambda_hist(1:step);
  omega_hist=omega_hist(1:step);  
    
end % function