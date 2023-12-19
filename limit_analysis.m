function [U,zeta_hist,alpha_hist]=limit_analysis...
                   (d_alpha_ini,n_alpha_max,d_zeta_min, step_max,...
                    WEIGHT,B,iD,jD,K_elast,U_elast,Q,f,c_bar,sin_phi,...
                    shear,bulk,lame)

%--------------------------------------------------------------------------
% The aim of this function is to solve the limit analysis problem by a
% regularization method. The regularization parameter alpha represents the
% work of external forces and enables us to control loading process.
% For any value of alpha, the following nonlinear system is solved:
%          find U, zeta:   F(alpha*U)=zeta*f, f'*U=1,
% where alpha*U is the displacement vector and zeta is the load parameter whose
% limit value zeta^* is unknown. It holds that zeta->zeta^* as
% alpha->+infty. The parameter alpha is enlarged adaptively. We start with 
% the constant increment d_alpha_ini. Increment d_alpha is enlarged 
% if the corresponding increment d_zeta of the load parameter is less than 
% given value d_zeta_min. We prescribe a fix number of enlargements of the
% increments (n_alpha_max) and a maximal number of the load steps (step_max) in
% order to terminate the algorithm. Let us note that d_alpha may decrease 
% if it isn't possible to find a solution of the system above by the Newton 
% method.
%
% Input data:
%   d_alpha_ini - initial increment of alpha
%   n_alpha_max - maximal number of enlargments of the increment d_alpha
%   d_zeta_min - critical value of the increment d_zeta
%   step_max - maximal number of load steps
%   WEIGHT - weight coefficients of integration points, size(WEIGHT)=(1,n_int)
%   B - the strain-displacement matrix, size(B)=(3*n_int,2*n_n)
%   D_elast=sparse(iD,jD,vD_elast) - elastic stress-strain matrix,
%                                    size(D_elast)=(3*n_int,3*n_int)
%   K_elast=B'*D_elast*B - the elastic stiffness matrix, size(K)=(2*n_n, 2*n_n)
%   U_elast - elastic solution, size(U_elast)=(2,n_n)
%   Q - logical 2xn_n array restricting nodes with Dirichlet boundary cond.
%   f - vector of external load
%   c_bar,sin_phi,shear,bulk,lame - material parameters at integration points
%
% Output data:
%   U - displacement field for a maximal achieved value of alpha
%   zeta_hist - history of the load parameters (1xstep array)
%   alpha_hist - history of the coontrol parameters (1xstep array)
%
%--------------------------------------------------------------------------
  
%
% Initialization 
%
  alpha_hist=zeros(1,step_max); % history of the parameter alpha
  zeta_hist=zeros(1,step_max);  % history of the load factor
  alpha=0;
  alpha_old=0;
  d_alpha=d_alpha_ini;          % increment of alpha
  zeta=0;                       % current load parameter
  U=U_elast/(f(Q)'*U_elast(Q)); % current displacement vector
  U_old=U; U_ini=U;
  n_alpha=0;

%
% loading process
%
  step=1;
  while true
      
      % Update of the control parameter
      alpha_it=alpha+d_alpha;
      
      % Newton's solver and its initialization by extrapolation using
      % solutions from previous steps
      if step>1
          U_ini=d_alpha*(U-U_old)/(alpha-alpha_old)+U;
      end
      [U_it,zeta_it,it,criterion]=newton(U_ini,alpha_it,...
                                 WEIGHT,B,iD,jD,K_elast,Q,f,...
                                 c_bar,sin_phi,shear,bulk,lame);
           
      % test on convergence: 
      if isnan(criterion)||((it==25)&&(criterion>1e-8)) 
          % the solver was not succesfull; reduction of the increment
          d_alpha=d_alpha/2;      
          n_alpha=n_alpha-1;
      else                                 
          % the solver was succefull; update of variables
          d_zeta=zeta_it-zeta;
          if (abs(d_zeta)<d_zeta_min)&&(criterion<1e-10)            
             d_alpha=2*d_alpha;
             n_alpha=n_alpha+1;
          end
          step=step+1;
          U_old=U;
          U=U_it;
          alpha_old=alpha;
          alpha=alpha_it;
          zeta=zeta_it;
          alpha_hist(step)=alpha;
          zeta_hist(step)=zeta;     
          disp(['  step=', num2str(step), ', alpha=', ...
                num2str(alpha), ', d_alpha=', num2str(d_alpha), ...
                ', zeta=', num2str(zeta),...
                ', d_zeta=',num2str(d_zeta)]) ;  
      end
      
      % stopping criteria for the loading process
      if d_alpha<=d_alpha_ini/10
          disp('Too small load increments.')
          break
      end      
      if step>=step_max
          disp('Maximal number of steps was achieved.')
          break
      end
      if n_alpha>=n_alpha_max
          disp(['Limit analysis solver was successfull:  ', ...
                'ell(lambda)=', num2str(zeta)])
          break
      end
            
  end %true  
  
  % clipping of the output arrays
  zeta_hist=zeta_hist(1:step);
  alpha_hist=alpha_hist(1:step);
  
end % function