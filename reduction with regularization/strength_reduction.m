function [U,lambda_hist,lagr_hist]=strength_reduction...
                         (alpha,lambda_init,d_lambda_init,d_lambda_min,...
                          WEIGHT,B,K_elast,Q,f,...
                          c0,phi,psi,Davis_type,shear,bulk,lame)

%
% Initialization 
%
  lambda_hist=zeros(1,1000);    % history of lambda-values
  lagr_hist=zeros(1,1000);      % history of Lagr-values
  %
  lambda=lambda_init;           % initial bound of SSR factor
  d_lambda=d_lambda_init;       % initial increment of lambda
  lagr=-1000;                   % initial estimate of Lagrangian
  %
  n_n=size(Q,2);                % number of nodes
  U = zeros(2,n_n);             % the size of the displacement field
  %
  it_newt_max=25;               % number of Newton's iterations
  it_damp_max=10;               % number of iterations within line search
  tol1=1e-10;                   % relative tolerance 1 <
  tol3=1e-7;                    % relative tolerance 3
  flag=0;

%
% while cyclus over the parameter lambda
%    
  step=0;
  while d_lambda>d_lambda_min 
                   
    % update of the parameter lambda 
    lambda_it=lambda+d_lambda;
    if lambda_it<lambda_init+d_lambda_init
        disp('Initial choice of lambda is too large')
        break
    end
               
    % update of strength parameters using lambda and Davis' approach
    [c_bar,sin_phi]=reduction(c0,phi,psi,lambda_it,Davis_type);
               
    % Newton's solver 
    [U_it,lagr_it,it,criterion]=newton...
                            (tol1,it_newt_max,it_damp_max,...
                             U,alpha,lambda_it,lagr,WEIGHT,B,K_elast,Q,f,...
                             c_bar,sin_phi,shear,bulk,lame);
    
    % test on convergence: 
    if isnan(criterion)||((it==it_newt_max)&&(criterion>tol3))||(lagr_it < min(lagr,0))
                    % the solver was not succesfull 
                    % or the Lagrangian value decreased
       % decrease of d_lambda             
       d_lambda=d_lambda/2;
       % display of outputs
       disp(['           lambda=', num2str(lambda), ...
                       ', d_lambda=', num2str(d_lambda), ...
                       ', lagr=', num2str(lagr)]) ;
    else            % the solver was successfull and 
                    % the Lagrangian value increased
       % update of variables    
       if lagr_it<=lagr
           flag=1;
       end
       U=U_it;
       lambda=lambda_it;  
       lagr=lagr_it;
       step=step+1;
       lambda_hist(step)=lambda;
       lagr_hist(step)=lagr;
       if flag==1
           break
       end
       % display of outputs
       fprintf('\n');   
       disp(['   step=', num2str(step), ...
                       ', lambda=', num2str(lambda), ...
                       ', d_lambda=', num2str(d_lambda), ...
                       ', lagr=', num2str(lagr)]) ;  
     end %if
  end % while lambda     
  
  % clipping of the output arrays
  lambda_hist=lambda_hist(1:step);
  lagr_hist=lagr_hist(1:step);  
    
end % function