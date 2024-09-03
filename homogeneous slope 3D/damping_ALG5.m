function alpha=damping_ALG5(it_damp_max,U_it,lambda_it,d_U,d_l,f,criterion,...
                            B,WEIGHT,Q,c0,phi,psi,Davis_type,shear,bulk,lame)

%--------------------------------------------------------------------------
% Computation of the damping coefficient for ALG5
%--------------------------------------------------------------------------
   
  if isnan(d_l)
      alpha=0;
      return
  end

  dim=size(U_it,1);          % dimension (dim=2 or dim=3)
  n_strain=dim*(dim+1)/2;    % number of strain (stress) components
  n_n=size(U_it,2);
  F_alpha = zeros(dim,n_n) ;   % vector of internal forces

  alpha=1; 
  it_damp=0;
  
  while true
    
    it_damp=it_damp+1;
    U_alpha = U_it + alpha*d_U ;
    E_alpha = reshape( B*U_alpha(:) , n_strain,[] ) ;
    lambda_alpha=lambda_it+alpha*d_l;
    [c_bar,sin_phi]=reduction(c0,phi,psi,lambda_alpha,Davis_type);
    if dim==2            % solution of the constitutive problem
         S_alpha=constitutive_problem_2D(E_alpha,c_bar,sin_phi,shear,bulk,lame);
     else
         S_alpha=constitutive_problem_3D(E_alpha,c_bar,sin_phi,shear,bulk,lame);
    end     
    F_alpha(:) = B'*reshape(repmat(WEIGHT,n_strain,1).*S_alpha(1:n_strain,:), [],1);
    crit_alpha = norm(F_alpha(Q)-f(Q));

    % stopping criteria
    if crit_alpha>=criterion
        alpha=alpha/2;
    else
        break  % the damping parameter was found
    end

    if it_damp>=it_damp_max
        alpha=0;
        break  % the damping parameter was not found
    end

  end % while

end % function
     