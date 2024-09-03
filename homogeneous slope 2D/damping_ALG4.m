function alpha=damping_ALG4(it_damp_max,U_it,dU,F,B,WEIGHT,c_bar,sin_phi,shear,bulk,lame)

%--------------------------------------------------------------------------
% Computation of the damping coefficient for Newton's iteration
%--------------------------------------------------------------------------
 

  decrease=F(:)'*dU(:);
  if isnan(decrease)||(decrease>=0)
      alpha=0;
      return
  end
  
  dim=size(U_it,1);          % dimension (dim=2 or dim=3)
  n_strain=dim*(dim+1)/2;    % number of strain (stress) components
  alpha=1; alpha_min=0; alpha_max=1;
  it_damp=0;
  
  while it_damp<it_damp_max
    
    it_damp=it_damp+1;
    U_alpha = U_it + alpha*dU ;
    E_alpha = reshape( B*U_alpha(:) , n_strain,[] ) ;
    if dim==2            % solution of the constitutive problem
         S_alpha=constitutive_problem_2D(E_alpha,c_bar,sin_phi,shear,bulk,lame);
     else
         S_alpha=constitutive_problem_3D(E_alpha,c_bar,sin_phi,shear,bulk,lame);
    end    
    F_alpha = B'*reshape(repmat(WEIGHT,n_strain,1).*S_alpha(1:n_strain,:), [],1);
    decrease = F_alpha'*dU(:) ;
    if decrease<0
        if alpha==1
           break
        end
        alpha_min=alpha;
    else
        alpha_max=alpha;
    end
    alpha=(alpha_min+alpha_max)/2;
  end

end % function
     