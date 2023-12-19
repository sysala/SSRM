function omega=damping(U_it,dU,F,B,WEIGHT,alpha,c_bar,sin_phi,shear,bulk,lame)

%--------------------------------------------------------------------------
% Computation of the damping coefficient for Newton's iteration
%--------------------------------------------------------------------------
 

  decrease=F(:)'*dU(:);
  if isnan(decrease)||(decrease>=0)
      omega=0;
      return
  end
  
  omega=1; omega_min=0; omega_max=1;
  it_damp=0;
  
  while it_damp<10
    
    it_damp=it_damp+1;
    U_omega = U_it + omega*dU ;
    E_omega = reshape( B*U_omega(:) , 3,[] ) ;
    [S_omega,DS]=constitutive_problem(E_omega,alpha,c_bar,sin_phi,shear,bulk,lame);
    F_omega = B'*reshape(repmat(WEIGHT,3,1).*S_omega(1:3,:), [],1);
    decrease = F_omega'*dU(:) ;
    if decrease<0
        if omega==1
           break
        end
        omega_min=omega;
    else
        omega_max=omega;
    end
    omega=(omega_min+omega_max)/2;
  end
  %omega=omega_min;

end % function
     