% ************************************************************************
%
%  constitutive_problem   -  solution of the constitutive 
%                            Drucker-Prager problem 
%
% ************************************************************************

function   Psi= potential_3D(E,c_bar,sin_phi,shear,bulk,lame)
% =========================================================================
%
% The aim of this function is to construct a potential to the constitutive 
% operator for the Mohr-Coulomb yield criterion and associative flow rule
% at integration points 1,2,...,n_int. 
%
% Input data:
%  E     - current strain tensor, size(E)=(6,n_int)
%  shear,bulk,lame - parameters (constants)
%  c_bar,sin_phi - 1 x n_int arrays with heterogeneous parameters
%
% Output data:
%  Psi    - potential to E at integration points, size(Psi)=(1,n_int)
%
% =========================================================================                                   

%
% Trial strain and other auxiliary arrays:
%   E_tr   - trial strain tensors, size(E_tr)=(6,n_int)
%
  n_int=size(E,2); % number of integration points
  E_trial = E   ;  % strain representation
  IDENT = diag([1, 1, 1, 1/2, 1/2, 1/2]) ; % identity operator
  E_tr = IDENT*E_trial ;     % stress representation
 
%  
% Invariants of the trial strain tensors (at integration points)
%
  I1 = E_tr(1,:)+E_tr(2,:)+E_tr(3,:); % trace of E_trial
  I2 = E_tr(1,:).*E_tr(2,:)+E_tr(1,:).*E_tr(3,:)+E_tr(2,:).*E_tr(3,:)- ...
       E_tr(4,:).^2 - E_tr(5,:).^2 - E_tr(6,:).^2;
  I3 = E_tr(1,:).*E_tr(2,:).*E_tr(3,:) - E_tr(3,:).*E_tr(4,:).^2 - ...
       E_tr(2,:).*E_tr(6,:).^2 - E_tr(1,:).*E_tr(5,:).^2 + ...
       2*E_tr(4,:).*E_tr(5,:).*E_tr(6,:) ;
   
  Q = max(0,(1/9)*((I1.^2) - 3*I2));
  R = (1/54)*(-2*(I1.^3) + 9*(I1.*I2) - 27*I3);
  test1 = (Q==0);
  theta0 = zeros(1,n_int);
  theta0(~test1) = R(~test1)./sqrt(Q(~test1).^3);
  theta = acos(min(max(theta0,-1),1))/3; % Lode's angle

%  
% Ordered eigenvalues of the trial strain tensors
%
  eig_1 = -2*sqrt(Q).*cos(theta+2*pi/3) + I1/3 ;
  eig_2 = -2*sqrt(Q).*cos(theta-2*pi/3) + I1/3 ;
  eig_3 = -2*sqrt(Q).*cos(theta) + I1/3;

%
% Critical values defining decision criteria
%
  f_tr=2*shear.*((1+sin_phi).*eig_1-(1-sin_phi).*eig_3)+ ...
       2*(lame.*sin_phi).*I1-c_bar;             % elast-plast interface             
  gamma_sl=(eig_1-eig_2)./(1+sin_phi);          % smooth-left interface
  gamma_sr=(eig_2-eig_3)./(1-sin_phi);          % smooth-right interface
  gamma_la=(eig_1+eig_2-2*eig_3)./(3-sin_phi);  % left-apex interface
  gamma_ra=(2*eig_1-eig_2-eig_3)./(3+sin_phi);  % right-apex interface
  
%  
% Candidates on plastic multipliers
%
  denom_s=4*lame.*sin_phi.*sin_phi+4*shear.*(1+sin_phi.*sin_phi);   
  denom_l=4*lame.*sin_phi.*sin_phi+shear.*(1+sin_phi).*(1+sin_phi)+ ...
                                 2*shear.*(1-sin_phi).*(1-sin_phi);
  denom_r=4*lame.*sin_phi.*sin_phi+2*shear.*(1+sin_phi).*(1+sin_phi)+ ...
                                     shear.*(1-sin_phi).*(1-sin_phi);
  denom_a=4*bulk.*sin_phi.*sin_phi;  % denominators for each type of return
  
  lambda_s=f_tr./denom_s ;
  lambda_l=(shear.*((1+sin_phi).*(eig_1+eig_2)-2*(1-sin_phi).*eig_3)+ ...
            2*lame.*sin_phi.*I1-c_bar)./denom_l ;
  lambda_r=(shear.*(2*(1+sin_phi).*eig_1-(1-sin_phi).*(eig_2+eig_3))+ ...
            2*lame.*sin_phi.*I1-c_bar)./denom_r ;
  lambda_a=(2*bulk.*sin_phi.*I1-c_bar)./denom_a ;

%
% Determination of the potential Psi
%
  
  % initialization of unknowns variables
  Psi = zeros(1,n_int); 
  trace_E=eig_1+eig_2+eig_3;        % trace of E_trial
  
  % elastic response
  test_el=(f_tr<=0);
  Psi(test_el)=0.5*lame(test_el).*trace_E(test_el).^2+shear(test_el).*(eig_1(test_el).^2+...
               eig_2(test_el).^2+eig_3(test_el).^2);
 
  % return to the smooth portion of the yield surface
  test_s=(lambda_s<=min(gamma_sl,gamma_sr))&(~test_el);
  Psi(test_s)=0.5*lame(test_s).*trace_E(test_s).^2+shear(test_s).*(eig_1(test_s).^2+...
              eig_2(test_s).^2+eig_3(test_s).^2)-...
              0.5*denom_s(test_s).*(lambda_s(test_s).^2);
                
  % return to the left edge of the yield surface             
  test_l=(gamma_sl<gamma_sr)&(lambda_l>=gamma_sl)&(lambda_l<=gamma_la)&...
         (~(test_el|test_s));
  Psi(test_l)=0.5*lame(test_l).*trace_E(test_l).^2+shear(test_l).*(eig_3(test_l).^2+...
              0.5*(eig_1(test_l)+eig_2(test_l)).^2)-...
              0.5*denom_l(test_l).*(lambda_l(test_l).^2);
            
  % return to the right edge of the yield surface             
  test_r=(gamma_sl>gamma_sr)&(lambda_r>=gamma_sr)&(lambda_r<=gamma_ra)&...
         (~(test_el|test_s));
  Psi(test_r)=0.5*lame(test_r).*trace_E(test_r).^2+shear(test_r).*(eig_1(test_r).^2+...
              0.5*(eig_2(test_r)+eig_3(test_r)).^2)-...
              0.5*denom_r(test_r).*(lambda_r(test_r).^2);  
              
  % return to the apex of the yield surface
  test_a=~(test_el|test_s|test_l|test_r);
  Psi(test_a)=0.5*bulk(test_a).*trace_E(test_a).^2-...
              0.5*denom_a(test_a).*(lambda_a(test_a).^2);
          
  end
