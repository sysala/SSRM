function [c_bar,sin_phi]=reduction(c0,phi,psi,lambda,Davis_type)

%--------------------------------------------------------------------------
% This function reduces strength parameters depending on the factor lambda
% and the prescribed type of Davis' approach
%
% Input data:
%   c_0 - effective cohesion at integration points
%   phi - effective frction angle at integration points
%   psi - dilatancy angle at integration points
%   lambda - strength reduction factor 
%   Davis_type - choice of Davis' approach; available choice 'A','B','C'
%
% Output data:
%   c_bar = 2*c0_lambda.*cos(phi_lambda) 
%   sin_phi = sin(phi_lambda)
%
%--------------------------------------------------------------------------

  switch(Davis_type)
    case 'A'
        beta=cos(phi).*cos(psi)./(1-sin(phi).*sin(psi));
        c0_lambda=beta.*c0/lambda;
        phi_lambda=atan(beta.*tan(phi)/lambda);
        c_bar= 2*c0_lambda.*cos(phi_lambda);   
        sin_phi = sin(phi_lambda) ; 
    case 'B'
        c01=c0/lambda;
        phi1=atan(tan(phi)/lambda);
        psi1=atan(tan(psi)/lambda);
        beta=cos(phi1).*cos(psi1)./(1-sin(phi1).*sin(psi1));
        c0_lambda=beta.*c01;
        phi_lambda=atan(beta.*tan(phi1));
        c_bar= 2*c0_lambda.*cos(phi_lambda);   
        sin_phi = sin(phi_lambda) ; 
    case 'C'
        c01=c0/lambda;
        phi1=atan(tan(phi)/lambda);
        if phi1>psi
            beta=cos(phi1).*cos(psi)./(1-sin(phi1).*sin(psi));
        else
            beta=1;
        end
        c0_lambda=beta.*c01;
        phi_lambda=atan(beta.*tan(phi1));
        c_bar= 2*c0_lambda.*cos(phi_lambda);   
        sin_phi = sin(phi_lambda) ; 
    otherwise
        disp('Incorrect choice of the Davis aprroach.');
  end        
  
end