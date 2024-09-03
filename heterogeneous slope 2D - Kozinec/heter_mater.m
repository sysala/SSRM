function [c0,phi,psi,gamma_sat,gamma_unsat]=heter_mater(mater,n_q,...
               c0_1,c0_2,c0_3,c0_4,c0_5,phi_1,phi_2,phi_3,phi_4,phi_5,...
               psi_1,psi_2,psi_3,psi_4,psi_5,gamma_sat_1,gamma_sat_2,...
               gamma_sat_3,gamma_sat_4,gamma_sat_5,gamma_unsat_1,...
               gamma_unsat_2,gamma_unsat_3,gamma_unsat_4,gamma_unsat_5)
                                 
%
% This function creates 1 x n_int arrays specifying material parameters at
% each integration point.
%
 
% a number of elements
  n_e=length(mater);
                    
% logical arrays specifying particular materials
  Q_mater1=(mater==1);
  Q_mater2=(mater==2);
  Q_mater3=(mater==3);
  Q_mater4=(mater==4);
  Q_mater5=(mater==5);
  Q_mater6=(mater==6);
  Q_mater7=(mater==7);
  
% the array c0
  c0=zeros(1,n_e);
  c0(Q_mater1)=c0_1;
  c0(Q_mater2)=c0_4;
  c0(Q_mater3)=c0_5;
  c0(Q_mater4)=c0_3;
  c0(Q_mater5)=c0_5;
  c0(Q_mater6)=c0_3;
  c0(Q_mater7)=c0_2;
  %
  c0=kron(c0,ones(1,n_q));
  
% the array phi
  phi=zeros(1,n_e);
  phi(Q_mater1)=phi_1;
  phi(Q_mater2)=phi_4;
  phi(Q_mater3)=phi_5;
  phi(Q_mater4)=phi_3;
  phi(Q_mater5)=phi_5;
  phi(Q_mater6)=phi_3;
  phi(Q_mater7)=phi_2;  
  %
  phi=kron(phi,ones(1,n_q));          

% the array psi
  psi=zeros(1,n_e);
  psi(Q_mater1)=psi_1;
  psi(Q_mater2)=psi_4;
  psi(Q_mater3)=psi_5;
  psi(Q_mater4)=psi_3;
  psi(Q_mater5)=psi_5;
  psi(Q_mater6)=psi_3;
  psi(Q_mater7)=psi_2;  
  %
  psi=kron(psi,ones(1,n_q));     

% saturated specific weight
  gamma_sat=zeros(1,n_e);
  gamma_sat(Q_mater1)=gamma_sat_1;
  gamma_sat(Q_mater2)=gamma_sat_4;
  gamma_sat(Q_mater3)=gamma_sat_5;
  gamma_sat(Q_mater4)=gamma_sat_3;
  gamma_sat(Q_mater5)=gamma_sat_5;
  gamma_sat(Q_mater6)=gamma_sat_3;
  gamma_sat(Q_mater7)=gamma_sat_2;
  %
  gamma_sat=kron(gamma_sat,ones(1,n_q));
  
% unsaturated specific weight
  gamma_unsat=zeros(1,n_e);
  gamma_unsat(Q_mater1)=gamma_unsat_1;
  gamma_unsat(Q_mater2)=gamma_unsat_4;
  gamma_unsat(Q_mater3)=gamma_unsat_5;
  gamma_unsat(Q_mater4)=gamma_unsat_3;
  gamma_unsat(Q_mater5)=gamma_unsat_5;
  gamma_unsat(Q_mater6)=gamma_unsat_3;
  gamma_unsat(Q_mater7)=gamma_unsat_2;
  %
  gamma_unsat=kron(gamma_unsat,ones(1,n_q));
                    
end