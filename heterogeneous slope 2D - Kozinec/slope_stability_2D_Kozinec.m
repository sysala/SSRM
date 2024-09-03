% =========================================================================
%
%  Heterogeneous slope in Doubrava-Kozinec and its stability
%
%  This program solves a 2D slope stability problem by the modified shear
%  strength reduction method suggested in (Sysala et al. 2021). It is
%  considered the Mohr-Coulomb yield criterion, 3 Davis approaches,
%  standard finite elements (P1, P2 or P4 elements) and uniform meshes 
%  with different densities. Gauss quadrature is used for numerical
%  integration. To find the safety factor of the SSR method, 2 continuation
%  techniques are available: the direct and the indirect techniques. The
%  indirect continuation is completed with 3 different Newton's algorithms.
% 
% ======================================================================
%

%
% The main input data 
%

  % elem_type - type of finite elements; available choices: 'P1', 'P2'          
  elem_type='P2';

  % Davis_type - choice of Davis' approach; available choices: 'A','B','C'
  Davis_type='B'; 
 
  % elastic material parameters (FoS should be independent of these param.)
  young = 16000;                     % Young's modulus
  poisson = 0.4;                     % Poisson's ratio
  shear = young/(2*(1+poisson)) ;    % shear modulus
  bulk = young/(3*(1-2*poisson)) ;   % bulk modulus    
  lame = bulk-2*shear/3 ;            % lame's coefficient (lambda)
  
  % 1. inelastic parameters for material #1
  c0_1 = 9 ;                         % cohesion
  phi_1 =26*pi/180;                  % friction angle  
  psi_1=0;                           % dilatation angle
  gamma_sat_1 = 20.3;                % specific weight - saturated
  gamma_unsat_1 = 20.7;              % specific weight - unsaturated  
  % 2. inelastic parameters for material #7
  c0_2 = 1 ;                         % cohesion
  phi_2 =45*pi/180;                  % friction angle  
  psi_2=0;                           % dilatation angle
  gamma_sat_2 = 20.5;                % specific weight - saturated
  gamma_unsat_2 = 20.6;              % specific weight - unsaturated  
  % 3. inelastic parameters for materials #4 and #6  
  c0_3 = 3 ;                         % cohesion
  phi_3 =13*pi/180;                  % friction angle  
  psi_3=0;                           % dilatation angle
  gamma_sat_3 = 20.0;                % specific weight - saturated
  gamma_unsat_3 = 20.5;              % specific weight - unsaturated  
  % 4. inelastic parameters for material #2
  c0_4 = 2 ;                         % cohesion
  phi_4 =33*pi/180;                  % friction angle 
  psi_4=0;                           % dilatation angle
  gamma_sat_4 = 19.0;                % specific weight - saturated
  gamma_unsat_4 = 20.5;              % specific weight - unsaturated  
  % 5. inelastic parameters for materials #3 and #5
  c0_5 = 5 ;                         % cohesion
  phi_5 =27*pi/180;                  % friction angle  
  psi_5=0;                           % dilatation angle
  gamma_sat_5 = 19.4;                % specific weight - saturated
  gamma_unsat_5 = 21.4;              % specific weight - unsaturated
     
%
% Data from the reference element
%
  
  % quadrature points and weights for volume integration
  [Xi, WF] = quadrature_volume_2D(elem_type);  
  % local basis functions and their derivatives 
  [HatP,DHatP1,DHatP2] = local_basis_volume_2D(elem_type, Xi);
  
%
%  Loading of the mesh imported from COMSOL
%
  
  % coordinates and Dirichlet nodes
  coord=load('coordinates3.txt')';
  
  % arrays over finite elements
  elem=load('elements3.txt')'+1; % vertices  
  mater=load('materials3.txt')';  % material type

  % extensions of the arrays elem and coord in the case of P2 or P4 elements 
  if strcmp(elem_type,'P2')==1
    [coord_mid, elem_mid]= midpoints_P2(coord,elem);
    coord=[coord, coord_mid];
    elem=[elem; elem_mid];
  end
  if strcmp(elem_type,'P4')==1
    [coord_mid, elem_mid]= midpoints_P4(coord,elem);
    coord=[coord, coord_mid];
    elem=[elem; elem_mid];
  end
  
  % logical array representing the Dirichlet boundary conditions
  x_max=max(coord(1,:));
  x_min=min(coord(1,:));
  n_n=size(coord,2);
  Q=false(2,n_n);
  Q(1,:) = (coord(1,:)>x_min+0.2)&(coord(2,:)>0.2)&(coord(1,:)<x_max-0.2) ;
  Q(2,:) = coord(2,:)>0.2 ;

  % number of nodes, elements and integration points + print
  n_unknown=length(coord(Q)); % number of unknowns
  n_e=size(elem,2);           % number of elements
  n_q=length(WF);             % number of quadratic points
  n_int = n_e*n_q ;           % total number of integrations points 
  % 
  fprintf('\n');   
  fprintf('Mesh data:'); 
  fprintf('  number of nodes =%d ',n_n);  
  fprintf('  number of unknowns =%d ',n_unknown);
  fprintf('  number of elements =%d ',n_e);
  fprintf('  number of integration points =%d ',n_int);
  fprintf('\n');   

%
% Material parameters at integration points
%

  % heterogeneous parameters
  [c0,phi,psi,gamma_sat,gamma_unsat]=heter_mater(mater,n_q,...
               c0_1,c0_2,c0_3,c0_4,c0_5,phi_1,phi_2,phi_3,phi_4,phi_5,...
               psi_1,psi_2,psi_3,psi_4,psi_5,gamma_sat_1,gamma_sat_2,...
               gamma_sat_3,gamma_sat_4,gamma_sat_5,gamma_unsat_1,...
               gamma_unsat_2,gamma_unsat_3,gamma_unsat_4,gamma_unsat_5);                         
   
  % specific weight at integration points depending on a given saturation
  % curve. This curve is prescribed within the function "gravity".
  gamma=gravity(gamma_sat,gamma_unsat,coord,elem,HatP);  

  % homogeneous (elastic) parameters 
  shear=shear*ones(1,n_int);
  bulk=bulk*ones(1,n_int);
  lame=lame*ones(1,n_int);

%  
% Assembling of the elastic stiffness matrix  
%
  [K_elast,B,WEIGHT]=elastic_stiffness_matrix_2D(elem,coord,...
                                              DHatP1,DHatP2,WF,shear,lame);  
 
%                                                
% Assembling of the vector of volume forces 
%

  % volume forces at integration points, size(f_V_int)=(2,n_int)
  f_V_int = [zeros(1,n_int);-gamma] ;  
  % vector of volume forces
  f_V=vector_volume_2D(elem,coord,f_V_int,HatP,WEIGHT);
  f=f_V;

%
% Computation of the factor of safety for the SSR method
%

  % input parameters for continuation
  lambda_init=0.7;     % initial lower bound of lambda
  d_lambda_init=0.1;   % initial increment of lambda 
  d_lambda_min=1e-4;   % minimal increment of lambda 
  step_max=100;        % maximal number of continuation steps
 
  % input parameters for Newton's solvers
  it_newt_max=50;               % number of Newton's iterations
  it_damp_max=10;               % number of iterations within line search
  tol=1e-6;                     % relative tolerance for Newton's solvers
  r_min=tol/100;                % basic regularization of the stiffness matrix
  r_damp=tol*100;               % regularization of the stiffness matrix if
                                % it is impossible to find descentdirection
  

  % direct continuation method - Algorithm 2
  fprintf('\n'); 
  fprintf('Direct continuation method - Algorithm 2'); 
  fprintf('\n'); 
  [U2,lambda_hist2,omega_hist2]=continuation_ALG2(...
               lambda_init,d_lambda_init,d_lambda_min,step_max,...
               it_newt_max,it_damp_max,tol,r_min,r_damp,...
               WEIGHT,B,K_elast,Q,f,c0,phi,psi,Davis_type,shear,bulk,lame);

  % indirect continuation method - Algorithm 3
  fprintf('\n'); 
  fprintf('Indirect continuation method - Algorithm 3'); 
  fprintf('\n'); 
  [U3,lambda_hist3,omega_hist3]=continuation_ALG3(...
               lambda_init,d_lambda_init,d_lambda_min,step_max,...
               it_newt_max,it_damp_max,tol,r_min,r_damp,...
               WEIGHT,B,K_elast,Q,f,c0,phi,psi,Davis_type,shear,bulk,lame);

%  
% Postprocessing - visualization of selected results for ALG2
%    

  % visualization of the mesh and the heterogeneity imported from COMSOL
  draw_mesh_Kozinec(coord,elem)
  draw_heterogeneity_Kozinec(coord,elem,mater)

  % visualization of total displacements and a deformed shape
  U_total = sqrt(U2(1,:).^2 + U2(2,:).^2);
  U_max=max(U_total);
  draw_quantity_Kozinec(coord,elem,10*U2/U_max,U_total)  
  
  % visualization of deviatoric strain
  IOTA=[1;1;0;1];                       % identity matrix in vector repres.
  VOL=IOTA*IOTA';                       % volumetric operator
  DEV=diag([1,1,1/2,1])-VOL/3;          % deviatoric operator
  E = zeros(3,n_int);                % strain tensors at integration points  
  E(:) = B*U2(:) ;                       % strain at integration points
  E_tr=[E;zeros(1,n_int)];              % trial strain
  dev_E=DEV*E_tr;                       % deviatoric part of E_tr
  norm_E=sqrt(max(0,sum(E_tr.*dev_E))); % norm of the deviatoric strain
  norm_E_node = transformation(norm_E,elem,WEIGHT); 
  draw_quantity_Kozinec(coord,elem,0*U2,norm_E_node)
    
  % visualization of the curve: omega -> lambda for Alg2
  figure
  hold on
  plot(omega_hist2, lambda_hist2,'-o');
  xlabel('control variable - \omega');
  ylabel('strength reduction factor - \lambda');
  hold off

%  
% Postprocessing - visualization of selected results for ALG3
%    

  % visualization of total displacements and a deformed shape
  U_total = sqrt(U3(1,:).^2 + U3(2,:).^2);
  U_max=max(U_total);
  draw_quantity_Kozinec(coord,elem,10*U3/U_max,U_total)    
  
  % visualization of deviatoric strain
  IOTA=[1;1;0;1];                       % identity matrix in vector repres.
  VOL=IOTA*IOTA';                       % volumetric operator
  DEV=diag([1,1,1/2,1])-VOL/3;          % deviatoric operator
  E = zeros(3,n_int);                % strain tensors at integration points  
  E(:) = B*U3(:) ;                      % strain at integration points
  E_tr=[E;zeros(1,n_int)];              % trial strain
  dev_E=DEV*E_tr;                       % deviatoric part of E_tr
  norm_E=sqrt(max(0,sum(E_tr.*dev_E))); % norm of the deviatoric strain
  norm_E_node = transformation(norm_E,elem,WEIGHT); 
  draw_quantity_Kozinec(coord,elem,0*U3,norm_E_node)

  % visualization of the curve: omega -> lambda for Alg3
  figure
  hold on
  plot(omega_hist3, lambda_hist3,'-o');
  xlabel('control variable - \omega');
  ylabel('strength reduction factor - \lambda');
  hold off

