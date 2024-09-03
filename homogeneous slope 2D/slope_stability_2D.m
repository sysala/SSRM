% =========================================================================
%
%  Homogeneous slope and its stability
%
%  This program solves a 2D slope stability problem by the modified shear
%  strength reduction method suggested in (Sysala et al. 2021). It is
%  considered the Mohr-Coulomb yield criterion, 3 Davis approaches,
%  standard finite elements (either P1 or P2 elements) and uniform meshes 
%  with different densities. For P2 elements, the 7-point Gauss quadrature 
%  is used. To find the safety factor of the SSR method, two continuation
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

  % geometrical parameters 
  x1 = 15 ;        % length of the body in front of the slope
  x3 = 15 ;        % length of the body behind the slope
  y1 = 10 ;        % hight of the body below the slope
  y2 = 10 ;        % height of the slope
  beta=pi/4;       % slope angle
  x2=y2/tan(beta); % length of the slope in x-direction
  
  % mesh data
  h  = 1;          % discretization parameter  

  % strength material parematers
  c0 = 6 ;                           % cohesion
  phi = 45*pi/180;                   % frictional angle  
  psi = 0;                           % dilatancy angle
  
  % elastic material parameters (FoS should be independent of these parameters)
  young = 40000;                     % Young's modulus
  poisson = 0.3;                     % Poisson's ratio
  shear = young/(2*(1+poisson)) ;    % shear modulus
  bulk = young/(3*(1-2*poisson)) ;   % bulk modulus    
  lame = bulk-2*shear/3 ;            % lame's coefficient (lambda)
 
  % specific weight of the material creating a slope
  gamma = 20 ;                      
     
%
% Data from the reference element
%
  
  % quadrature points and weights for volume integration
  [Xi, WF] = quadrature_volume_2D(elem_type);  
  % local basis functions and their derivatives 
  [HatP,DHatP1,DHatP2] = local_basis_volume_2D(elem_type, Xi);
  
%
% Creation of the uniform finite element mesh
%  
    
  switch(elem_type)
    case 'P1'
        [coord,elem,ELEM_ED,EDGE_EL,Q]=mesh_P1_2D(h,x1,x2,x3,y1,y2);
        fprintf('P1 elements: \n')
    case 'P2'
        [coord,elem,ELEM_ED,EDGE_EL,Q]=mesh_P2_2D(h,x1,x2,x3,y1,y2);
        fprintf('P2 elements: \n')
    otherwise
        disp('bad choice of element type');
  end      
  draw_mesh_2D(coord,elem,x1+x2+x3,y1+y2) 

  % number of nodes, elements and integration points + print
  n_n=size(coord,2);          % number of nodes
  n_unknown=length(coord(Q)); % number of unknowns
  n_e=size(elem,2);           % number of elements
  n_ed=size(EDGE_EL,2);       % number of edges
  n_q=length(WF);             % number of quadratic points
  n_int = n_e*n_q ;           % total number of integrations points 
  % 
  fprintf('\n');   
  fprintf('The mesh data:'); 
  fprintf('  number of nodes =%d ',n_n);  
  fprintf('  number of unknowns =%d ',n_unknown);
  fprintf('  number of elements =%d ',n_e);
  fprintf('  number of edges =%d ',n_ed);  
  fprintf('  number of integration points =%d ',n_int);
  fprintf('\n'); 

%
% Material parameters at integration points 
% (for a unified treatment of homogeneous and heterogeneous slopes)
%
  c0=c0*ones(1,n_int);
  phi=phi*ones(1,n_int);
  psi=psi*ones(1,n_int);
  shear=shear*ones(1,n_int);
  bulk=bulk*ones(1,n_int);
  lame=lame*ones(1,n_int);
  gamma=gamma*ones(1,n_int);

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
  lambda_init=0.9;     % initial lower bound of lambda
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

  % visualization of total displacements and a deformed shape
  U_total = sqrt(U2(1,:).^2 + U2(2,:).^2);
  U_max=max(U_total);
  draw_quantity_2D(coord,elem,0.2*y2*U2/U_max,U_total,x1,x2,x3,y1,y2)  
  
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
  draw_quantity_2D(coord,elem,0*U2,norm_E_node,x1,x2,x3,y1,y2)  
    
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
  draw_quantity_2D(coord,elem,0.2*y2*U3/U_max,U_total,x1,x2,x3,y1,y2)  
  
  % visualization of deviatoric strain
  IOTA=[1;1;0;1];                       % identity matrix in vector repres.
  VOL=IOTA*IOTA';                       % volumetric operator
  DEV=diag([1,1,1/2,1])-VOL/3;          % deviatoric operator
  E = zeros(3,n_int);                % strain tensors at integration points  
  E(:) = B*U3(:) ;                       % strain at integration points
  E_tr=[E;zeros(1,n_int)];              % trial strain
  dev_E=DEV*E_tr;                       % deviatoric part of E_tr
  norm_E=sqrt(max(0,sum(E_tr.*dev_E))); % norm of the deviatoric strain
  norm_E_node = transformation(norm_E,elem,WEIGHT); 
  draw_quantity_2D(coord,elem,0*U3,norm_E_node,x1,x2,x3,y1,y2)  

  % visualization of the curve: omega -> lambda for Alg3
  figure
  hold on
  plot(omega_hist3, lambda_hist3,'-o');
  xlabel('control variable - \omega');
  ylabel('strength reduction factor - \lambda');
  hold off

