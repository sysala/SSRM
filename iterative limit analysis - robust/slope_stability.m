% =========================================================================
%
%  Slope stability
%
%  This program solves a 2D slope stability problem by the modified shear
%  strength reduction method suggested in (Sysala et al. 2021). It is
%  considered the Mohr-Coulomb yield criterion, 3 Davis approaches and
%  standard finite elements (either P1 or P2 elements). Numerical
%  quadrature is used according to a chosen type of elements. For P2 
%  elements, the 7-point Gauss quadrature is used. We apply an iterative 
%  limit analysis approach to find the safety factor. In particular, we use
%  a fixed limit analysis solver (as a black-box). Such a solver is robust
%  but it can be slower than the remaining investigated solvers. Its detail
%  description can be found within the limit analysis method.
% 
% ======================================================================
%

%
% The main input data 
%

  % elem_type - the type of finite elements; available choices: 'P1', 'P2'          
  elem_type='P2';

  % Davis_type - the choice of Davis' approach; available choice 'A','B','C'
  Davis_type='B';

  % geometrical parameters 
  x1 = 15 ;        % length of the body in front of the slope
  x3 = 15 ;        % length of the body behind the slope
  y1 = 10 ;        % hight of the body below the slope
  y2 = 10 ;        % height of the slope
  beta=pi/4;       % slope angle
  x2=y2/tan(beta); % length of the slope in x-direction
  
  % mesh data
  h  = 0.25;          % discretization parameter  

  % basic inelastic material parematers
  c0 = 6 ;                          % cohesion
  phi = 45*pi/180;                   % frictional angle  
  psi = 0; 
  
  % elastic material parameters (not necessary in limit analysis)
  young = 40000;                     % Young's modulus
  poisson = 0.3;                     % Poisson's ratio
  shear = young/(2*(1+poisson)) ;    % shear modulus
  bulk = young/(3*(1-2*poisson)) ;   % bulk modulus    
  lame = bulk-2*shear/3 ;            % lame's coefficient (lambda)
 
  % specific weight of a soil material
  gamma = 20 ;                      
     
%
% Data from the reference element
%
  
  % quadrature points and weights for volume integration
  [Xi, WF] = quadrature_volume_2D(elem_type);  
  % local basis functions and their derivatives 
  [HatP,DHatP1,DHatP2] = local_basis_volume_2D(elem_type, Xi);
  
%
% Creation of the finite element mesh
%  
    
  switch(elem_type)
    case 'P1'
        [coord,elem,ELEM_ED,EDGE_EL,Q]=mesh_P1(h,x1,x2,x3,y1,y2);
        fprintf('P1 elements: \n')
    case 'P2'
        [coord,elem,ELEM_ED,EDGE_EL,Q]=mesh_P2(h,x1,x2,x3,y1,y2);
        fprintf('P2 elements: \n')
    otherwise
        disp('bad choice of element type');
  end           
%   draw_mesh(coord,elem,x1+x2+x3,y1+y2) 

  % number of nodes, elements and integration points + print
  n_n=size(coord,2);          % number of nodes
  n_unknown=length(coord(Q)); % number of unknowns
  n_e=size(elem,2);           % number of elements
  n_ed=size(EDGE_EL,2);       % number of edges
  n_q=length(WF);             % number of quadratic points
  n_int = n_e*n_q ;           % total number of integrations points 
  % 
  fprintf('\n');   
  fprintf('The coarsest mesh:'); 
  fprintf('  number of nodes =%d ',n_n);  
  fprintf('  number of unknowns =%d ',n_unknown);
  fprintf('  number of elements =%d ',n_e);
  fprintf('  number of edges =%d ',n_ed);  
  fprintf('  number of integration points =%d ',n_int);
  fprintf('\n'); 

%
% Material parameters at integration points 
% (it unifies a treatment for homogeneous and heterogeneous structures)
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
  [K_elast,B,WEIGHT,iD,jD]=elastic_stiffness_matrix(elem,coord,...
                                              DHatP1,DHatP2,WF,shear,lame);  
 
%                                                
% Assembling of the vector of volume forces 
%

  % volume forces at integration points, size(f_V_int)=(2,n_int)
  f_V_int = [zeros(1,n_int);-gamma] ;  
  % vector of volume forces
  f_V=vector_volume(elem,coord,f_V_int,HatP,WEIGHT);
  f=f_V;

%
% Elastic solution and the corresponding initialization
%
  
  U_el = zeros(2,n_n);             % elastic displacement vector
  U_el(Q)=K_elast(Q,Q)\f(Q);
  alpha_el=f(Q)'*U_el(Q);          % work of external forces

  % initial increment of the regularization parameter alpha which is 
  % important for the limit analysis solver. The default factor 5 can be
  % changed if necessary
  d_alpha_ini=alpha_el/5;

%
% Computation of the SSRM factor of safety by iterative limit analysis
%

  % We solve the nonlinear equation ell(lambda*)=1, where the function ell
  % is nonincreasing and its value at lambda represents FoS for the limit
  % analysis method w.r.t. the strength parameters multiplied by lambda.
  % For any lambda, two functions are used: reduction - reduces (enlarges) 
  % the strength parameters by the factor lambda according to chosen Davis'
  % approach; limit_analysis - a limit analysis solver. In order to solve
  % the nonlinear equation, we find lower and upper bounds of lambda* at
  % first and then we apply the secant method.

  % other input parameters for the used limit analysis solver
  n_alpha_max=3;       % maximal of enalargements of alpha
  tol=1e-3;            % tolerance for the safety factors lamda^*, zeta^*
  step_max=100;        % maximal number of load steps 

  % limit analysis for lambda=1 (initialization) 
  lambda=1;           
  d_lambda=0.1;
  fprintf('\n'); 
  fprintf('Current value of lambda = %3.3f', lambda); 
  fprintf('\n'); 
  [c_bar,sin_phi]=reduction(c0,phi,psi,lambda,Davis_type);
  [U,zeta_hist,alpha_hist]=limit_analysis(...
                  d_alpha_ini,n_alpha_max,tol, step_max,...
                  WEIGHT,B,iD,jD,K_elast,U_el,Q,f,c_bar,sin_phi,...
                  shear,bulk,lame);
  ell=zeta_hist(end);

  % finding lower and upper bounds of lambda
  if ell>=1+tol
      lambda_min=1; % lower bound
      ell_max=ell;
      % we are looking for an upper bound of lambda*
      while true
          lambda=lambda+d_lambda;
          fprintf('\n'); 
          fprintf('Current value of lambda = %3.3f', lambda); 
          fprintf('\n'); 
          [c_bar,sin_phi]=reduction(c0,phi,psi,lambda,Davis_type);
          [U,zeta_hist,alpha_hist]=limit_analysis(...
                  d_alpha_ini,n_alpha_max,tol, step_max,...
                  WEIGHT,B,iD,jD,K_elast,U_el,Q,f,c_bar,sin_phi,...
                  shear,bulk,lame);
          ell=zeta_hist(end);
          if ell>=1+tol 
              lambda_min=lambda;
              ell_max=ell;
          elseif ell<=1-tol
              lambda_max=lambda;
              ell_min=ell;
              break
          else
              lambda_min=lambda;
              lambda_max=lambda;
              ell_max=ell;
              ell_min=ell;
              break
          end
      end
  elseif ell<=1-tol
      lambda_max=1; % upper bound
      ell_min=ell;
      % we are looking for a lower bound of lambda*
      while true
          lambda=lambda-d_lambda;
          fprintf('\n'); 
          fprintf('Current value of lambda = %3.3f', lambda); 
          fprintf('\n'); 
          [c_bar,sin_phi]=reduction(c0,phi,psi,lambda,Davis_type);
          [U,zeta_hist,alpha_hist]=limit_analysis(...
                  d_alpha_ini,n_alpha_max,tol, step_max,...
                  WEIGHT,B,iD,jD,K_elast,U_el,Q,f,c_bar,sin_phi,...
                  shear,bulk,lame);
          ell=zeta_hist(end);
          if ell<=1-tol
              lambda_max=lambda;
              ell_min=ell;
          elseif ell>=1+tol
              lambda_min=lambda;
              ell_max=ell;
              break
          else
              lambda_min=lambda;
              lambda_max=lambda;
              ell_max=ell;
              ell_min=ell;
              break
          end
      end
  else
      lambda_min=1;
      lambda_max=1;
      ell_max=1;
      ell_min=1;
  end

  % secant method for finding the strength reduction factor lambda*
  while lambda_max-lambda_min>tol
      lambda=lambda_min+(ell_max-1)*(lambda_max-lambda_min)/(ell_max-ell_min);
      fprintf('\n'); 
      fprintf('Current value of lambda = %3.3f', lambda); 
      fprintf('  lambda_min =%3.3f ',lambda_min);  
      fprintf('  lambda_max =%3.3f ',lambda_max);
      fprintf('  ell_min =%3.3f ',ell_min);
      fprintf('  ell_max =%3.3f ',ell_max);  
      fprintf('\n'); 
      [c_bar,sin_phi]=reduction(c0,phi,psi,lambda,Davis_type);
      [U,zeta_hist,alpha_hist]=limit_analysis(...
                  d_alpha_ini,n_alpha_max,tol, step_max,...
                  WEIGHT,B,iD,jD,K_elast,U_el,Q,f,c_bar,sin_phi,...
                  shear,bulk,lame);
      ell=zeta_hist(end);
      if ell>=1+tol 
         lambda_min=lambda;
         ell_max=ell;
      elseif ell<=1-tol
         lambda_max=lambda;
         ell_min=ell;
      else
         lambda_min=lambda;
         lambda_max=lambda;
         ell_max=ell;
         ell_min=ell;
         break
      end
  end  
  fprintf('\n'); 
  fprintf('Value of lambda* = %3.3f', lambda); 
  fprintf('  lambda_min =%3.3f ',lambda_min);  
  fprintf('  lambda_max =%3.3f ',lambda_max);
  fprintf('  ell_min =%3.3f ',ell_min);
  fprintf('  ell_max =%3.3f ',ell_max);  
  fprintf('\n'); 

%  
% Postprocessing - visualization of selected results
%      
  
  % visualization of total displacement rates + deformed shape
  U_total = sqrt(U(1,:).^2 + U(2,:).^2);
  draw_quantity(coord,elem,100*U,U_total,x1,x2,x3,y1,y2)  
  
  % visualization of deviatoric strain
  IOTA=[1;1;0;1];                       % identity matrix in vector repres.
  VOL=IOTA*IOTA';                       % volumetric operator
  DEV=diag([1,1,1/2,1])-VOL/3;          % deviatoric operator
  E = zeros(3,n_int);                % strain tensors at integration points  
  E(:) = B*U(:) ;                       % strain at integration points
  E_tr=[E;zeros(1,n_int)];              % trial strain
  dev_E=DEV*E_tr;                       % deviatoric part of E_tr
  norm_E=sqrt(max(0,sum(E_tr.*dev_E))); % norm of the deviatoric strain
  norm_E_node = transformation(norm_E,elem,WEIGHT); 
  draw_quantity(coord,elem,0*U,norm_E_node,x1,x2,x3,y1,y2)  
    
  % visualization of the curve alpha -> zeta for lambda*
  figure
  hold on
  plot(alpha_hist,zeta_hist);
  xlabel('work of external forces - \alpha'); ylabel('load parameter - \zeta');
  hold off

