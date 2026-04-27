%%  Homogeneous slope and its stability  (via SSR methods)
% ======================================================================
%  This program solves a 3D slope stability problem by the modified shear
%  strength reduction (SSR) method described in (Sysala et al., CAS 2025). 
%  The Mohr-Coulomb yield criterion, 3 Davis approaches (denoted by A, B, C),
%  standard finite elements (either P1 or P2 elements) and meshes with 
%  different densities are considered. For P2 elements, the 11-point 
%  Gauss quadrature is used. To find the safety factor of the SSR method, 
%  two continuation techniques are available: direct and indirect. 
%  A benchmark with a homogeneous slope from (Sysala et al., CAS 2025) is 
%  considered.
%
% ======================================================================

%% The main input data
% elem_type - type of finite elements; available choices: 'P2'
elem_type = 'P2';

% Davis_type - choice of Davis' approach; available choices: 'A','B','C'
Davis_type = 'B';

% Material parameters for each subdomain. In the following table, we
% specify in each column the following material parameters, respectively:
% [c0, phi, psi, young, poisson, gamma_sat, gamma_unsat], where
%    c0 ... Cohesion (c)
%    phi ... Friction angle (phi in degrees)
%    psi ... Dilatancy angle (psi in degrees)
%    young ... Young's modulus (E)
%    poisson ...  Poisson's ratio (nu)
%    gamma_sat ...   Specific weight - saturated (gamma_sat in kN/m^3)
%    gamma_unsat ... Specific weight - unsaturated (gamma_unsat in kN/m^3)
% If gamma_sat and gamma_unsat are not distinguished, use the same values 
% for these parameters. Each row of the table represents one subdomain. If 
% a homogeneous body is considered, only one row is prescribed.
mat_props = [6, 45, 0, 40000, 0.3, 20, 20]; 

% constant traction force prescribed in the normal direction
traction_normal = -50; 

% Geometrical parameters
x1 = 15;           % Length of the body in front of the slope
x3 = 15;           % Length of the body behind the slope
y1 = 10;           % Height of the body below the slope
y2 = 10;           % Height of the slope
beta = 45*pi/180;  % Slope angle
x2 = y2/tan(beta); % Length of the slope in the x-direction
z = 5;             % Thickness of the slope

% Mesh data
h = 1/1;         % Discretization parameter
N_h=round(1/h);

%% Data from the reference element
% quadrature points and weights for volume integration
[Xi, WF] = ASSEMBLY.quadrature_volume_3D(elem_type);
% local basis functions and their derivatives
[HatP, DHatP1, DHatP2, DHatP3] = ASSEMBLY.local_basis_volume_3D(elem_type, Xi);
% Quadrature points and weights for surface integration
[Xi_s, WF_s] = ASSEMBLY.quadrature_surface_3D(elem_type);
% Local basis functions and their derivatives for surface
[HatP_s, DHatP1_s, DHatP2_s] = ASSEMBLY.local_basis_surface_3D(elem_type, Xi_s);

%% Creation/loading of the finite element mesh

switch(elem_type)
    case 'P1'
        [coord, elem, surf, neumann, Q] = MESH.mesh_P1_3D_n(N_h, x1, x2, x3, y1, y2, z);
        fprintf('P1 elements: \n')
    case 'P2'
        [coord, elem, surf, neumann, Q] = MESH.mesh_P2_3D_n(N_h, x1, x2, x3, y1, y2, z);
        fprintf('P2 elements: \n')
    otherwise
        error('Bad choice of element type');
end

% number of nodes, elements and integration points + print
n_n = size(coord,2);          % number of nodes
n_unknown = length(coord(Q)); % number of unknowns
n_e = size(elem,2);           % number of elements
n_q = length(WF);             % number of quadratic points
n_int = n_e * n_q;            % total number of integration points
%
fprintf('\n The mesh data:');
fprintf('  number of nodes = %d ', n_n);
fprintf('  number of unknowns = %d ', n_unknown);
fprintf('  number of elements = %d ', n_e);
fprintf('  number of integration points = %d \n', n_int);

% The array material_identifier for a homogeneous body
material_identifier = zeros(1,n_e);

%% Material parameters at integration points
% Fields with prescribed material properties
fields = {'c0',      ... % Cohesion (c)
          'phi',     ... % Friction angle (phi in degrees)
          'psi',     ... % Dilatancy angle (psi in degrees)
          'young',   ... % Young's modulus (E)
          'poisson', ... % Poisson's ratio (nu)
          'gamma_sat', ... % Specific weight - saturated (gamma_sat in kN/m^3)
          'gamma_unsat'};  % Specific weight - unsaturated (gamma_unsat in kN/m^3)

% Convert properties to structured format.
materials = cellfun(@(x) cell2struct(num2cell(x), fields, 2), num2cell(mat_props, 2), 'UniformOutput', false);

% saturation - a prescribed logical array indicating integration points 
%              where the body is saturated. If gamma_sat and gamma_unsat 
%              are the same, set saturation=true(1,n_int). Otherwise,
%              this logical array is derived from a given phreatic surface.
saturation = true(1,n_int);

% Material parameters at integration points.
[c0, phi, psi, shear, bulk, lame, gamma] = ...
      ASSEMBLY.heterogenous_materials(material_identifier, saturation, n_q, materials);

%% Assembly of the elastic matrix 

[K_elast, B, WEIGHT] = ASSEMBLY.elastic_stiffness_matrix_3D(elem, coord, shear, bulk, DHatP1, DHatP2, DHatP3, WF);

%% Assembling of the vector of volume forces

% volume forces at integration points, size(f_V_int) = (3, n_int)
f_V_int = [zeros(1, n_int); -gamma; zeros(1, n_int)];
% vector of volume forces
f_V = ASSEMBLY.vector_volume_3D(elem, coord, f_V_int, HatP, WEIGHT);

%% Assembing of the vector of surface forces

% surface force in the normal direction at integration points
n_e_s=size(neumann,2); % number of surface elements
n_q_s=length(WF_s);   % number of quadrature points on a surface element
n_int_s=n_e_s*n_q_s ; % number of integration points on the surface
                      % (on the upper side of the body)
f_n_int=traction_normal*ones(1,n_int_s); % size(f_n_int)=(1,n_int_s)

% assembling of the vector of traction (surface) forces
f_t=ASSEMBLY.vector_traction_n_3D(neumann,coord,f_n_int,HatP_s,DHatP1_s,DHatP2_s,WF_s);

%% Input parameters for the continuation methods
lambda_init = 0.7;        % initial lower bound of lambda
d_lambda_init = 0.1;      % initial increment of lambda
d_lambda_min = 1e-3;      % minimal increment of lambda
d_lambda_diff_scaled_min = 0.001; % minimal rate of increment of lambda
omega_max_stop = 3000;    % maximal omega, then stop
step_max = 100;           % maximal number of continuation steps

%% Input parameters for Newton's solvers
it_newt_max = 100;               % number of Newton's iterations
it_damp_max = 10;                % number of iterations within line search
tol = 1e-4;                      % relative tolerance for Newton's solvers
r_min = 1e-4;                    % basic minimal regularization of the stiffness matrix

%% Defining linear solver
agmg_folder = "agmg"; % Check for AGMG in specified folder
solver_type = 'DFGMRES_AGMG'; % Type of solver: "DIRECT", "DFGMRES_ICHOL", "DFGMRES_AGMG"

linear_solver_tolerance = 1e-1;
linear_solver_maxit = 100;
deflation_basis_tolerance = 1e-3;
linear_solver_printing = 0;

[linear_system_solver] = LINEAR_SOLVERS.set_linear_solver(agmg_folder, solver_type, ...
    linear_solver_tolerance, linear_solver_maxit, deflation_basis_tolerance, linear_solver_printing, Q);


%% Constitutive problem and matrix builder
dim = 3;
n_strain = dim * (dim + 1) / 2;
constitutive_matrix_builder = CONSTITUTIVE_PROBLEM.CONSTITUTIVE(B, c0, phi, psi, Davis_type, shear, bulk, lame, WEIGHT, n_strain, n_int, dim);

%--------------------------------------------------------------------------
%% Computation of the factor of safety for the SSR method

direct_on = 1; % Whether direct continuation method should be used
indirect_on = 1; % Whether indirect continuation method should be used

if direct_on  % Direct continuation method
    fprintf('\n Direct continuation method\n');
    tic;
    [U2, lambda_hist2, omega_hist2, Umax_hist2] = CONTINUATION.SSR_direct_continuation(...
        lambda_init, d_lambda_init, d_lambda_min, d_lambda_diff_scaled_min, step_max, ...
        it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f_V+f_t, ...
        constitutive_matrix_builder, linear_system_solver.copy());
    time_run = toc;
    fprintf("Running_time = %f \n", time_run);
end
if indirect_on     % Indirect continuation method
    fprintf('\n Indirect continuation method\n');
    tic;
    [U3, lambda_hist3, omega_hist3, Umax_hist3] = CONTINUATION.SSR_indirect_continuation(...
        lambda_init, d_lambda_init, d_lambda_min, d_lambda_diff_scaled_min, step_max, ...
        omega_max_stop, it_newt_max, it_damp_max, tol, r_min, K_elast, Q, f_V+f_t, ...
        constitutive_matrix_builder, linear_system_solver.copy());
    time_run = toc;
    fprintf("Running_time = %f \n", time_run);
end

%% Postprocessing - visualization of selected results for direct continuation
if direct_on
    VIZ.plot_displacements_3D(U2, coord, elem, surf, 0.05);
    VIZ.plot_deviatoric_strain_3D(U2, coord, elem, surf, B);
    % Visualization of the curve: omega -> lambda for direct continuation.
    figure; hold on; box on; grid on;
    plot(omega_hist2, lambda_hist2, '-o');
    title('Direct continuation method', 'Interpreter', 'latex')
    xlabel('variable - $\xi$', 'Interpreter', 'latex');
    ylabel('strength reduction factor - $\lambda$', 'Interpreter', 'latex');
end

%% Postprocessing - visualization of selected results for indirect continuation
if indirect_on
    VIZ.plot_displacements_3D(U3, coord, elem, surf, 0.05);
    VIZ.plot_deviatoric_strain_3D(U3, coord, elem, surf, B);
    % Visualization of the curve: omega -> lambda for indirect continuation.
    figure; hold on; box on; grid on;
    plot(omega_hist3, lambda_hist3, '-o');
    title('Indirect continuation method', 'Interpreter', 'latex')
    xlabel('control variable - $\omega$', 'Interpreter', 'latex');
    ylabel('strength reduction factor - $\lambda$', 'Interpreter', 'latex');
end
