function f_t = vector_traction_n_2D(ELEM_s,COORD,f_n_int,HatP_s,DHatP1_s,WF_s)

% =========================================================================
%
% Assembling of the vector of traction forces acting on the upper side of
% the 2D body
%
%    output: 
%      f_t - vector of traction forces, size(f_V)=(2,n_n), where n_n is 
%            the number of nodes
%
%    input data:
%      ELEM_s   - to indicate nodes belonging to each surface element 
%                 size(ELEM_s)=(n_p_s,n_e_s) where n_e_s is a number of 
%                 surface elements n_p_s is a number of the nodes within one
%                 surface element                           
%      COORD    - coordinates of the nodes, size(COORD)=(2,n_n)
%      f_n_int  - values of normal traction forces at integration points
%                 size(f_n_int)=(1,n_int_s), where n_int_s=n_e_s*n_q_s is a 
%                 number of surface integration points and n_q_s is a number  
%                 of quadrature points on the surface
%      HatP_s   - values of the surface basis functions at quadrature points
%      DHatP1_s - xi_1-derivatives of the surface basis functions at q. p.
%                 size(HatP_s)=size(DHatP1_s)=(n_p_s,n_q_s)
%      WF_s     - weight factors at surface quadrature points, 
%                 size(WF)=(1,n_q_s)
%
% =========================================================================

%
% Auxilliary notation
%

  n_n=size(COORD,2);    % number of nodes including midpoints
  n_e_s=size(ELEM_s,2); % number of surface elements
  n_p_s=size(ELEM_s,1); % number of nodes within one surface element
  n_q_s=length(WF_s);   % number of quadrature points on a surface element
  n_int_s=n_e_s*n_q_s ; % number of integration points on the surface
                        % (on the upper side of the body)

%
% Jacobians and their determinants at surface integration points
%                         

  % extension of the input arrays HatP_s,DHatP1_s by replication
  % size(HatPhi_s)=size(DHatPhi_s)=(n_p_s,n_int_s)
  DHatPhi_s=repmat(DHatP1_s,1,n_e_s);
   HatPhi_s=repmat(HatP_s,1,n_e_s)  ;
  
  % coordinates of nodes defining each surface element
  % size(COORDs1)=size(COORDs2)=(n_p_s,n_e_s)
  COORDs1=reshape(COORD(1,ELEM_s(:)),n_p_s,n_e_s);
  COORDs2=reshape(COORD(2,ELEM_s(:)),n_p_s,n_e_s);
  
  % coordinates of nodes around each surface integration point
  % size(COORDint1)=size(COORDint2)=(n_p_s,n_int_s)
  COORDint1=kron(COORDs1,ones(1,n_q_s)); 
  COORDint2=kron(COORDs2,ones(1,n_q_s)); 
  
  % derivatives of the isoparametric function at integration points
  % size(Dchi1)=size(Dchi2)=(1,n_int_s)
  Dchi1=sum(COORDint1.*DHatPhi_s);
  Dchi2=sum(COORDint2.*DHatPhi_s);

  % transformation ratio between real and reference surfaces at each int.p.
  % size(ratio)=(1,n_int_s)
  ratio=sqrt(Dchi1.^2+Dchi2.^2);

  % components of outward normal vectors at integration points
  % size(normal1)=size(normal2)=(1,n_int_s)
  normal1=Dchi2./ratio;
  normal2=-Dchi1./ratio;  
  
  % weight coefficients: size(WEIGHT_s)=(1,n_int_s)
  WEIGHT_s = ratio.*repmat(WF_s,1,n_e_s);
  
%
% Assembling of the vector of traction forces, size(f_t)=(2,n_n)
%

  % auxilliary values at surface integration points, 
  % size(vF1)=size(vF2)=(n_p_s,n_int_s)   
  vF1 = HatPhi_s.*(ones(n_p_s,1)*(WEIGHT_s.*normal1.*f_n_int));    
  vF2 = HatPhi_s.*(ones(n_p_s,1)*(WEIGHT_s.*normal2.*f_n_int));    
  % row and column indices, size(iF)=size(jF)=(n_p_s,n_int_s)   
  iF = ones(n_p_s,n_int_s);
  jF = kron(ELEM_s,ones(1,n_q_s));
  % the asssembling by using the sparse command - values v for duplicate
  % doubles i,j are automatically added together
  f_t = [ sparse(iF(:), jF(:), vF1(:), 1, n_n);
          sparse(iF(:), jF(:), vF2(:), 1, n_n) ];  

 end