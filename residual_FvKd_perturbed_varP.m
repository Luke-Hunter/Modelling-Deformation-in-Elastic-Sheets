%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to calculate finite differenced residuals of the displacement
% formulation Foeppl von Karman equations for solved for a circular domain, 
% subject to a radial tension mu with Poisson ratio nu
% for P that can vary with r
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#####################################################
# Edited to include altered nonlinear strain tensor #
#####################################################
function [R,J,sigma,d_sigma_dr]=residual_FvKd_perturbed_varP(guess,F_funct,Q_funct,nlc,mu,nu,xi)

% Nondimensionalisation means that E is 1
E=1;

% Set-up vectors of the unknowns
N=length(guess)/2;
h=1/(N-1);
W=guess(1:N);     % Deflection
U=guess(N+1:2*N); % Stress function
r=linspace(0,1,N); % Radial coordinate
F=F_funct(r);
Q=Q_funct(r);
% The strain tensor in polar coordinates
epsilon_fct = {...
  @(n,fd_fct)(fd_fct(1)(U,n,h)+nlc*0.5*fd_fct(1)(W,n,h).^2+xi*0.5*fd_fct(1)(W,n,h).^2),...
  @(n,fd_fct) 0.0,...
  @(n,fd_fct) (U(n)./r(n)+xi*0.5*(U(n)./r(n)).^2)...
 };

d_epsilon_fct_dr = {...
  @(n,fd_fct)(fd_fct(2)(U,n,h)+nlc*fd_fct(1)(W,n,h).*fd_fct(2)(W,n,h)+xi*fd_fct(2)(U,n,h).*fd_fct(1)(U,n,h)),...
  @(n,fd_fct) 0.0,...
  @(n,fd_fct) (r(n).*fd_fct(1)(U,n,h)-U(n))./(r(n).^2+xi*(r(n).*fd_fct(1)(U,n,h)-U(n))./(r(n).^2).*(U(n)./r(n)))...
  };

% The stress tensor in polar coordinates
sigma_fct = kirchhoff_st_venant_constitutive_cell(E,nu,epsilon_fct);
d_sigma_fct_dr = kirchhoff_st_venant_constitutive_cell(E,nu,d_epsilon_fct_dr);

% The expression at r=0 - NOT needed for now so left commented out
% The stain at r=0 - this needs a seperate expression due to the r=0 terms
% epsilon_atzero_fct = epsilon_fct;
% At zero we have 0/0 in epsilon_tt so use L'Hopitals rule
% epsilon_atzero_fct{3}  @(n,fd_fct) fd_fct(1)(u,n,h);
% sigma_atzero_fct = kirchhoff_st_venant_constitutive_cell(E,nu,...
%  epsilon_atzero_fct);

sigma   = {zeros(1,N),zeros(1,N),zeros(1,N)};
d_sigma_dr = {zeros(1,N),zeros(1,N),zeros(1,N)};

% Loop over tensorial components and precompute epsilon and sigma
for (i =1:3)
  % Fill in epsilon and sigma using backward differences at r=1 NOT needed for
  % now so left commented out
  % n=1;
  % sigma{i}(n) = sigma_atzero_fct{i}(n,@ForwardDifference);
  % Derivatives of strain and stress tensor
  % d_sigma_dr{i}(n) = d_sigma_atzero_fct_dr{i}(n,@ForwardDifference);
  n=2:N-1;
  % Fill in epsilon and sigma using central differences at r=2:N-1
  sigma{i}(n) = sigma_fct{i}(n,@CentralDifference);
  % Derivatives of strain and stress tensor
  d_sigma_dr{i}(n) = d_sigma_fct_dr{i}(n,@CentralDifference);
  % Fill in epsilon and sigma using backward differences at r=N
  n=N;
  sigma{i}(n) =sigma_fct{i}(n,@BackwardDifference);
  % Derivatives of strain and stress tensor
  d_sigma_dr{i}(n) =d_sigma_fct_dr{i}(n,@BackwardDifference);
end

% Deflection Equation
n=3:N-2;
R1=zeros(1,N);

% A shorthand wrapper around CentralDifference
cd=@(order)CentralDifference(order);

% Biharmonic operator and Forcing (linear bending equation)
R1(n)= cd(4)(W, n, h) +2*cd(3)(W, n, h)./r(n) -cd(2)(W, n, h)./r(n).^2 ...
     + cd(1)(W, n ,h)./r(n).^3  - F(n) ... 
      % The non linear term, with a homotopy parameter nlc
     -nlc*(cd(2)(W,n,h).*sigma{1}(n)+cd(1)(W,n,h).*sigma{3}(n)./r(n));

% In-plane Equation
n=2:N-1;
R2=zeros(1,N);

%  sigma_ij,j = 0
R2(n)= d_sigma_dr{1}(n) + (sigma{1}(n)-sigma{3}(n))./r(n)- Q(n);
      
 
% Boundary Conditions
% Odd derivatives zero at centre
R1(1)=ForwardDifference(1)(W,1,h); 
R1(2)=ForwardDifference(3)(W,1,h);

R1(N-1)=BackwardDifference(1)(W,N,h); % "Built--in" clamping
R1(N)=W(N);

R2(1)=U(1);  % U continuous across centre
% R2(1)=sigma(1,1)-sigma(3,1);
R2(N)=sigma{1}(N) - mu;  % Derivative at rim sets rim tension

R=[R1,R2];
 
% If we are outputting more than one argument
if(nargout>1)
 % Indexing function for the matrix
 % NB: would expect  index=@(r,c)(r-1)*N+c; but matlab (and octave) do the 
 % other way around. Consider using ind2sub as it has guards in it.
 index=@(r,c)(c-1)*2*N+r;

 % Kronecker delta
 delta=@(m,n) double(m==n);

% Now for physical quantities, defined for central points
 % Strains:
 d_epsilon_dw_fct = {...
  @(n,m,fd_fct,dfd_fct) nlc*dfd_fct(1)(n,m,h).*fd_fct(1)(W,n,h),...  
  @(n,m,fd_fct,dfd_fct) 0.0, ...
  @(n,m,fd_fct,dfd_fct) 0.0};

 d_epsilon_du_fct = {...
  @(n,m,fd_fct,dfd_fct) dfd_fct(1)(n,m,h)+xi*nlc*dfd_fct(1)(n,m,h).*fd_fct(1)(U,n,h),...
  @(n,m,fd_fct,dfd_fct) 0.0, ...
  @(n,m,fd_fct,dfd_fct) (delta(n,m)./r(n)).*(1+xi*(U(n)./r(n)))};

 % Strain derivatives:
 d2_epsilon_dr_dw_fct ={...
  @(n,m,fd_fct,dfd_fct) nlc*(dfd_fct(1)(n,m,h).*fd_fct(2)(W,n,h) + dfd_fct(2)(n,m,h).*fd_fct(1)(W,n,h)),...
  @(n,m,fd_fct,dfd_fct) 0.0 ,...
  @(n,m,fd_fct,dfd_fct) 0.0};

 d2_epsilon_dr_du_fct ={...
  @(n,m,fd_fct,dfd_fct) dfd_fct(2)(n,m,h)+xi*nlc*dfd_fct(2)(n,m,h).*fd_fct(1)(U,n,h)+xi*nlc*dfd_fct(1)(n,m,h).*fd_fct(2)(U,n,h),...
  @(n,m,fd_fct,dfd_fct) 0.0,... 
  @(n,m,fd_fct,dfd_fct)-delta(n,m)./(r(n).^2)+dfd_fct(1)(n,m,h)./r(n)+xi*((dfd_fct(1)(n,m,h).*U(n))./(r(n).^2))+xi*((delta(n,m).*fd_fct(1)(U,n,h))./(r(n).^2))-2*xi*((U(n).*delta(n,m))./(r(n).^3))};

 % Get d_sigma_dUi and d2_sigma_dr_dUi function handles
 d_sigma_dw_fct = kirchhoff_st_venant_constitutive_cell(E,nu,d_epsilon_dw_fct);
 d2_sigma_dr_dw_fct = kirchhoff_st_venant_constitutive_cell(E,nu,d2_epsilon_dr_dw_fct);

 d_sigma_du_fct = kirchhoff_st_venant_constitutive_cell(E,nu,d_epsilon_du_fct);
 d2_sigma_dr_du_fct = kirchhoff_st_venant_constitutive_cell(E,nu,d2_epsilon_dr_du_fct);
 
 %% Now fill in the Jacobian
 J=sparse(2*N,2*N);

 % Shorthand for cd
 dcd_df=@(order)DCentralDifferenceDf(order);

 %% Out of plane discrete equations 
 %% Do first oop equation - boundary condition #1
 n=1;
 % Unknown index
 m=1:3;
 J(n,m)= DForwardDifferenceDf(1)(1, m, h); % Forward difference about pt 1 (r=0)

 %% Do Second oop equation - boundary condition #2
 n=2;
 % Unknown index
 m=1:6;
 % Equation stored at n = 2 but the finite difference starts on point 1
 J(n,m)= DForwardDifferenceDf(3)(1, m , h); % Forward difference about pt 1 (r=0)

 % Loop oop equation index - this also coincides with the position in the domain
 n=3:N-2;
 % Now fill in five central diagonals
 for(l=-2:2)
  % Unknown index
  m=n+l;
  % Fill in the Jacobian (top right linear terms) for the central differences. NB Consider using sub2ind.
  J(index(n,m))= J(index(n,m))+dcd_df(4)(n, m, h) +2*dcd_df(3)(n, m, h)./r(n)...
     -dcd_df(2)(n, m, h)./r(n).^2 + dcd_df(1)(n, m, h)./r(n).^3;
  end

 % Non linear terms
 for(l=-1:1)
  % Unknown index
  m=n+l;
  % Fill in
  J(index(n,m))= J(index(n,m)) - nlc*(dcd_df(2)(n, m, h).*sigma{1}(n)...
     + dcd_df(1)(n, m, h).*sigma{3}(n)./r(n));
  J(index(n,m))= J(index(n,m)) - nlc*(cd(2)(W,n,h) .* d_sigma_dw_fct{1}(n,m,cd,dcd_df)...
     + cd(1)(W,n,h) .* d_sigma_dw_fct{3}(n, m, cd,dcd_df)./r(n));
  % Fill in the top right
  J(index(n,m+N))= J(index(n,m+N)) - nlc*(cd(2)(W,n,h) .* d_sigma_du_fct{1}(n,m,cd,dcd_df)...
     + cd(1)(W,n,h).*d_sigma_du_fct{3}(n,m,cd,dcd_df)./r(n));
 end  

 % Do second-to-last oop equation - wboundary condition #3
 % Equation stored at n = N-1 but the finite difference starts on point N
 n=N-1;
 % Unknown index
 m=N-3:N;
 J(n,m)= DBackwardDifferenceDf(1)(N, m ,h); % backward difference about pt N (r=1) 

 %% Do last oop equation - boundary condition #4
 n=N;
 % Unknown index
 m=N;
 J(n,m)=1.0;

 %% In plane discrete equations 
 %% Do first ip equation - u  boundary condition #1
 n=N+1;
 % Unknown index
 m=1;
 J(n,m+N) =1.0;
 
 % Range of ip equation index - this doesn't coincide with the position in the domain
 n=2:N-1;
 % Now fill in three central diagonals
 for(l=-1:1)
  % Unknown index
  m=n+l;
  % Fill in the Jacobian for the central differences. NB Consider using sub2ind.
  % First fill in the bottom right
  J(index(n+N,m+N))= J(index(n+N,m+N)) + d2_sigma_dr_du_fct{1}(n,m,cd,dcd_df)...
     + (d_sigma_du_fct{1}(n,m,cd,dcd_df) -  d_sigma_du_fct{3}(n,m,cd,dcd_df))./r(n);
  % Now fill in the bottom left
  J(index(n+N,m))= J(index(n+N,m)) + d2_sigma_dr_dw_fct{1}(n,m,cd,dcd_df)...
     + (d_sigma_dw_fct{1}(n,m,cd,dcd_df) -  d_sigma_dw_fct{3}(n,m,cd,dcd_df))./r(n);
 end  
 
 %% Do last ip equation - u boundary condition #2
 n=2*N;
 % Unknown index
 m=N-2:N;
 J(n,m+N)=d_sigma_du_fct{1}(N,m,@BackwardDifference,@DBackwardDifferenceDf);
 % Boundary condition depend on w too
 m=N-2:N;
 J(n,m)  =d_sigma_dw_fct{1}(N,m,@BackwardDifference,@DBackwardDifferenceDf);
end
% End of file
end 
