####################################################################################################
## Function to calculate finite differenced residuals of the displacement                         ##
## formulation Foppl von Karman equations (WITH BENDING TERMS) for solved for a circular domain,  ##
## subject to a radial tension mu with Poisson ratio nu and thickness parameter xi = h/R          ##
####################################################################################################

function [R] = residual_bending_axisymmetric(guess,P,nu,xi)
  
# non-dimensionalisation
E=1;
xi2=xi^2;
xi4=xi^4;
  
# set up vectors of unknowns
N=length (guess)/2; # number of points
h=1/(N-1);         # distance between points
W=guess(1:N);      # out-of-plane displacement
U=guess(N+1:2*N);  # in-plane displacement
r=linspace(0,1,N); # mesh of points

# strain tensor in polar coordinates (xi terms from nonlinear strain contribution)
epsilon_fct = {...
  @(n,fd_fct)(fd_fct(1)(U,n,h)+0.5*fd_fct(1)(W,n,h).^2+xi2*0.5*fd_fct(1)(W,n,h).^2),...
  @(n,fd_fct) 0.0,...
  @(n,fd_fct) (U (n)./r(n)+xi2*0.5*(U (n)./r(n)).^2)...
 };
 
d_epsilon_fct_dr = {...
  @(n,fd_fct) (fd_fct(2)(U,n,h)+fd_fct(1)(W,n,h).*fd_fct(2)(W,n,h)+xi2*fd_fct(2)(U,n,h).*fd_fct(1)(U,n,h)),...
  @(n,fd_fct) 0.0,...
  @(n,fd_fct) (r(n).*fd_fct(1)(U,n,h)-U(n))./(r(n).^2+xi2*(r(n).*fd_fct(1)(U,n,h)-U(n))./(r(n).^2).*(U (n)./r(n)))...
  };

# stress tensor from strain using constitutive law
sigma_fct = kirchhoff_st_venant_constitutive_cell(E,nu,epsilon_fct);
d_sigma_fct_dr = kirchhoff_st_venant_constitutive_cell(E,nu,d_epsilon_fct_dr);

# christoffel symbols determined from strain tensor
G_rrr_fct = @(n,fd_fct) fd_fct(2)(U,n,h) + fd_fct(1)(W,n,h).*fd_fct(2)(W,n,h) + xi2.*fd_fct(1)(U,n,h).*fd_fct(2)(U,n,h);
G_rtt_fct = @(n,fd_fct) -(1./r(n)).*fd_fct(1)(U,n,h) - U (n)./(r (n).^2) - xi2.*((U (n)./(r (n).^2)).*fd_fct(1)(U,n,h) - (U (n).^2)./(r (n).^3));

d_G_rrr_fct_dr = @(n,fd_fct) fd_fct(3)(U,n,h) + fd_fct(2)(W,n,h).^2 + fd_fct(1)(W,n,h).*fd_fct(3)(W,n,h) + xi2.*fd_fct(2)(U,n,h).^2 + xi2.*fd_fct(1)(U,n,h).*fd_fct(3)(U,n,h);
d_G_rtt_fct_dr = @(n,fd_fct) -((1./r(n)).*fd_fct(2)(U,n,h) - 2*U (n)./(r (n).^3) + (xi2./(r (n).^2)).*(fd_fct(1)(U,n,h).^2 + U (n).*fd_fct(2)(U,n,h) - (4.*U (n).*fd_fct(1)(U,n,h))./r(n) + 3.*(U (n).^2)./(r (n).^2)));

# determinant of metric tensor as function of displacements
# (used to determine normal)
g = @(n,fd_fct) (1+xi2.*fd_fct(1)(U,n,h)).^2+xi2.*(fd_fct(1)(W,n,h)).^2;
d_g_dr = @(n,fd_fct) 2.*xi2.*(fd_fct(2)(U,n,h).*(1 + xi2.*fd_fct(1)(U,n,h))) + 2.*xi2.*(fd_fct(1)(W,n,h).*fd_fct(2)(W,n,h));
d2_g_dr2 = @(n,fd_fct) 2.*xi2.*(fd_fct(3)(U,n,h).*(1 + xi2.*fd_fct(1)(U,n,h))) + 2.*xi2.*(fd_fct(1)(W,n,h).*fd_fct(3)(W,n,h)) + 2.*xi2.*(fd_fct(2)(W,n,h)).^2 + 2.*xi4.*(fd_fct(2)(U,n,h)).^2;

# normal vector as function of displacements
# (used to determine curvature tensor and bending moments)
n_r = @(n,fd_fct) (-xi.*fd_fct(1)(W,n,h))./sqrt(g(n,fd_fct));
n_z = @(n,fd_fct) (1+xi2.*fd_fct(1)(U,n,h))./sqrt(g(n,fd_fct));

d_n_r_dr = @(n,fd_fct) (-xi.*fd_fct(2)(W,n,h))./sqrt(g(n,fd_fct)) + (0.5*xi.*fd_fct(1)(W,n,h).*d_g_dr(n,fd_fct))./sqrt(g(n,fd_fct).^3);
d_n_z_dr = @(n,fd_fct) (xi2.*fd_fct(2)(U,n,h))./sqrt(g(n,fd_fct)) + (0.5*(1+xi2.*fd_fct(1)(U,n,h)).*d_g_dr(n,fd_fct))./sqrt(g(n,fd_fct).^3);

d2_n_r_dr2 = @(n,fd_fct) (-xi.*fd_fct(3)(W,n,h))./sqrt(g(n,fd_fct)) + (xi.*fd_fct(2)(W,n,h).*d_g_dr(n,fd_fct))./sqrt(g(n,fd_fct).^3) + (0.5*xi.*fd_fct(1)(W,n,h).*d2_g_dr2(n,fd_fct)./sqrt(g(n,fd_fct).^3)) - (0.75*xi.*fd_fct(1)(W,n,h).*(d_g_dr(n,fd_fct).^2)./sqrt(g(n,fd_fct).^5));
d2_n_z_dr2 = @(n,fd_fct) (xi2.*fd_fct(3)(U,n,h))./sqrt(g(n,fd_fct)) + (0.5*(1+xi2.*fd_fct(1)(U,n,h)).*d2_g_dr2(n,fd_fct)./sqrt(g(n,fd_fct).^3)) - (0.75*(1+xi2.*fd_fct(1)(U,n,h)).*(d_g_dr(n,fd_fct).^2)./sqrt(g(n,fd_fct).^5));

# curvature tensor in polar coordinates
b_fct = {...
  @(n,fd_fct) xi.*n_r (n,fd_fct).*fd_fct(2)(U,n,h) + n_z (n,fd_fct).*fd_fct(2)(W,n,h),...
  @(n,fd_fct) 0.0,...
  @(n,fd_fct) (-1./r(n)).*(xi.*n_r (n,fd_fct).*fd_fct(1)(U,n,h) + n_z (n,fd_fct).*fd_fct(1)(W,n,h))...
 };

d_b_fct_dr = {...
  @(n,fd_fct) xi.*d_n_r_dr (n,fd_fct).*fd_fct(2)(U,n,h) + d_n_z_dr (n,fd_fct).*fd_fct(2)(W,n,h) + xi.*n_r (n,fd_fct).*fd_fct(3)(U,n,h) + n_z (n,fd_fct).*fd_fct(3)(W,n,h),...
  @(n,fd_fct) 0.0,...
  @(n,fd_fct) (1./(r(n).^2)).*(xi.*n_r (n,fd_fct).*fd_fct(1)(U,n,h) + n_z (n,fd_fct).*fd_fct(1)(W,n,h)) - (1./r(n)).*(xi.*d_n_r_dr (n,fd_fct).*fd_fct(1)(U,n,h) + d_n_z_dr (n,fd_fct).*fd_fct(1)(W,n,h) + xi.*n_r (n,fd_fct).*fd_fct(2)(U,n,h) + n_z (n,fd_fct).*fd_fct(2)(W,n,h)) ...
 };
 
d2_b_fct_dr2 = {...
  @(n,fd_fct) xi.*d2_n_r_dr2 (n,fd_fct).*fd_fct(2)(U,n,h) + 2.*xi.*d_n_r_dr (n,fd_fct).*fd_fct(3)(U,n,h) + d2_n_z_dr2 (n,fd_fct).*fd_fct(2)(W,n,h) + 2.*d_n_z_dr (n,fd_fct).*fd_fct(3)(W,n,h) + xi.*n_r (n,fd_fct).*fd_fct(4)(U,n,h) + n_z (n,fd_fct).*fd_fct(4)(W,n,h),...
  @(n,fd_fct) 0.0,...
  @(n,fd_fct) (2./(r(n).^2)).*(xi.*d_n_r_dr (n,fd_fct).*fd_fct(1)(U,n,h) + d_n_z_dr (n,fd_fct).*fd_fct(1)(W,n,h) + xi.*n_r (n,fd_fct).*fd_fct(2)(U,n,h) + n_z (n,fd_fct).*fd_fct(2)(W,n,h)) - (2./(r(n).^3)).*(xi.*n_r (n,fd_fct).*fd_fct(1)(U,n,h) + n_z (n,fd_fct).*fd_fct(1)(W,n,h)) - (1./r(n)).*(xi.*d2_n_r_dr2 (n,fd_fct).*fd_fct(1)(U,n,h) + 2.*xi.*d_n_r_dr (n,fd_fct).*fd_fct(2)(U,n,h) + xi.*n_r (n,fd_fct).*fd_fct(3)(U,n,h) + d2_n_z_dr2 (n,fd_fct).*fd_fct(1)(W,n,h) + 2.*d_n_z_dr (n,fd_fct).*fd_fct(2)(W,n,h) + n_z (n,fd_fct).*fd_fct(3)(W,n,h))...
 };
 
# The bending moments are determined from the curvature tensor
M_r_fct = kirchhoff_st_venant_constitutive_cell(E,nu,b_fct);
M_z_fct = kirchhoff_st_venant_constitutive_cell(E,nu,b_fct);

d_M_r_fct_dr = kirchhoff_st_venant_constitutive_cell(E,nu,d_b_fct_dr);
d_M_z_fct_dr = kirchhoff_st_venant_constitutive_cell(E,nu,d_b_fct_dr);

d2_M_r_fct_dr2 = kirchhoff_st_venant_constitutive_cell(E,nu,d2_b_fct_dr2);
d2_M_z_fct_dr2 = kirchhoff_st_venant_constitutive_cell(E,nu,d2_b_fct_dr2);

# Loop over tensorial components
for (i =1:3)
  
  n=3:N-2;
  
  # Fill in sigma using central differences at r=2:N-1
  sigma{i}(n) = sigma_fct{i}(n,@CentralDifference);
  # Derivatives of stress tensor
  d_sigma_dr{i}(n) = d_sigma_fct_dr{i}(n,@CentralDifference);
  
  # Fill in Christoffel symbols using central differences at r=2:N-1
  G_rrr(n) = G_rrr_fct(n,@CentralDifference);
  G_rtt(n) = G_rtt_fct(n,@CentralDifference);

  # Derivatives of Christoffel symbol tensors
  d_G_rrr_dr(n) = d_G_rrr_fct_dr(n,@CentralDifference);
  d_G_rtt_dr(n) = d_G_rtt_fct_dr(n,@CentralDifference);
  
  # Fill in normal using central differences at r=2:N-1
  nr(n) = n_r(n,@CentralDifference);
  nz(n) = n_z(n,@CentralDifference);
  
  # Derivatives of the normal
  d_nr_dr(n) = d_n_r_dr(n,@CentralDifference);
  d_nz_dr(n) = d_n_z_dr(n,@CentralDifference);
  d2_nr_dr2(n) = d2_n_r_dr2(n,@CentralDifference);
  d2_nz_dr2(n) = d2_n_z_dr2(n,@CentralDifference);
  
  # Fill in bending moments using central differences at r=2:N-1
  M_r{i}(n) = M_r_fct{i}(n,@CentralDifference);
  M_z{i}(n) = M_z_fct{i}(n,@CentralDifference);

  # Derivatives of bending moment tensors
  d_M_r_dr{i}(n) = d_M_r_fct_dr{i}(n,@CentralDifference);
  d_M_z_dr{i}(n) = d_M_z_fct_dr{i}(n,@CentralDifference);
  d2_M_r_dr2{i}(n) = d2_M_r_fct_dr2{i}(n,@CentralDifference);
  d2_M_z_dr2{i}(n) = d2_M_z_fct_dr2{i}(n,@CentralDifference);
  
end

# A shorthand wrapper around CentralDifference
cd = @(order)CentralDifference(order);

# Deflection equation
n = 3:N-2;
R1 = zeros(1,N);

# d_M_i_dr{1}(n)
# (d_ni_dr(n).*M_i{1}(n) + ni(n).*d_M_i_dr{1}(n))
# d2_M_i_dr2{1}(n)
# (d2_ni_dr2(n).*M_i{1}(n) + 2.*d_ni_dr(n).*d_M_i_dr{1}(n) + ni(n).*d2_M_i_dr2{1}(n))

# Missing -1/12.*d2Mzrr/dr2
R1(n) = cd(1)(W,n,h).*d_sigma_dr{1}(n) + sigma{1}(n).*cd(2)(W,n,h) - P - (1./12).*(d2_nz_dr2(n).*M_z{1}(n) + 2.*d_nz_dr(n).*d_M_z_dr{1}(n) + nz(n).*d2_M_z_dr2{1}(n)) ...
       - xi2*(1/12).*((d_nz_dr(n).*M_z{1}(n) + nz(n).*d_M_z_dr{1}(n)).*G_rrr(n) + nz(n).*M_z{1}(n).*d_G_rrr_dr(n) + (d_nz_dr(n).*M_z{3}(n) + nz(n).*d_M_z_dr{3}(n)).*G_rtt(n) + nz(n).*M_z{3}(n).*d_G_rtt_dr(n)) - xi2.*P.*cd(1)(U,n,h);

# In-plane Equation
n = 3:N-2;
R2 = zeros(1,N);

# Missing -1/12.*d2Mrrr/dr2
R2(n) = d_sigma_dr{1}(n)...
     + xi2.*(d_sigma_dr{1}(n).*cd(1)(U,n,h) + sigma{1}(n).*cd(2)(U,n,h) + P.*cd(1)(W,n,h) - (1./12).*(d2_nr_dr2(n).*M_r{1}(n) + 2.*d_nr_dr(n).*d_M_r_dr{1}(n) + nr(n).*d2_M_r_dr2{1}(n)))...
     - (xi4/12).*((d_nr_dr(n).*M_r{1}(n) + nr(n).*d_M_r_dr{1}(n)).*G_rrr(n) + nr(n).*d_G_rrr_dr (n).*M_r{1}(n) + (d_nr_dr(n).*M_r{3}(n) + nr(n).*d_M_r_dr{3}(n)).*G_rtt(n) + nr(n).*M_r{3}(n).*d_G_rtt_dr(n));

# Boundary Conditions
R1(1) = ForwardDifference(1)(W,1,h);
R1(2) = ForwardDifference(3)(W,1,h);
R1(N-1) = BackwardDifference(1)(W,N,h);
R1(N) = W(N);
R2(1) = U(1);
R2(2) = ForwardDifference(2)(U,1,h);
R2(N-1) = BackwardDifference(1)(U,N,h);
R2(N) = U(N);

R = [R1,R2];

endfunction
