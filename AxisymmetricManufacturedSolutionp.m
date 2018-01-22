%% Script to use a Manufactured solution to test the perturbed Foeppl von Karman
%  Equations. This uses a homegenous axisymmetic solution so it does not test
%  the axisymmetric solution. 
%  Author: David Robinson
%  Date 9/10/2017
% Test parameters
# edited to account for extra terms from modified strain tensor
nlc = 1.0;
N=800;
xi = 0.00001; # xi = h^2/R^2

% Manufactured parameters
nu = 1/3.;
mu = -9 / 4.;

# terms including xi are added by accounting for the modified strain
p = @(r) 64 - nlc * 6*(-2+r.^2.*(8 - 9*r.^2 + 4*(-1+r.^2).^2.*(-2+5*r.^2)*nlc))...
          - xi * nlc.^2 *3*r.*(-1 + r.^2).*(-3 + 20*r.^2 - 51*r.^4 + 46*r.^6);
q = @(r) 3*r.*(-3+4*nlc*(2-7*r.^2+5*r.^4))...
          + (3/8)*nlc*xi*(-3 +71*r.^2 -249*r.^4 +205*r.^6);
w = @(r) (1-r.^2).^2;
u = @(r) r.*(1-r.^2);

% Now Loop
[wapx,~,info,infoarray,~]=fsolve(@(x)residual_FvKd_perturbed_varP(x,p,q,nlc,mu,nu,xi), zeros(1,2*N),optimset('Jacobian','on'));
info
GetL2Norm(wapx(1:N),w)
GetL2Norm(wapx(N+1:2*N),u)
r=linspace(0,1,N);
%% Test the manufactured solution
R=residual_FvKd_perturbed_varP([w(r),u(r)],p,q,nlc,mu,nu,xi);
newwapx = wapx
% Plot the two solutions
plot(r,wapx(1:N))
hold all;
plot(r,w(r),'.')
xlabel 'r'
ylabel 'w(r)'

figure
plot(r,wapx(N+1:2*N))
hold all;
plot(r,u(r),'.')
xlabel 'r'
ylabel 'u(r)'

figure
plot(R)
xlabel 'n'
ylabel 'R_n'
hold all;

%% Second part: now loop to find the scaling of the norm with h
narray=50:50:500;
for(i=1:length(narray))
 % Get N
 N=narray(i);
 % Get the solution
[wapx,~,info,infoarray,~]=fsolve(@(x)residual_FvKd_perturbed_varP(x,p,q,nlc,mu,nu,xi),...
  zeros(1,2*N),optimset('Jacobian','on'));
 info

 % Get the norms 
 [~,normW2(i)]=GetL2Norm(wapx(1:N),w);
 [~,normU2(i)]=GetL2Norm(wapx(N+1:2*N),u);

end

% Now plot the resulting error graph
figure
param=polyfit(log10(1./(narray-1)),0.5*log10(normW2+normU2),1);
fit=@(x)param(1)*x+param(2);

plot(log10(1./(narray-1)),0.5*log10(normW2+normU2),'.b')
hold all
plot(log10(1./(narray-1)),fit(log10(1./(narray-1))),'-b')
legend_entry = sprintf('y = %f x + %f',param(1),param(2));
legend('data',legend_entry)
xlabel 'log10(h)'
ylabel '||e||'
