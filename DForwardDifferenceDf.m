%% Forward Finite Difference Functions for equally spaced 1D finite differences
% Forward Differences
function [dfd_df]=DForwardDifferenceDf(Order)
% Kronecker delta
delta=@(m,n) double(m==n);
% Now for the switch
switch Order
 case 0
   dfd_df=@(n,m,h) delta(n,m,h);
 case 1
   dfd_df=@(n,m,h)-(3*delta(n,m) - 4*delta(1 + n,m) + delta(2 + n,m))/(2*h);
 case 2
   dfd_df=@(n,m,h)(2*delta(n,m) - 5*f(1 + n,m) + 4*delta(2 + n,m)...
       - delta(3 + n,m))/h^2;
 case 3
   dfd_df=@(n,m,h)-(5*delta(n,m) - 18*delta(1 + n,m) + 24*delta(2 + n,m)...
     - 14*delta(3 + n,m) + 3*delta(4 + n,m))/(2*h^3);
 case 4
   dfd_df=@(n,m,h)(3*delta(n,m) - 14*delta(1 + n,m) + 26*delta(2 + n,m)...
      - 24*delta(3 + n,m) + 11*delta(4 + n,m) - 2*delta(5 + n,m))/h^4;
 otherwise
  warning('Unexpected derivative order in ForwardDifference function.')
end
