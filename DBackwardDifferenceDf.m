%% Backward Finite Difference Functions for equally spaced 1D finite differences
% Backward Differences
function [dbd_df]=DBackwardDifferenceDf(Order)
% Kronecker delta
delta=@(m,n) double(m==n);
% Now for the switch
switch Order
 case 0
   dbd_df=@(n,m,h) delta(n,m,h);
 case 1
   dbd_df=@(n,m,h) (delta(-2 + n,m) - 4*delta(-1 + n,m) + 3*delta(n,m))/(2*h);
 case 2
   dbd_df=@(n,m,h)-((delta(-3 + n,m) - 4*delta(-2 + n,m) + 5*delta(-1 + n,m)...
         - 2*delta(n,m))/h^2);
 case 3
   dbd_df=@(n,m,h)(-delta(-2 + n, m) + 2*delta(-1 + n,m) - 2*delta(1 + n,m) ...
                +delta(2 + n,m))/(2*h^3);

 case 4
   dbd_df=@(n,m,h)(3*delta(-4 + n,m) - 14*delta(-3 + n,m) + 24*delta(-2 + n,m)...
         - 18*delta(-1 + n,m) + 5*delta(n,m))/(2*h^3);
 otherwise
  warning('Unexpected derivative order in BackwardDifference function.')
end
