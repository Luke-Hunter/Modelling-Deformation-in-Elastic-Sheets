%% Central Finite Difference Functions for equally spaced 1D finite differences
% Central Differences
function [dcd_df]=DCentralDifferenceDf(Order)
% Kronecker delta
delta=@(m,n) double(m==n);
% Now for the switch
switch Order
 case 0
   dcd_df=@(n,m,h) delta(n,m,h);
 case 1
   dcd_df=@(n,m,h)(-delta(-1 + n,m) + delta(1 + n,m))/(2*h);
 case 2
   dcd_df=@(n,m,h)( delta(-1 + n, m) - 2*delta(n,m) + delta(1 + n,m))/h^2;
 case 3
   dcd_df=@(n,m,h)(-delta(-2 + n, m) + 2*delta(-1 + n,m) - 2*delta(1 + n,m) ...
                +delta(2 + n,m))/(2*h^3);
 case 4
   dcd_df=@(n,m,h)(delta(-2 + n, m) - 4*delta(-1 + n,m) + 6*delta(n,m)...
                -4*delta(1 + n,m) + delta(2 +n,m))/h^4;
 otherwise
  warning('Unexpected derivative order in CentralDifference function.')
end
