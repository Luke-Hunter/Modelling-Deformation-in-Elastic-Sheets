%% Central Finite Difference Functions for equally spaced 1D finite differences
% Central Differences
function [cd]=CentralDifference(Order)
switch Order
 case 0
   cd=@(f,n,h)(f(n));
 case 1
   cd=@(f,n,h)(-f(-1 + n) + f(1 + n))/(2*h);
 case 2
   cd=@(f,n,h)(f(-1 + n) - 2*f(n) + f(1 + n))/h^2;
 case 3
  cd=@(f,n,h)(-f(-2 + n) + 2*f(-1 + n) - 2*f(1 + n) + f(2 + n))/(2*h^3);
 case 4
  cd=@(f,n,h)(f(-2 + n) - 4*f(-1 + n) + 6*f(n) - 4*f(1 + n) + f(2 + n))/h^4;
 otherwise
  warning('Unexpected derivative order in CentralDifference function.')
end
