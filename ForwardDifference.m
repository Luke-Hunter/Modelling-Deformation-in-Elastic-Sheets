%% Forward Finite Difference Functions for equally spaced 1D finite differences
% Forward Differences
function [fd] = ForwardDifference(Order)
switch Order
 case 0
   fd=@(f,n,h)(f(n));
 case 1
   fd=@(f,n,h)-(3*f(n) - 4*f(1 + n) + f(2 + n))/(2*h);
 case 2
   fd=@(f,n,h)(2*f(n) - 5*f(1 + n) + 4*f(2 + n) - f(3 + n))/h^2;
 case 3
   fd=@(f,n,h)-(5*f(n) - 18*f(1 + n) + 24*f(2 + n) - 14*f(3 + n) +...
                 3*f(4 + n))/(2*h^3);
 case 4
   fd=@(f,n,h)(3*f(n) - 14*f(1 + n) + 26*f(2 + n) - 24*f(3 + n) +...
             11*f(4 + n) - 2*f(5 + n))/h^4;
 otherwise
  warning('Unexpected derivative Order in ForwardDifference function.')
end
