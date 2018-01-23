%% Backward Finite Difference Functions for equally spaced 1D finite differences
% Backward Differences
function [bd] = BackwardDifference(Order)
switch Order
 case 0
   bd=@(f,n,h)(f(n));
 case 1
   bd=@(f,n,h)(f(-2 + n) - 4*f(-1 + n) + 3*f(n))/(2*h);
 case 2
   bd=@(f,n,h)-((f(-3 + n) - 4*f(-2 + n) + 5*f(-1 + n) - 2*f(n))/h^2);
 case 3
   bd=@(f,n,h)(3*f(-4 + n) - 14*f(-3 + n) + 24*f(-2 + n) - 18*f(-1 + n) +...
                5*f(n))/(2*h^3);
 case 4
   bd=@(f,n,h)(-2*f(-5 + n) + 11*f(-4 + n) - 24*f(-3 + n) + 26*f(-2 + n) -...
                14*f(-1 + n) + 3*f(n))/h^4;
 otherwise
  warning('Unexpected derivative Order in BackwardDifference function.')
end
