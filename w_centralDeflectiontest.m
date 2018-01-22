function [w_cd] = w_centralDeflectiontest
  
  
  
  for (i=1:100)
    
    F = 650*i;
    
    [w,~,info,~,~] = fsolve(@(x)residual_FvKd_perturbed(x,F,1,-2.25,0.33333,0.00000000001),zeros(1,800),optimset('Jacobian','on','TolX',1E-8,'TolFun',1E-8))
    
    ww(i) = w(1)
    P(i) = F
  end
  
  plot(P,ww)
  
endfunction