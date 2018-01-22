function [w_cd] = bendingw_centralDeflectiontest
  
  
  
  for (i=1:10)
    
    F = -6000*i;
    
    [w,~,info,~,~] = fsolve(@(x)residual_bending_axisymmetric(x,F,0.5,0.01),zeros(1,200))
    
    ww(i) = w(1)
    P(i) = -F
  end
  
  plot(P,ww)
  
  w_cd = ww
  
endfunction