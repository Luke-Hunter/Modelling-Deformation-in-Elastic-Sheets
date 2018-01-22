function [norm,norm2] =GetL2Norm(soln,exact)
  % Extract exact
  wt=exact;
  % Extact aproximate
  N=length(soln);
  wa=soln;
  % Construct r and w_inter
  r=linspace(0,1,N);
  w_inter=@(x)interp1(r,wa,x,'*linear');
  % Now integrate the diff
  norm2=quad(@(x)(exact(x)-w_inter(x)).^2,0,1);
  norm=sqrt(norm2);
endfunction
