% Function that returns the stress from the strain using a linear constitutive 
% relation. Using cell array structure with function handles.
% This works as a generic interface for handles of up to three variables.

function A = kirchhoff_st_venant_constitutive_cell(E,nu,B)

% A_{jk} = E ( nu B_{ii}\delta_{jk}  + (1-\nu)B_{jk})
% Or :
% A_{ij} = \Gamma_{ijkl} B_{kl}
% with \Gamma_{ijkl} = E/(1-\nu^2) (\nu \delta_{ij}\delta_{jk} \
%      + (1-\nu) \delta_{ik} \delta_{jl})

% Test that the input arguments are in the correct format
if(nargin(B{1}) == nargin(B{2}) && nargin(B{1})== nargin(B{3}))

 if(nargin(B{1}) == 1) 
  % Now compute the 'stress-like' matrix from the 'strain-like' matrix
  A1= @(n)(B{1}(n)+nu*B{3}(n)) ;
  A2= @(n)(1-nu)*B{2}(n) ;
  A3= @(n)(B{3}(n)+nu*B{1}(n)) ;
 
 elseif(nargin(B{1}) == 2)
  % Now compute the 'stress-like' matrix from the 'strain-like' matrix (2 args)
  A1= @(n,m)(B{1}(n,m)+nu*B{3}(n,m)) ;
  A2= @(n,m)(1-nu)*B{2}(n,m) ;
  A3= @(n,m)(B{3}(n,m)+nu*B{1}(n,m)) ;
 
 elseif(nargin(B{1}) == 3)
  % Now compute the 'stress-like' matrix from the 'strain-like' matrix (3 args)
  A1= @(n,m,l)(B{1}(n,m,l)+nu*B{3}(n,m,l)) ;
  A2= @(n,m,l)(1-nu)*B{2}(n,m,l) ;
  A3= @(n,m,l)(B{3}(n,m,l)+nu*B{1}(n,m,l)) ;

 elseif(nargin(B{1}) == 4)
  % Now compute the 'stress-like' matrix from the 'strain-like' matrix (3 args)
  A1= @(n,m,l,k)(B{1}(n,m,l,k)+nu*B{3}(n,m,l,k)) ;
  A2= @(n,m,l,k)(1-nu)*B{2}(n,m,l,k) ;
  A3= @(n,m,l,k)(B{3}(n,m,l,k)+nu*B{1}(n,m,l,k)) ;
 
 else 
  error('Not implemented for handles with more than four arguments.');
 end
 
% Otherwise throw an error
else
 error('Function handles in cell array have different numbers of input arguments');
end

A={A1,A2,A3};
end
