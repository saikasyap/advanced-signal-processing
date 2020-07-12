function cossq=gencos(a,b,C)
%% GENCOS:   Fxn to compute generalized cosine between two vectors
%   cossq=gencos(a,b,C)
%
%   Variables:  a = vector 1
%               b = vector 2
%               C = Scaling matrix [OPTIONAL]
%

%--------------------------------------------------------------------------
if nargin<3
  cossq=abs(a'*b).^2/((a'*a)*(b'*b));
else  
  cossq=abs(a'*C*b).^2/((a'*C*a)*(b'*C*b));
end
%--------------------------------------------------------------------------
return
