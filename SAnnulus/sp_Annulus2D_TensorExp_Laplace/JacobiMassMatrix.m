%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


 function mass = JacobiMassMatrix(degree,alpha, beta)
 
 i = [0:degree];
 
 mass = (2^(alpha+beta+1)./(2*i+alpha+beta+1)).* ...
     gamma(i+alpha+1).*gamma(i+beta+1)./(gamma(i+1).*gamma(i+alpha+beta+1));
 
 mass = mass';