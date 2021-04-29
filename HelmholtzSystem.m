% SHESystem: Implements the second derivative of the vector V= [ U , Uy ] 
% where U is the solution to the sideways heat equation. This function is 
% supposed to be used with the MATLAB ODE solvers, e.g. ode45.
%
function [Vy ] = HelmholtzSystem(y,V, events ,k2,regu_method,regu_param);

%
% Extract U(:,y) and Uy(:,y). Also use the fact that both U and Uy are
% supposed to be zero at the boundary. 
%
 n=length(V)/2;
 U=V(1:n);
 Uy=V(n+1:n*2);
%
% Compute the vector Vy 
%
 if strcmp(regu_method,'SSP')
  Uxx=SSPDeriv(U,regu_param,2);
 else
  error('HelmholtzSystem: Unknown regularization method');
 end
 
  if isnumeric(k2),
      k2U=k2*U;;
  else,
      dx=1/(n-1);x=dx*(0:n-1)';
      k2U=zeros(size(U));
      for i = 1:n;
        k2U=k2(x(i),y)*U(i);
      end
  end 


    
%
% Compute the vector Vy = [ Uy ; D^2*U ] using SSPDeriv2 to implement the
% product D^2*U.
 Uyy=-k2U-Uxx;
 Uyy(1)=0;Uyy(end)=0;Uy(1)=0;Uy(end)=0;
 Vy = [ Uy ; Uyy ]; 
 
 

 
 
