% SHESystem: Implements the derivative of the vector V= [ T , kappa*Tx ] 
% where U is the solution to the sideways heat equation. This function is 
% supposed to be used with the MATLAB ODE solvers, e.g. ode45.
%
function [ Vx ] = SHESystem( x , V , events , a_fun , b_fun , regu_method,regu_param);

%
% Extract U(x,:) and Ux(x,:).
%
 n=length(V)/2;
 T=V(1:n);Tx=V(n+1:2*n)./a_fun(T);
 
 
%
% Compute the vector Vx = [ Ux ; D*U ] using FFTDeriv to implement the
% product D*U.
%
 if strcmp(regu_method,'FFT')
  Tt=FFTDeriv(T,regu_param);
 elseif strcmp(regu_method,'SSP')
  Tt=SSPDeriv(T,regu_param);
 else
   error('SHESystem: Unknown regularization method');
 end
 
 Vx = [ Tx ; b_fun(T).*Tt ];
 

 
 
