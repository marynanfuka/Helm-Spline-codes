%SSPDeriv: Compute the derivative of a function using the Cubic Smoothing spline.  
% The regularization parameter is lambda.
% Usage:
%   >> dF = SSPDeriv(F , lamnbda );
%
% where F is a vector the values of a function for 0<=x<1. 
%  
function [dF,F]=SSPDeriv(F,lambda);

 [n1,n2]=size(F); if n2>n1,F=F';,end; % Need column vectors
 
 n=length(F);t=(0:n-1)'/(n-1);h=t(2)-t(1);
 
 % Compute smoothing spline
  
 pp=csaps(t,F,h/(lambda+h));
 dF=ppval( fnder(pp) , t ); 
 F=ppval(pp,t);
 
 if n2>n1,dF=dF';,end; % Original F was row vector
