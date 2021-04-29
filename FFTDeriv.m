%FFTDeriv: Compute the derivative of a function using the Fourier 
% transform. The number of frequency components used is decided 
% by the parameter xi_c.
%
% Usage:
%   >> dF = FFTDeriv(F , xi_c);
%
% where F is a vector the values of a function for 0<=x<1. Gibbs
% phenomena are avoided by first extending the vector to double
% length using a cubic spline.
%  
function [dF,F]=FFTDeriv(F,xi_c);

 [n1,n2]=size(F); if n2>n1,F=F';,end; % Need column vectors
 
 % Extend the vector to double length first
 F=ExtendVector(F);
 
 N=length(F); ii=sqrt(-1);
 
 % Create the symbol for the derivative with high frequencies removed.
 
 symbol=(0:N-1)'; identity=ones(N,1);
 symbol(ceil(N/2):N)=-(N-ceil(N/2)+1:-1:1);
 ind=find(abs(symbol)>xi_c);symbol(ind)=0;
 identity(ind)=0;
 
 % Apply symbol 
 
 dF=pi*real(ifft(ii*symbol.*fft(F)));F=real(ifft(identity.*fft(F)));
 
 % Reduce the size of F again
 
 dF=dF(1:N/2);F=F(1:N/2);
 
 if n2>n1,dF=dF';F=F';,end; % Original F was row vector
 
%ExtendVector: Extend the function to the interval [0,2] smoothly
% using a spline function. 
function [F]=ExtendVector(F);

 N=length(F);x=(0:N-1)'/N;
 width=min(length(F)/2,6);
%
% Find a least squares fit y=c1+c2*x to F(1:width) and a fit
% y=c1+c2*x to F(N-width+1:N).
%
ca=[x(1:width).^0,x(1:width).^1]\F(1:width);
cb=[x(N-width+1:N).^0,x(N-width+1:N).^1]\F(N-width+1:N);

%
% Now find a cubic spline that match the linear fits.
% 
 xval=[x(1:width);x(N-width+1:N)];
 fval=[cb(1)+cb(2)*(1+x(1:width));ca(1)+ca(2)*(-x(width+1:-1:2))];
 F2=spline(xval,fval,x);
 F=[F;F2];
