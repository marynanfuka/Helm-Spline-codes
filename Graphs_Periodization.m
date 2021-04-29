
 Tend=1;N=500;t=Tend*(0:N-1)'/(N-1);

%
% Create the original non-periodic function.
%
t1=[0 0.15 0.3 0.45 0.48 0.52 0.6 0.7  0.74 0.78  0.8 0.85 0.95 1.0];
y1=[0 0    0.3 1.3  1.4  1.25 1.1 2.73 2.86 2.80  2.3 1.3  1.1 1.8];
pp=spline(t1,y1);
F=ppval(pp,t);

%
% Add noise.
%
 epsilon=0.05;
 F=F+epsilon*randn(size(F))

%
% Extend to double lengh 
%
 N=length(F);
 x=(0:N-1)'/N;
 width=min(length(F)/2,6);
 ca=[x(1:width).^0,x(1:width).^1]\F(1:width);
 cb=[x(N-width+1:N).^0,x(N-width+1:N).^1]\F(N-width+1:N);
 xval=[x(1:width);x(N-width+1:N)];
 fval=[cb(1)+cb(2)*(1+x(1:width));ca(1)+ca(2)*(-x(width+1:-1:2))];
 F2=spline(xval,fval,x);
 F2=[F;F2];
 t2=[t;1+t(2)-t(1)+t];
 
% Compute "exact derivative of the spline.
 dF=ppval(fnder(pp),x); 
 
% Use FFTDeriv on the non-perodic function
xi_c=35;
dF2=FFTDeriv(F,xi_c);


% Make plots

plot(t2,F2,'b--',t,F,'k','LineWidth',1.3)
axis([0 2 -0.5 3.5])
xlabel('Time: t','FontSize',14);
ylabel('Function: f(t)','FontSize',14);
print -depsc F1-Periodization-1.eps

plot(t,dF,'k',t,dF2,'--','LineWidth',1.3)
xlabel('Time: t','FontSize',14);
ylabel('Derivatives: f^\prime(t) and (D_{\xi_c}f)(t)','FontSize',14);
print -depsc F1-Periodization-2.eps

 
 