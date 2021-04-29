
%
% 
% First Create the test problem.
%

 % Set parameter values
 a=1;b=5;               % The domain is 0<x<a and 0<t<b.
 M=200;x=a*(0:M-1)'/(M-1);
 N=500; t=b*(0:N-1)'/(N-1);
 
 % Set thermal diffusivity
% kappa= @(x)(1+ sin(pi*(a-x)/a).^0.8 + 4.2*exp(-27.3*(a/3-x).^2) );
 kappa= @(u)(1+ sin(pi*u).^0.8 + 4.2*exp(-27.3*(a/3-u).^2) );
 % Step function as exact solution
 
 F=zeros(size(t))+(t>b/2).*ones(size(t));
 ind=find( 8*b/17 < t & t < 9*b/17 );
 F(ind)=sin( pi*(t(ind)-t(ind(1)))/(t(ind(end))-t(ind(1)))/2 ).^2; 
 G=zeros(size(F));
 
 % Solve the problem and create the Cauchy data at x=a.
 [H,W]=SHEForwardSolv(x,t,F,G,kappa,1);

 % Plot the exact solution. That is both the temperature at x=0 and the
 % heat-flux at x=a.
 
 plot(t,F,'b-',t,W,'k-','LineWidth',1.4);
 xlabel('Time: t','FontSize',14);
 ylabel('Data: u(0,t)','FontSize',14)

 %=========================================================================%
 %
 % Find the Error as a function of lambda for a range of values of lambda. 
 %
 %=========================================================================%
 
 % Add noise of the prescribed noiselevel
 NoiseLevel=0.01;
 We=W+randn(size(W))*NoiseLevel;
 Ge=G+randn(size(G))*NoiseLevel;
% Decide the range of lambda to use. Solve the problem for
 % each value of lambda and record the error.

  lambda = 10.^-(6:0.1:11);
  Error = zeros(size(lambda));
  
  for i = 1:length(lambda)
      
  [T,Tx,x2]=SHESolver( x , t , Ge , We ,  kappa ,1,'SSP', lambda(i) );
  Fode=T(end,:)';
  Error(i) = norm(Fode -F)/sqrt(N);
 end
 
 semilogx(lambda,Error,'LineWidth',1.4);
 xlabel('Lambda: \lambda','FontSize',14)
 ylabel('Error: || f(t)-f_{m,\lambda}||_2','FontSize',14); 
 %print -depsc F4-Error-depends-on-lambda.eps
 
 % Decide the range if cut off frequencies to use. Solve the problem for
 % each frequency and record the error.
 xi_c= (0.1:0.2:20);
 Error = zeros(size(xi_c));
 for i = 1:length(xi_c)
      
    [T,Tx,x1]=SHESolver( x , t , Ge , We ,  kappa,1,'FFT' ,xi_c(i));
    Fode=T(end,:)';
    Error(i) = norm(Fode -F)/sqrt(N);
 end
 
 plot(xi_c,Error,'LineWidth',1.4);
 xlabel('Frequency: \xi_c','FontSize',14)
 ylabel('Error: || f(t)-f_{m,\xi_c}||_2','FontSize',14); 
 %print -depsc F4-Error-depends-on-Xic.eps
 


%=========================================================================%
%
% Finally plot two solutions for noise level 10^-2 and for a good
% lambda and \xi_c
%
%=========================================================================%

 % Add noise of the prescribed noiselevel
 NoiseLevel=0.01;
 We=W+randn(size(W))*NoiseLevel;
 Ge=G+randn(size(G))*NoiseLevel;
 
 % Set lambda to the flat part after the worst instabilities are removed
 
 lambda=1e-8
 
 % Solve the problem and plot the solution. Also compute the error
 
 [T,Tx,x1]=SHESolver( x , t , Ge , We ,  kappa,1,'SSP' ,lambda);
 Fode=T(end,:)';
 Error = norm(Fode -F)/sqrt(N);
 fprintf(1,'The error for lambda=%e is %e\n',lambda,Error)
 
 plot(t,F,'b--',t,Fode,'k-','LineWidth',1.4);
 xlabel('Time: t','FontSize',14);
 ylabel('Solutions: f(t) and f_{m,\lambda}(t)','FontSize',14)
 %print -depsc F4-Step-solution-SSP-lambda-1e-8.eps

 
 % Set xi_c to the flat part of the error curve where the error is low.
 
 xi_c=7;
 
 % Solve the problem and plot the solution. Also compute the error
 
 [T,Tx,x1]=SHESolver( x , t , Ge , We ,  kappa,1,'FFT' ,xi_c);
 Fode=T(end,:)';
 Error = norm(Fode -F)/sqrt(N);
 fprintf(1,'The error for xi_c=%f is %e\n',xi_c,Error)
 
 plot(t,F,'b--',t,Fode,'k-','LineWidth',1.4);
 xlabel('Time: t','FontSize',14);
 ylabel('Solutions: f(t) and f_{m,\xi_c}(t)','FontSize',14)
 %print -depsc F4-Step-solution-FFT-xic7.eps
 
 
%=========================================================================%
%
% Use noise 10^-2 and run experiments and measure the error for 0<t<1
% lambda and \xi_c
%
%=========================================================================%

 m=100;
 Error=zeros(m,2);
 lambda=1e-8;xi_c=7;
 ind=find(t<1);
 
 for i=1:m
   i
   % Add noise of the prescribed noiselevel
   NoiseLevel=0.01;
   We=W+randn(size(W))*NoiseLevel;
   Ge=G+randn(size(G))*NoiseLevel;
 

  % Solve the problem and compute the errors
 
  [T,Tx,x1]=SHESolver( x , t , Ge , We ,  kappa,1,'SSP' ,lambda);
  Fode=T(end,:)';
  Error(i,1)=norm(Fode(ind) -F(ind))/sqrt(N);

 
  [T,Tx,x1]=SHESolver( x , t , Ge , We ,  kappa,1,'FFT' ,xi_c);
  Fode=T(end,:)';
  Error(i,2) = norm(Fode(ind) -F(ind))/sqrt(N);

 end
 
 plot(1:m,Error(:,1),'kx',1:m,Error(:,2),'bo','MarkerSize',8,'LineWidth',1.4)
 xlabel('Experiment number: k','FontSize',14)
 ylabel('Errors ||f-f_{m_k,\lambda}||_{2,I} and ||f-f_{m_k,\xi_c}||_{2,I}','FontSize',14)
 %print -depsc F4-Step-Errors-same-noise-1e-2.eps
 
 
 % Redo but with a smaller noiselevel
 
 m=100;
 Error=zeros(m,2);
 lambda=2e-9;xi_c=8.3;
 ind=find(t<1);
 
 for i=1:m
   i
   % Add noise of the prescribed noiselevel
   NoiseLevel=0.001;
   We=W+randn(size(W))*NoiseLevel;
   Ge=G+randn(size(G))*NoiseLevel;
 

  % Solve the problem and compute the errors
 
  [T,Tx,x1]=SHESolver( x , t , Ge , We ,  kappa,1,'SSP' ,lambda);
  Fode=T(end,:)';
  Error(i,1)=norm(Fode(ind) -F(ind))/sqrt(N);

 
  [T,Tx,x1]=SHESolver( x , t , Ge , We ,  kappa,1,'FFT' ,xi_c);
  Fode=T(end,:)';
  Error(i,2) = norm(Fode(ind) -F(ind))/sqrt(N);

 end
 
 plot(1:m,Error(:,1),'kx',1:m,Error(:,2),'bo','MarkerSize',8,'LineWidth',1.4)
 xlabel('Experiment number: k','FontSize',14)
 ylabel('Errors ||f-f_{m_k,\lambda}||_{2,I} and ||f-f_{m_k,\xi_c}||_{2,I}','FontSize',14)
 %print -depsc F4-Step-Errors-same-noise-1e-3.eps
 