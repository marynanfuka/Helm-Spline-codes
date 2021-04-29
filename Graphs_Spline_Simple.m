%
% 
% First Create the test problem.
%

 % Set parameter values
 a=1;b=5;               % The domain is 0<x<a and 0<t<b.
 M=200;x=a*(0:M-1)'/(M-1);
 N=500; t=b*(0:N-1)'/(N-1);
 kappa=1.2;              % Thermal diffusivity.

 F=[0 0 0 0.1 0.73 0.98 0.83 0.22 0.37 0.77 0.55 0.67 0.40 0.35 0.21 0.23 0.29 0.37 0.45];
 F=spline(linspace(0,b,length(F)),F,t);
 G=zeros(size(F));
 
 % Solve the problem and create the Cauchy data at x=a.
 [H,W]=SHEForwardSolv(x,t,F,G,kappa,1);

 % Plot the exact solution. That is both the temperature at x=0 and the
 % heat-flux at x=a.
 
 plot(t,F,'b--',t,W,'k','LineWidth',1.4);
 xlabel('Time: t','FontSize',14);
 ylabel('Data: u(0,t) and u_x(a,t)','FontSize',14)
 print -depsc F2-Exact-Solution.eps

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
  lambda = 10.^-(5:0.1:11);
  Error = zeros(size(lambda));
  
  for i = 1:length(lambda)
      
  [T,Tx,x2]=SHESolver( x , t , Ge , We ,  kappa ,1,'SSP', lambda(i) );
  Fode=T(end,:)';
  Error(i) = norm(Fode -F)/sqrt(N);
 end
 
 semilogx(lambda,Error,'LineWidth',1.4);
 xlabel('Lambda: \lambda','FontSize',14)
 ylabel('Error: || f(t)-f_{m,\lambda}||_2','FontSize',14); 
% print -depsc F3-Error-depends-on-lambda.eps

 
 %=========================================================================%
 %
 % Find the optimum lambda as a function of the noise level.
 %
 %=========================================================================%

 
  NoiseLevel=10.^-(1:0.2:5);
  Optimallambda = zeros(size(NoiseLevel));
  OptimalError=zeros(size(NoiseLevel));
   
   for i = 1:length(NoiseLevel)
       
       
       % Add the current level of noise to the data
         We=W+randn(size(W))*NoiseLevel(i);
         Ge=G+randn(size(G))*NoiseLevel(i);
  
        % Set the range of values of lambda to test for 
        %lambda= (0.2:0.2:20);
        lambda = 10.^-(4:0.2:11);
        Error = zeros(size(lambda));
        
        % Compute error as a function of lambda
        for j = 1:length(lambda)
     
          [T,Tx,x1]=SHESolver( x , t , Ge , We ,  kappa,1,'SSP' ,lambda(j));
          Fode=T(end,:)';
          Error(j) = norm(Fode -F)/sqrt(N);
        end
        % Find the optimal lambda for this noiselevel
        [m,k]=min(Error);
        Optimallambda(i)=lambda(k(1));
        OptimalError(i)=m;
        fprintf('Testing noise level %e. Optimal lambda=%f\n',NoiseLevel(i),Optimallambda(i))
   end
loglog(NoiseLevel,OptimalError,'LineWidth',1.4);
xlabel('Noiselevel: \epsilon','FontSize',14);
ylabel('Optimal error ||f(t)-f_{m,\lambda}(t)||_2','FontSize',14);
print -depsc F3-Optimal-Lambda-Error.eps

loglog(NoiseLevel,Optimallambda,'LineWidth',1.4);
xlabel('Noiselevel: \epsilon','FontSize',14);
ylabel('Optimal \lambda','FontSize',14);
print -depsc F3-Optimal-lambda.eps

%=========================================================================%
%
% Finally plot one solution for noise level 10^-2 and for a good lambda
% 
%
%=========================================================================%

 % Add noise of the prescribed noiselevel
 NoiseLevel=0.01;
 We=W+randn(size(W))*NoiseLevel;
 Ge=G+randn(size(G))*NoiseLevel;
 
 % Set lambda to the flat part of the error curve where the error is low.
 
 lambda=2e-8;
 
 % Solve the problem and plot the solution. Also compute the error
 
 [T,Tx,x1]=SHESolver( x , t , Ge , We ,  kappa,1,'SSP' ,lambda);
 Fode=T(end,:)';
 Error = norm(Fode -F)/sqrt(N);
 fprintf(1,'The error for lambda=%f is %e\n',lambda,Error)
 
 plot(t,F,'b--',t,Fode,'k-','LineWidth',1.4);
 xlabel('Time: t','FontSize',14);
 ylabel('Solutions: f(t) and f_{m,\lambda}(t)','FontSize',14)
 print -depsc F3-Solutions-Good-lambda-2e-8.eps
 
  