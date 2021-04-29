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
 % Find the Error as a function of xi_c for a range of frequencies. 
 %
 %=========================================================================%
 
 % Add noise of the prescribed noiselevel
 NoiseLevel=0.01;
 We=W+randn(size(W))*NoiseLevel;
 Ge=G+randn(size(G))*NoiseLevel;

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
 %print -depsc F2-Error-depends-on-xi_c.eps
 
 %=========================================================================%
 %
 % Find the optimum xi_c as a function of the noise level.
 %
 %=========================================================================%

 
  NoiseLevel=10.^-(1:0.2:5);
  OptimalXi = zeros(size(NoiseLevel));
  OptimalError=zeros(size(NoiseLevel));
   
   for i = 1:length(NoiseLevel)
       
       
       % Add the current level of noise to the data
         We=W+randn(size(W))*NoiseLevel(i);
         Ge=G+randn(size(G))*NoiseLevel(i);
  
        % Set the frequency range to test for 
        xi_c= (0.2:0.2:20);
        Error = zeros(size(xi_c));
        
        % Compute error as a function of xi_c
        for j = 1:length(xi_c)
     
          [T,Tx,x1]=SHESolver( x , t , Ge , We ,  kappa,1,'FFT' ,xi_c(j));
          Fode=T(end,:)';
          Error(j) = norm(Fode -F)/sqrt(N);
        end
        % Find the optimal xi_c for this noiselevel
        [m,k]=min(Error);
        OptimalXi(i)=xi_c(k(1));
        OptimalError(i)=m;
        fprintf('Testing noise level %e. Optimal xi_c=%f\n',NoiseLevel(i),OptimalXi(i))
   end
loglog(NoiseLevel,OptimalError,'LineWidth',1.4);,
xlabel('Noiselevel: \epsilon','FontSize',14);
ylabel('Optimal error ||f(t)-f_{m,\xi_c}(t)||_2','FontSize',14);
%print -depsc F2-Optimal-Error.eps

loglog(NoiseLevel,OptimalXi,'LineWidth',1.4);,
xlabel('Noiselevel: \epsilon','FontSize',14);
ylabel('Optimal \xi_c','FontSize',14);
%print -depsc F2-Optimal-Xi_c.eps

%=========================================================================%
%
% Finally plot two solutions for noise level 10^-2 and for different \xi_c
%
%=========================================================================%

 % Add noise of the prescribed noiselevel
 NoiseLevel=0.01;
 We=W+randn(size(W))*NoiseLevel;
 Ge=G+randn(size(G))*NoiseLevel;
 
 % Set xi_c to the flat part of the error curve where the error is low.
 
 xi_c=6;
 
 % Solve the problem and plot the solution. Also compute the error
 
 [T,Tx,x1]=SHESolver( x , t , Ge , We ,  kappa,1,'FFT' ,xi_c);
 Fode=T(end,:)';
 Error = norm(Fode -F)/sqrt(N);
 fprintf(1,'The error for xi_c=%f is %e\n',xi_c,Error)
 
 plot(t,F,'b--',t,Fode,'k-','LineWidth',1.4);
 xlabel('Time: t','FontSize',14);
 ylabel('Solutions: f(t) and f_{m,\xi_c}(t)','FontSize',14)
 print -depsc F2-Solutions-Good-Xi-6.eps
 
  
 % Set xi_c to the flat part of the error curve where the error is low.
 
 xi_c=18;
 
 % Solve the problem and plot the solution. Also compute the error
 
 [T,Tx,x1]=SHESolver( x , t , Ge , We ,  kappa,1,'FFT' ,xi_c);
 Fode=T(end,:)';
 Error = norm(Fode -F)/sqrt(N);
 fprintf(1,'The error for xi_c=%f is %e\n',xi_c,Error)
 
 plot(t,F,'b--',t,Fode,'k-','LineWidth',1.4);
 xlabel('Time: t','FontSize',14);
 ylabel('Solutions: f(t) and f_{m,\xi_c}(t)','FontSize',14)
 print -depsc F2-Solutions-Bad-Xi-18.eps

 
