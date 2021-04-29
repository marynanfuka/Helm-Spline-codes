%=========================================================================%
% Test: Do a simple numerical test with a known solution and a known coefficient
% a(T), b(T). Note that the randon noiuse is different every time and thus
% the graphs may differ from the ones in the paper.
%
 clear,close all,clf

 % Create the data for the text problem and also visualize the solution.
 
 a=2;b=5;                  % The domain is 0<x<a and 0<t<b.
 M=200;x=a*(0:M-1)'/(M-1);
 N=500; t=b*(0:N-1)'/(N-1);
 [X,T]=meshgrid(x,t);
 
 kappa=@(T)1+sin(T/10)/3;                 % Thermal diffusivity.

 % Create the exact data function
 F=[10 10 10 10 9 6 5 5 7 6 6 9 23 29 29 30 28 27 29 29 30 29 13 11 10 10 10 11 10 9 10 10 10 10 ];
 F=spline(linspace(0,b,length(F)),F,t);
 G=10*ones(size(F));
 
 % Solve and compute a surface representation of the solution
 [~,~,Tsurf]=SHEForwardSolv(x,t,F,G,kappa,1);


 % Now obtain numerical data at x=1 and reset space grid.
 [~,k]=min(abs(x-1));G=Tsurf(:,k);a=1;x=a*(0:M-1)'/(M-1);
  
 % Solve and compute a surface representation of the solution
 [~,H]=SHEForwardSolv(x,t,F,G,kappa,1);
 
 
  % Add noise of the prescribed noiselevel
 NoiseLevel=0.1;
 He=H+randn(size(H))*NoiseLevel;
 Ge=G+randn(size(G))*NoiseLevel;
 
 %===== Figure 1: Display the numerical solution and the data ============= 
 
 % Display a surface of the solution 
 mesh(X,T,Tsurf);
 zlabel('Temperature: T(x,t)','FontSize',14);
 xlabel('Position: x','FontSize',14)
 ylabel('Time: t','FontSize',14)
 %print -depsc F1-Test-Temperature-Surface.eps
 

 % Plot the data T,Tx for x=a.
 
 plot(t,Ge,'b-',t,He,'k-','LineWidth',1.4);
 xlabel('Time: t','FontSize',14);
 ylabel('Cauchy Data: T(a,t) and T_x(a,t)','FontSize',14)
 %print -depsc F1-Test-Cauchy-data.eps
 
 
 
 % Now compute the error as a function of the regularization parameter
 % \lambda and also the residual; and make a plot of the L-curve
 
  
  lambda= 10.^-(6:0.05:11);
  Errors = zeros(size(lambda));
  Residuals = zeros(size(lambda));
  SolNorm = zeros(size(lambda));
    
  for i = 1:length(lambda)      
   [T,Tx,x1]=SHESolver( x , t , Ge , He ,  kappa ,1,'SSP', lambda(i) );
   Fode=T(end,:)';
   Errors(i) = norm(Fode -F)/sqrt(N);
   SolNorm(i)=norm(Fode(3:end)-2*Fode(2:end-1)+ Fode(1:end-2))/sqrt(N);
   [~,Gtmp]=SSPDeriv(Ge,lambda(i));
   [~,Htmp]=SHEForwardSolv(x,t,Fode,Gtmp,kappa,1);
   Residuals(i) = norm(He-Htmp)/sqrt(N);
  end
  
%========= Figure 2: Plot the error as a function of \lambda and also the L-curve =================== 

semilogx(lambda,Errors,'b','LineWidth',1.4);
xlabel('Parameter: \lambda','FontSize',14)
ylabel('Error: ||f-f_\lambda^\delta||_2','FontSize',14); 
hold on,k=61;plot(lambda(k),Errors(k),'ro'),hold off
%print -depsc F2-Error-vs-Lambda.eps


loglog(Residuals, SolNorm,'LineWidth',1.4);axis([0.11 0.4 0.01 1])
xlabel('Residual norm: ||\partial_xv^\delta_\lambda(a,\cdot)-h_\delta||_2','FontSize',14);
ylabel('Solution norm: |f^\delta_\lambda|_2','FontSize',14);
hold on,k=60;plot(Residuals(k),SolNorm(k),'ro'),hold off
%print -depsc F2-L-curve.eps
 
  
 % Solve the problem and plot the solution that corresponds to the minimum
 % error.
 
 [T,Tx,x1]=SHESolver( x , t , Ge , He ,  kappa,1,'SSP' ,lambda(61)); %Minimum error value
 Fode=T(end,:)';
 plot(t,F,'b',t,Fode,'k-','LineWidth',1.4);
 xlabel('Time: t','FontSize',14);
 ylabel('Solutions: f(t) and f _\lambda^\delta(t)','FontSize',14)
 %print -depsc F2-Solution-Min-Error.eps
 
 
 [T,Tx,x1]=SHESolver( x , t , Ge , He ,  kappa,1,'SSP' ,lambda(60));%L-curve value
 Fode=T(end,:)';
 plot(t,F,'b',t,Fode,'k-','LineWidth',1.4);
 xlabel('Time: t','FontSize',14);
 ylabel('Solutions: f(t) and f _\lambda^\epsilon(t)','FontSize',14)
 %print -depsc F2-Solution-L-curve.eps

 
 [T,Tx,x1]=SHESolver( x , t , Ge , He ,  kappa,1,'SSP' ,1e-11);
 Fode=T(end,:)';
 plot(t,F,'b',t,Fode,'k-','LineWidth',1.4);
 xlabel('Time: t','FontSize',14);
 ylabel('Solutions: f(t) and f _\lambda^\delta(t)','FontSize',14)
 %print -depsc F2-Solution-Small-lambda.eps
 
 [T,Tx,x1]=SHESolver( x , t , Ge , He ,  kappa,1,'SSP' ,1e-6);
 Fode=T(end,:)';
 plot(t,F,'b',t,Fode,'k-','LineWidth',1.4);
 xlabel('Time: t','FontSize',14);
 ylabel('Solutions: f(t) and f _\lambda^\epsilon(t)','FontSize',14)
 %print -depsc F2-Solution-Large-lambda.eps
 

 
 %=========================================================================%
 %
 % Compute the optimal regularization parameter for various different noise
 % levels and illustrate the ill-posedness.
 %
 %
 close all,clear,clear global

 % Create the data for the text problem and also visualize the solution.
 
 a=2;b=5;                  % The domain is 0<x<a and 0<t<b.
 M=200;x=a*(0:M-1)'/(M-1);
 N=500; t=b*(0:N-1)'/(N-1);
 [X,T]=meshgrid(x,t);
 
 kappa=@(T)1+ sin(T/10)/3;        % Thermal diffusivity.

 % Create the exact data function
 F=[10 10 10 10 9 6 5 5 7 6 6 9 23 29 29 30 28 27 29 29 30 29 13 11 10 10 10 11 10 9 10 10 10 10 ];
 F=spline(linspace(0,b,length(F)),F,t);
 G=10*ones(size(F));
 
 % Solve and compute a surface representation of the solution
 [~,~,Tsurf]=SHEForwardSolv(x,t,F,G,kappa,1);

 % Now obtain numerical data at x=1 and reset space grid.
 [~,k]=min(abs(x-1));G=Tsurf(:,k);a=1;x=a*(0:M-1)'/(M-1);
 [~,H]=SHEForwardSolv(x,t,F,G,kappa,1);
 
 
% Vary the noise level and find both the optimal lambda as well as the
% resulting errors.

  NoiseLevel=10.^-(1:0.2:5);
  Optimallambda = zeros(size(NoiseLevel));
  OptimalError=zeros(size(NoiseLevel));
   
   for i = 1:length(NoiseLevel)
              
     He=H+randn(size(H))*NoiseLevel(i); % Add noise
     Ge=G+randn(size(G))*NoiseLevel(i);
  
     lambda= 10.^-(6:0.04:15);
     Error = zeros(size(lambda));
         
     % Compute error as a function of lambda. Also find optimal lambda.
     for j = 1:length(lambda) 
        [T,Tx,x1]=SHESolver( x , t , Ge , He ,  kappa,1,'SSP' ,lambda(j));
        Fode=T(end,:)';
        Error(j) = norm(Fode -F)/sqrt(N);
     end
     [m,k]=min(Error);
     Optimallambda(i)=lambda(k(1));
     OptimalError(i)=m;
     fprintf('Testing noise level %e. Optimal lambda=%f\n',NoiseLevel(i),Optimallambda(i))
  end
loglog(NoiseLevel,Optimallambda,'LineWidth',1.4);
xlabel('Noiselevel: \delta','FontSize',14);
ylabel('Optimal parameter: \lambda','FontSize',14);
%print -depsc F3-Optimal-lambda.eps

loglog(NoiseLevel,OptimalError,'LineWidth',1.4);
xlabel('Noise level: \delta','FontSize',14);
ylabel('Optimal error: ||f-f^\delta_\lambda||_2','FontSize',14);
%print -depsc F3-Optimal-Error.eps

 
 
 
%===========================================================================
% Application. A heat treatment experiment.
%

%
% The file contains data from the industrial plant. The data is TC1, TC2,
% and TC3.  They correspond to temperatures 5, 11, and 17 mm below the surface of the steel. 
% 
 clear,close all
 load FurnaceHeating.mat

 % Set parametrars and define constants. To save time we don't use all the
 % data. 
 a=4;
 TC1=TC1(1:a:end);TC2=TC2(1:a:end);TC3=TC3(1:a:end);
 N=length(TC1);t=linspace(0,a*N/10,N); % The sampling rate is 10 Hz. 

 L=[5 11 17]*1e-3;            % Distance between measurement points and the 
                              % surface in mm.
 n=500;x=linspace(0,L(3),n)'; % Create the space discretization. We use L(3)
 ind=zeros(size(L));for i=1:length(L),[~,k]=min(abs(L(i)-x));,ind(i)=k(1);,end   
 
 % Plot the measurements. We are going to use TC1 and TC2 to compute the
 % surface temperature. T
 plot(t,TC1,'b-',t,TC2,'k-','LineWidth',1.4)
 xlabel('Time: t [s]','FontSize',14)
 ylabel('Measurements: TC_1, TC_2 [^oC]','FontSize',14);
 %print -depsc Application-Data-TC1-2.eps
 
 % Display the coefficients as functions of T.
 
 T=0:5:800;
 plot(T,ThermalConductivity(T),'LineWidth',1.2)
 ylabel('Thermal conductivity: k [W/m ^oC]','FontSize',14)
 xlabel('Temperature: T [^oC]','FontSize',14)
 %print -depsc Application-Coef-Kappa.eps
 
 plot(T,Density_SpecificHeat(T),'LineWidth',1.2)
 ylabel('Coefficient: \rho{\cdot}c_p [J/m^3 {}^oC]','FontSize',14)
 xlabel('Temperature: T [^oC]','FontSize',14)
 %print -depsc Application-Coef-Rho-Cp.eps
  
  
 % Solve the direct problem in the interval [L1,L2] and compute the
 % temperature derivative at L1. 
 [H1,~,~]=SHEForwardSolv(x(ind(1):ind(2)),t,TC1,TC2,@ThermalConductivity,@Density_SpecificHeat);
 
 plot(t,-ThermalConductivity(TC1).*H1/10^3,'LineWidth',1.4),axis([0 6600 0 120])
 xlabel('Time: t [s]','FontSize',14)
 ylabel('Heat flux: -\kappa T_x(L_1,t) [kW/m^2]','FontSize',14)

 
 % Solve the inverse problem in the interval [0,L1] and compute the surface
 % temperature. 

 lambda=1e-19; % First use too little regularization
 
 [T,Tx,x1]=SHESolver( x(1:ind(1)) , t , TC1 , H1 ,@ThermalConductivity,@Density_SpecificHeat,'SSP' ,lambda);
 Fode=T(end,:)';
 Hode=Tx(end,:);
 plot(t,Fode,'k-','LineWidth',1.4);
 xlabel('Time: t','FontSize',14);
 ylabel('Solutions: T_\lambda(0,t)','FontSize',14)
 %print -depsc Application-Solution-T-lambda-1e-19.eps
 
 plot(t,-ThermalConductivity(Fode)'.*Hode/10^3,'LineWidth',1.4),axis([0 6600 0 130])
 xlabel('Time: t [s]','FontSize',14)
 ylabel('Heat flux: -\kappa T_x(0,t) [kW/m^2]','FontSize',14)
 %print -depsc Application-Solution-Tx-lambda-1e-19.eps
 
 lambda=1e-7; % Second use too much regularization
 
 [T,Tx,x1]=SHESolver( x(1:ind(1)) , t , TC1 , H1 ,@ThermalConductivity,@Density_SpecificHeat,'SSP' ,lambda);
 Fode=T(end,:)';
 Hode=Tx(end,:);
 plot(t,Fode,'k-','LineWidth',1.4);
 xlabel('Time: t','FontSize',14);
 ylabel('Solutions: T_\lambda(0,t)','FontSize',14)
 %print -depsc Application-Solution-T-lambda-1e-7.eps
 
 plot(t,-ThermalConductivity(Fode)'.*Hode/10^3,'LineWidth',1.4),axis([0 6600 0 130])
 xlabel('Time: t [s]','FontSize',14)
 ylabel('Heat flux: -\kappa T_x(0,t) [kW/m^2]','FontSize',14)
 %print -depsc Application-Solution-Tx-lambda-1e-7.eps
 
 
 lambda=1e-13; % Third use just enough regularization
 
 [T,Tx,x1]=SHESolver( x(1:ind(1)) , t , TC1 , H1 ,@ThermalConductivity,@Density_SpecificHeat,'SSP' ,lambda);
 Fode=T(end,:)';
 Hode=Tx(end,:);
 plot(t,Fode,'k-','LineWidth',1.4);
 xlabel('Time: t','FontSize',14);
 ylabel('Solutions: T_\lambda(0,t)','FontSize',14)
 %print -depsc Application-Solution-T-lambda-1e-13.eps
 
 plot(t,-ThermalConductivity(Fode)'.*Hode/10^3,'LineWidth',1.4),axis([0 6600 0 130])
 xlabel('Time: t [s]','FontSize',14)
 ylabel('Heat flux: -\kappa T_x(0,t) [kW/m^2]','FontSize',14)
 %print -depsc Application-Solution-Tx-lambda-1e-13.eps
 

