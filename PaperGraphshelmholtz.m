%==========================================================================
%
% Test0: Set parameters and illustrate the testproblem to be used in the 
%  rest of the experiments section. We create Dirichlet data and also
%  compute both the exact solution U and obtain Dirichletdata F1,F2, and
%  Neumann data H1,H2. We use L=0.2 and k^2=12 for this test.
%
 close all,clear,clear global
 
 % Set parameters for the problem and create the grid.
 L=0.2; N=501; 
 dx=1/(N-1);x=dx*(0:N-1)';
 M=round(L/dx)+1;y=dx*(0:M-1)'; 
 [X,Y]=meshgrid(x,y);
 k2=12;
 [F1,H1,F2,H2,U]=AnalyticSolution( x , y , L , k2 ); 

 % Produce graphs
 mesh(X(1:6:end,1:6:end),Y(1:6:end,1:6:end),U(1:6:end,1:6:end));
 xlabel('x','FontSize',14);
 ylabel('y','FontSize',14);
 zlabel('Solution: u(x,y)','FontSize',14); 
 view([320 40]);
 %print -depsc F1-Exact-Mesh-1.eps
 
 plot(x,F2,'b--',x,H1,'k-','LineWidth',1.4);
 xlabel('x','FontSize',14);
 ylabel('Data: u(x,0) and u_y(x,a)','FontSize',14)
 print -depsc F1-Exact-Cauchy-Data.eps
 
%==========================================================================
%
% Test1: Solve the test problem using noisy data and illustrate the
% dependence on \lambda. Also compute en L-curve and plot the error
% as a function of \lambda.
%
 close all,clear,clear global
 
 % Set parameters for the problem and create the grid.
 L=0.2; N=501;k2=12;
 dx=1/(N-1);x=dx*(0:N-1)';
 M=round(L/dx)+1;y=dx*(0:M-1)'; 
 [X,Y]=meshgrid(x,y);
 [F1,H1,F2,H2,U]=AnalyticSolution( x , y , L , k2 );

 % Add noise for a prescribed level and solve the inverse problem

 NoiseLevel=0.001;
 H1e=H1+randn(size(H1))*NoiseLevel;
 F1e=F1+randn(size(F1))*NoiseLevel;
 
 % Compute the Solutions, Errors, and Residuals as a function of \lambda
 
  lambda = 10.^-(5:0.05:7);
  Errors = zeros(size(lambda));
  Residuals = zeros(size(lambda));
  SolNorm = zeros(size(lambda));
  
  for i = 1:length(lambda)
       [U]=HelmholtzCauchy(x , y , F1e , H1e , k2 ,'SSP',lambda(i));
        Fout=U(end,:)';
        Errors(i) = norm(Fout -F2)/sqrt(N);
        SolNorm(i)=norm(Fout(3:end)- 2*Fout(2:end-1) + Fout(1:end-2))/sqrt(N);
        %[Htmp]=SSPDeriv(H1e,lambda(i),0);
        [U,outVals] = DirectHelmholtzSolv(L,k2,[1,0],[H1e,Fout],[0,0]);
        F3=outVals(:,1);  
        Residuals(i) = norm(F1e -F3)/sqrt(N);
  end


 % Plot an L-curve  
 loglog(Residuals, SolNorm,'LineWidth',1.4);
 xlabel('Residual norm:||v_\lambda^\delta-g^\delta||_2','FontSize',14);
 ylabel('Solution norm :|| u_\lambda^\delta||_2','FontSize',14); 
 axis([0.0024    0.0201    0.0018    0.0056])
 hold on,k=33;plot(Residuals(k),SolNorm(k),'ro'),hold off
 %print -depsc F2-L-curve.eps

 % Plot the error as a function of \lambda

 semilogx(lambda,Errors,'LineWidth',1.4);
 xlabel('Regularization parameter value: \lambda','FontSize',14);
 ylabel('Error: ||u(x,0)- u_\lambda^\delta(x,a)||_2','FontSize',14)
 hold on,k=24;plot(lambda(k),Errors(k),'ro'),hold off
 %print -depsc F2-Error-vs-Lambda.eps
 
 % Plot the solutions for the different values of \lambda

[U]=HelmholtzCauchy(x , y , F1e , H1e , k2 ,'SSP',lambda(33)); % L-curve value
Fout=U(end,:)';
plot(x,F2,'b--',x,Fout,'k-','LineWidth',1.4);
xlabel('x','FontSize',14);
ylabel('Solutions: u(x,a) and u_\lambda^\delta(x,a)','FontSize',14)
print -depsc F2-Solution-L-Lambda.eps


[U]=HelmholtzCauchy(x , y , F1e , H1e , k2 ,'SSP',lambda(24)); % Minimum error value
Fout=U(end,:)';
plot(x,F2,'b--',x,Fout,'k-','LineWidth',1.4);
xlabel('x','FontSize',14);
ylabel('Solutions: u(x,a) and u^\delta_\lambda(x,a)','FontSize',14)
print -depsc F2-Solution-Min-Error.eps


tmp_lambda= 2e-8; 
[U]=HelmholtzCauchy(x , y , F1e , H1e , k2 ,'SSP',tmp_lambda); % too noisy value
Fout=U(end,:)';
plot(x,F2,'b--',x,Fout,'k-','LineWidth',1.4);
xlabel('x','FontSize',14);
ylabel('Solutions: u(x,a) and  u_\lambda^\delta(x,a)','FontSize',14)
print -depsc F2-Solution-Too-Little.eps

tmp_lambda=1e-4;
[U]=HelmholtzCauchy(x , y , F1e , H1e , k2 ,'SSP',tmp_lambda); % too smooth value
Fout=U(end,:)';
plot(x,F2,'b--',x,Fout,'k-','LineWidth',1.4);
xlabel('x','FontSize',14);
ylabel('Solutions: u(x,a) and u_\lambda^\delta(x,a)','FontSize',14)
print -depsc F2-Solution-Too-much.eps
 

 
%==========================================================================
%
% Test2: Find the optimal \lambda as a function of the noise level.
%
 close all,clear,clear global


  NoiseLevel=10.^-(1:0.1:4);
  Optimallambda = zeros(size(NoiseLevel));
  OptimalError=zeros(size(NoiseLevel));
   
   for i = 1:length(NoiseLevel)
              
        % Add the current level of noise to the data
         We=H1+randn(size(H1))*NoiseLevel(i);
         Ge=F1+randn(size(F1))*NoiseLevel(i);
  
        % Set the range of lambda values to test for 
        lambda= 10.^-(5:0.02:7);
        Error = zeros(size(lambda));
        
        % Compute error as a function of lambda
        for j = 1:length(lambda)
          [Uout]=HelmholtzCauchy(x , y , Ge , We , k2,'SSP',lambda(j));
          Fout=Uout(end,:)';
          Error(j) = norm(Fout -F2)/sqrt(N);
        end
        
        % Find the optimal lambda for this noiselevel
        [m,k]=min(Error);
        Optimallambda(i)=lambda(k(1));
        OptimalError(i)=m;
        fprintf('Testing noise level %e. Optimal \lambda=%f\n',NoiseLevel(i),Optimallambda(i))
   end
 
loglog(NoiseLevel,Optimallambda,'LineWidth',1.4);
xlabel('Noiselevel: \delta','FontSize',14);
ylabel('Optimal \lambda','FontSize',14);
%print -depsc F3-Optimal-lambda.eps


loglog(NoiseLevel,OptimalError,'LineWidth',1.4);
xlabel('Noiselevel: \delta','FontSize',14);
ylabel('Optimal error ||u(x,0)-u ^\delta_\lambda(x,a)||_2','FontSize',14);
%print -depsc F3-Optimal-Error.eps

 
 
%==========================================================================
%
% Test3: Find the optimal \lambda as a function of the noise level.
%
 close all,clear,clear global


