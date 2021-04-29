%
% Test with the modified Helmholtz equation where the wavenumber is a
% function.
%
%
 close all,clear,clear global
 
 % Set parameters for the problem and create the grid.
 L=0.2; N=501; 
 dx=1/(N-1);x=dx*(0:N-1)';
 M=round(L/dx)+1;y=dx*(0:M-1)'; 
 [X,Y]=meshgrid(x,y);
 k2=@(x,y)7+sin(x/3)*cos(y/2)+3*x*y^2;
 
 [F1,~,F2]=AnalyticSolution( x , y , L , k2(1/2,L/2)); 
 F1=(sin(pi*x).^4).*F1;F2=(sin(pi*x).^4).*F2;
 
 [U,outVals] = DirectHelmholtzSolv(L,k2,[0,0],[F1,F2],[1,1]);
 H1=outVals(:,1);H2=outVals(:,2);
 
 % Produce graphs
 mesh(X(1:4:end,1:4:end),Y(1:4:end,1:4:end),U(1:4:end,1:4:end));
 xlabel('x','FontSize',14);
 ylabel('y','FontSize',14);
 zlabel('Solution: u(x,y)','FontSize',14); 
 view([320 40]);
 print -depsc F8-Exact-Mesh-L022.eps 
 
 plot(x,H1,'k',x,F1,'b--','LineWidth',1.4);
 xlabel('x','FontSize',14);
 ylabel('Data: g(x)=u(x,0) and \eta(x)=u_y(x,0)','FontSize',14)
 print -depsc F8-Exact-Cauchy-Data.eps
 

 % Add noise for a prescribed level and solve the inverse problem

 NoiseLevel=0.001;
 H1e=H1+randn(size(H1))*NoiseLevel;
 F1e=F1+randn(size(F1))*NoiseLevel;


 % Compute the Solutions, Errors, and Residuals as a function of \lambda
 
  lambda = 10.^-(4:0.05:8);
  Errors = zeros(size(lambda));
  Residuals = zeros(size(lambda));
  SolNorm = zeros(size(lambda));
  
  for i = 1:length(lambda)
       [U]=HelmholtzCauchy(x , y , F1e,H1e, k2,'SSP',lambda(i));
        Fout=U(end,:)';
        Errors(i) = norm((Fout-F2))/sqrt(N);
        SolNorm(i)=norm(Fout(3:end)- 2*Fout(2:end-1) + Fout(1:end-2))/sqrt(N);
  
        [U,outVals] = DirectHelmholtzSolv(L,k2,[1,0],[H1e,Fout],[0,0]);
        F3=outVals(:,1);  
        Residuals(i) = norm(F1e -F3)/sqrt(N);
  end
  

 % Plot the error as a function of \lambda

 semilogx(lambda,Errors,'LineWidth',1.4);
 xlabel('Regularization parameter value: \lambda','FontSize',14);
 ylabel('Error: ||u(x,a)-u_\lambda^\delta(x,a)||_2','FontSize',14)
 axis([5*10^-8 10^-5 0.05 0.6])
 hold on,k=55;plot(lambda(k),Errors(k),'ro'),hold off
 %print -depsc F4-Errors-vs-Lambda.eps

% Plot an L-curve  
 loglog(Residuals, SolNorm,'LineWidth',1.4);
 xlabel('Residual norm:||v_\lambda^\delta-g^\delta||_2','FontSize',14);
 ylabel('Solution norm :|| u_\lambda^\delta||_2','FontSize',14); 
 %axis([5*10^-8 10^-5 0.05 0.6])
 hold on,k=53;plot(Residuals(k),SolNorm(k),'ro'),hold off
 %print -depsc F8-L-curve.eps


%=========================================================================%
%Uses the optimal lambda of error curve, computes the solution and        %
% the solution Mesh.                                                      %
%=========================================================================%
 
lambda=2e-7;
[U]=HelmholtzCauchy(x , y , F1e , H1e , k2,'SSP',lambda); % Minimum error value
Fout=U(end,:)';
plot(x,F2,'b--',x,Fout,'k','LineWidth',1.4);
xlabel('x','FontSize',14);
ylabel('Solutions: u(x,a) and u^\delta_\lambda(x,a)','FontSize',14)
print -depsc F4-Solution-error-curve-lambda.eps

mesh(X(1:6:end,1:6:end),Y(1:6:end,1:6:end),U(1:6:end,1:6:end));
 xlabel('x','FontSize',14);
 ylabel('y','FontSize',14);
 zlabel('Regularized Solution: u_\lambda^\delta(x,a)','FontSize',14); 
 view([320 40]);
 %print -depsc F8-Regularized-Mesh-Errors-curve.eps
 
 %========================================================================%
%                                                                         %
%Computes the L-curve value of lambda solution and Display its  mesh.     %
%=========================================================================%
 
lambda=2.5119e-7;
[U]=HelmholtzCauchy(x , y , F1e , H1e , k2,'SSP',lambda); % Minimum error value
Fout=U(end,:)';
plot(x,F2,'b--',x,Fout,'k-','LineWidth',1.4);
xlabel('x','FontSize',14);
ylabel('Solutions: u(x,a) and u^\delta_\lambda(x,a)','FontSize',14)
print -depsc F8-Solution-L-lambda.eps


mesh(X(1:6:end,1:6:end),Y(1:6:end,1:6:end),U(1:6:end,1:6:end));
 xlabel('x','FontSize',14);
 ylabel('y','FontSize',14);
 zlabel('Regularized Solution: u_\lambda^\delta(x,a)','FontSize',14); 
 view([320 40]);
 print -depsc F8-Regularized-Mesh-Lcurve.eps

%=========================================================================%
% Computes the solution with too-small lambda and display the mesh for the%                                                                        %
%corresponding solution                                                   %
%=========================================================================% 

lambda=5e-8;
[U]=HelmholtzCauchy(x , y , F1e , H1e , k2,'SSP', lambda); % L-curve value
Fout=U(end,:)';
plot(x,F2,'b--',x,Fout,'k-','LineWidth',1.4);
xlabel('x','FontSize',14);
ylabel('Solutions: u(x,a) and u_\lambda^\delta(x,a)','FontSize',14)
print -depsc F8-Solution-Too-Little.eps

 
 mesh(X(1:6:end,1:6:end),Y(1:6:end,1:6:end),U(1:6:end,1:6:end));
 xlabel('x','FontSize',14);
 ylabel('y','FontSize',14);
 zlabel('Regularized Solution: u_\lambda^\delta(x,a)','FontSize',14); 
 view([320 40]);
 print -depsc F8-Too-little-regularized-Mesh.eps


%=========================================================================%
% Computes the solution with too-large lambda and display the mesh for the%                                                                        %
%corresponding solution                                                   %
%=========================================================================% 

lambda=8e-6;
[U]=HelmholtzCauchy(x , y , F1e , H1e , k2,'SSP', lambda); % Minimum error value
Fout=U(end,:)';
plot(x,F2,'b--',x,Fout,'k-','LineWidth',1.4);
xlabel('x','FontSize',14);
ylabel('Solutions: u(x,a) and u^\delta_\lambda(x,a)','FontSize',14)
print -depsc F8-Solution-Too-much.eps

mesh(X(1:6:end,1:6:end),Y(1:6:end,1:6:end),U(1:6:end,1:6:end));
xlabel('x','FontSize',14);
ylabel('y','FontSize',14);
zlabel('Regularized Solution: u_\lambda^\delta(x,a)','FontSize',14); 
view([320 40]);
print -depsc F8-Too-much-regularized-Mesh.eps

%=========================================================================%
 %
 % Find the optimum lambda as a function of the noise level.
 %
 %=========================================================================%

 %NoiseLevel=10.^-(-1:0.06:3);
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
       % [U]=HelmholtzCauchy(x , y , F1e , H1e , k2 ,lambda(32));    
       [Uout]=HelmholtzCauchy(x , y , Ge , We , k2,'SSP',lambda(j));
        Fout=Uout(end,:)';
        Error(j) = norm(Fout -F2)/sqrt(N);
        fprintf(1,'The error for lambda=%e is %e\n',lambda,Error)
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
print -depsc F8-Optimal-lambda.eps


loglog(NoiseLevel,OptimalError,'LineWidth',1.4);
xlabel('Noiselevel: \delta','FontSize',14);
ylabel('Optimal error ||u(x,0)-u ^\delta_\lambda(x,a)||_2','FontSize',14);
print -depsc F8-Optimal-Error.eps

 


