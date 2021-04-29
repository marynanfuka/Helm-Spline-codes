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
% print -depsc F4-Solution-Mesh.eps
 
 plot(x,H1,'k',x,F1,'b--','LineWidth',1.4);
 xlabel('x','FontSize',14);
 ylabel('Data: g(x)=u(x,0) and \eta(x)=u_y(x,0)','FontSize',14)
 print -depsc F4-Cauchy-Data.eps
 

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
  end
  

 % Plot the error as a function of \lambda

 semilogx(lambda,Errors,'LineWidth',1.4);
 xlabel('Regularization parameter value: \lambda','FontSize',14);
 ylabel('Error: ||u(x,a)-u_\lambda^\delta(x,a)||_2','FontSize',14)
 axis([5*10^-8 10^-5 0.05 0.6])
 hold on,k=55;plot(lambda(k),Errors(k),'ro'),hold off
% print -depsc F4-Errors-vs-Lambda.eps

 % Plot the solutions for the different values of \lambda
 
 lambda=2e-7;
[U]=HelmholtzCauchy(x , y , F1e , H1e , k2,'SSP',lambda); % Minimum error value
Fout=U(end,:)';
plot(x,F2,'k-',x,Fout,'b--','LineWidth',1.4);
xlabel('x','FontSize',14);
ylabel('Solutions: u(x,a) and u^\delta_\lambda(x,a)','FontSize',14)
print -depsc F4-Solution-lambda-2e-7.eps

 
 lambda=5e-8;
[U]=HelmholtzCauchy(x , y , F1e , H1e , k2,'SSP',lambda); % Minimum error value
Fout=U(end,:)';
plot(x,F2,'k-',x,Fout,'b--','LineWidth',1.4);
xlabel('x','FontSize',14);
ylabel('Solutions: u(x,a) and u^\delta_\lambda(x,a)','FontSize',14)
print -depsc F4-Solution-lambda-5e-8.eps


 
 lambda=8e-6;
[U]=HelmholtzCauchy(x , y , F1e , H1e , k2,'SSP',lambda); % Minimum error value
Fout=U(end,:)';
plot(x,F2,'k-',x,Fout,'b--','LineWidth',1.4);
xlabel('x','FontSize',14);
ylabel('Solutions: u(x,a) and u^\delta_\lambda(x,a)','FontSize',14)
print -depsc F4-Solution-lambda-8e-6.eps

