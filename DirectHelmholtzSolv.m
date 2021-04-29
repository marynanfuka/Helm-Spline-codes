%DirectHelmholtzSolv: 
%
% Explain the problem....(I will complete this tomorrow.)
%
% Usage:
%  >>[U,outVals] = DirectHelmholtzSolv(L,k2,bndType,bndVals,outType)
%
%
%
function[U,outVals,K] = DirectHelmholtzSolv(L,k2,bndType,bndVals,outType)


%
% The grid parameter dx is determined by the length of the boundary value vectors 
% and the assumption that the domain is 0<x<1 and 0<y<L. Might have to
% change L slightly to get dx=dy.
%

 [n1,m1]=size(bndVals);,if n1<m1,bndVals=bndVals';,end;
 [n,tmp] = size(bndVals);
 dx=1/(n-1);
 x=dx*(0:n-1)';
 m=round(L/dx)+1;
 y=dx*(0:m-1)';
 
 
%
% The wavenumber k2 can be either a constant or a function. If its a
% function we compute Lambda=k2(1/2,L/2) and use it to check if a matrix
% already exists.
%
 if isnumeric(k2),Lambda=k2;,else,Lambda=k2(1/2,L/2);,end 

% 
% Check if an appropriate matrix exists. Otherwise create the matrix and save it.
% 
 if  ~CheckMatrix(n,m,Lambda,bndType)
   % Need to create the matrix.
   [A]=CreateMatrix(n,m,k2,bndType);
   StoreMatrix(A,n,m,Lambda,bndType); 
 end;

%
% Always need to recompute the boudnary conditions.  
%      
 [b]=CreateRHS(n,m,bndType,bndVals); 
 
%
% We have the sparse LU decomposition already done in global workspace. 
%
 global MatrixPrecompLU
 L=MatrixPrecompLU.L{ bndType(1)+2*bndType(2)+1 };
 U=MatrixPrecompLU.U{ bndType(1)+2*bndType(2)+1 };
 P=MatrixPrecompLU.P{ bndType(1)+2*bndType(2)+1 };
 U=U\( L\(P*b) );
 U=reshape(U,n,m);


%
% Compute output at y=0
%
 if outType(1)==0, % Want u(x,0)
    outVals(:,1)=U(:,1)';
 elseif outType(1)==1,%Want u_x(x,0)    
    outVals(:,1)=(-3*U(:,1)+4*U(:,2)-U(:,3))'/2/dx; 
 else,
    error('Invalid choice for outType(1)');
 end;

%
% Compute output at y=L
%
 if outType(2)==0, % Want u(x,0)
    outVals(:,2)=U(:,m)';
 elseif outType(2)==1,%Want u_x(x,0)
    outVals(:,2)=(3*U(:,m)-4*U(:,m-1)+U(:,m-2))'/2/dx;
 else,
    error('Invalid choice for outType(2)');
 end;

if nargout > 2,
    %
    % Create a linear operator Kf=g, where K is defined by using the
    % boundary conditions Uy(x,0)=0 and U(x,L)=F. Only works properly
    % if Dirichlet-Neumann boundary conditions were picked.
    %
    A=MatrixPrecompLU.A{ bndType(1)+2*bndType(2)+1 };
    E=speye(n*m);
    K=A\E(:,(m-1)*n+1:m*n);
    K=full(K(1:n,:));
end;
 
U=U'; % This means the dimensions will be the same as the matrices created 
      % by meshgrid() to simplify plotting.
      
%=========================================================================%
%                                                                         %
%                     Local Subroutines                                   %
%                                                                         %
%=========================================================================%
%=========================================================================%
%
% Function that checks if the matrix stored in global memory exists and 
% is correct for the current problem.
%
%=========================================================================%
function [bool]=CheckMatrix(n,m,Lambda,bndType); 

  global MatrixPrecompLU

  bool=true; % Assume correct unless proven otherwise.
  
  if prod(size(MatrixPrecompLU))==0, 
     %
     % No Matrix exists in global memory
     %
      bool=false;
  
  elseif ~( ( MatrixPrecompLU.n == n ) & ( MatrixPrecompLU.m == m ) & (MatrixPrecompLU.Lambda == Lambda) ),
    % 
    % Matrix exists but isn't computed for the correct dimensions. Clear current matrix and 
    % recompute.
    %
     bool=false;
     
elseif prod(size(MatrixPrecompLU.A{ bndType(1)+2*bndType(2)+1}))==0,
    %
    % This particular combinations of boundary conditions isn't yet
    % treated.
    %
     bool=false;
end;
  
  
%=========================================================================%
%
% Function that stores the created matrix in a structure in the global  
% workspace. 
%=========================================================================%
function []=StoreMatrix(A,n,m,Lambda,bndType); 

  global MatrixPrecompLU
  
 %
 % Check if the data structure already exists and that the grid parameters
 % are the same. Othervise clear old and create a new one.
 %
  if prod( size(MatrixPrecompLU) )== 0  % Didn't exist in global workspace before cerated by global...
    MatrixPrecompLU.n=n;
    MatrixPrecompLU.m=m;
    MatrixPrecompLU.Lambda=Lambda;
    MatrixPrecompLU.A=cell(4,1);  % 4 different combinations of boundary values we can use.
    MatrixPrecompLU.L=cell(4,1);
    MatrixPrecompLU.U=cell(4,1);
    MatrixPrecompLU.P=cell(4,1);
  end;
  
  
  bool = ( MatrixPrecompLU.n == n ) & ( MatrixPrecompLU.m == m ) & (MatrixPrecompLU.Lambda == Lambda); 
  if ~bool,
    %
    % Existing data structure needs to be wiped. Create empty structure.
    %
     MatrixPrecompLU.n=n;
     MatrixPrecompLU.m=m;
     MatrixPrecompLU.Lambda=Lambda;
     MatrixPrecompLU.A=cell(4,1);  % 4 different combinations of boundary values we can use.
     MatrixPrecompLU.L=cell(4,1);
     MatrixPrecompLU.U=cell(4,1);
     MatrixPrecompLU.P=cell(4,1);
  end;
       
  % 
  % Now we compute the sparse LU factorization of A and store the result in
  % the structure. There are 4 different combinations of boundary
  % conditions 
  % 
  % 1 - Dirichlet - Dirichlet
  % 2 - Neumann - Dirichlet
  % 3 - Dirichlet - Neumann
  % 4 - Neumann - Neumann
  %
   MatrixPrecompLU.A{ bndType(1)+2*bndType(2)+1 }=A;
   [L,U,P]=lu(A);
   MatrixPrecompLU.L{ bndType(1)+2*bndType(2)+1 }=L;
   MatrixPrecompLU.U{ bndType(1)+2*bndType(2)+1 }=U;
   MatrixPrecompLU.P{ bndType(1)+2*bndType(2)+1 }=P;
  
  
      
%=========================================================================%
%
% Create the matrix that represent the differential operator on the grid.
% The matrix only depends on the grid, the types of boundary conditions and Lambda. 
%
%=========================================================================%
function [A]=CreateMatrix(n,m,k2,bndType); 
  if isnumeric(k2),Lambda=k2;,else,Lambda=k2(1/2,(m-1)/2/(n-1));,end 
  fprintf(1,'Creating Matrix for n=%i, Lambda=%f, bdnType=[%f,%f].\n',n,Lambda,bndType(1),bndType(2));
  A = sparse(n*m,n*m);
  dx = 1/(n-1);
 %
 % Set the boundary condition u(0,y)=0.
 %
  i = 1;
  for j = 1:m,
    ind = i+(j-1)*n;
    A(ind,ind) = 1;    
  end
 %
 % Les points intérieurs por lequels l'équation de Helmholtz est résolu.
 %
  for j = 2:m-1
    for i = 2:n-1 ;     
         ind = i +(j-1)*n; 
          % Check if the wavenumber is constant or a function
          if isnumeric(k2),Lambda=k2;,else,Lambda=k2((i-1)*dx,(j-1)*dx);,end 
         A(ind,[ind-n,ind-1,ind,ind+1,ind+n]) = [-1,-1,4-Lambda*dx^2,-1,-1];
    end
  end
 %
 % Set boundary condition u(1,y)=0
 %
  i = n;
  for j = 1:m,
     ind = i+(j-1)*n;
     A(ind,ind) = 1;
  end


  %
  % Points pour lequels on a les conditions initiales for (x,0) 0<x<1
  %
   j = 1;
   for i = 1:n-1
     ind = i+(j-1)*n;
     if bndType(1)==0,    % Dirichlet conditions for y=0.
       A(ind,ind) = 1;
      elseif bndType(1)==1,  % Neumann boundary conditions at y=0.
          A(ind,[ind,ind+n,ind+2*n]) = [-3,4,-1]/2;
      else      % Invalid type of boundary condition
         error('Wrong type of boundary condition at y=0');
      end;
   end
 %
 % Points pour lequels on a les conditions initiales u(b,y) = 0
 %
  j = m;
  for i= 1:n-1
     ind = i+(j-1)*n; 
    
     %
     % Set boundary conditios at y=L!
     %
      if bndType(2)==0,    % Dirichlet conditions for y=L.
         A(ind,ind) = 1;
      elseif bndType(2)==1,  % Neumann boundary conditions at y=L.
          A(ind,[ind,ind-n,ind-2*n]) = [3,-4,1]/2;  
      else      % Invalid type of boundary condition
         error('Wrong type of boundary condition at y=0');
      end;
  end
  
  
      
%=========================================================================%
%
% Create the Right-handside vector for the differential operator on the grid.
% The vector only depends on the grid, the types of boundary conditions. 
%
%=========================================================================%
function [b]=CreateRHS(n,m,bndType,bndVals); 
  dx = 1/(n-1);  
  b = zeros(n*m,1);
 %
 % Points pour lequels on a les conditions initiales for (x,0) 0<x<1
 %
  j = 1;
  for i = 1:n-1
     ind = i+(j-1)*n;
     if bndType(1)==0,    % Dirichlet conditions for y=0.
        b(ind,1) = bndVals(i,1);
     elseif bndType(1)==1,  % Neumann boundary conditions at y=0.
        b(ind,1) = dx*bndVals(i,1);
     else      % Invalid type of boundary condition
        error('Wrong type of boundary condition at y=0');
     end;
   end
 %
 % Points pour lequels on a les conditions initiales u(b,y) = 0
 %
  j = m;
  for i= 1:n-1
     ind = i+(j-1)*n;     
     %
     % Set boundary conditios at y=1!
     %
      if bndType(2)==0,    % Dirichlet conditions for y=0.
        b(ind,1) = bndVals(i,2);
      elseif bndType(2)==1,  % Neumann boundary conditions at y=0.
       b(ind,1) = dx*bndVals(i,2);
      else      % Invalid type of boundary condition
       error('Wrong type of boundary condition at y=0');
      end;   
  end


 %=========================================================================%

