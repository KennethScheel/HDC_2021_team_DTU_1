function [x,f_vec] = FISTA_TVsmooth(r,b,lambda,x0)

%This function solves a smooth approximation to the TV regularized least squares problem

%   min_x 1/2||Ax-b||_2^2+alpha TV_{rho}(x)
%   s.t. 0<=x<=1

%where the regularization term TV_{rho} is a smooth approximation to the TV
%of x. The problem is solved using accelerated projected gradient (FISTA).

%INPUT:
%r: Radius parameter for blurring kernel (Positive real)
%b: Noisy blurry image (R^(m x n))
%lambda: Regularization paramter (Positive real)
%x0: Initial guess for deblurred image (R^(m x n))

%OUTPUT:
%x: Deblurred Image (R^(m x n))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialize
[m,n] = size(x0);           %Size of image
rho = 10^(-2);              %Smoothing parameter
h = 1;                      %Fineness of derivative discretization
epsilon = 10^(-6);          %Stopping parameter
maxiters = 100;             %Maximum number of iterations
N = m*n;
b = reshape(b,N,1);       
x0 = reshape(x0,N,1);

k = 0;
x = x0;
y = x;
xold = x;
converged = 0;
told = 1;


%Power iteration to determine Lipschitz constant
n_lip = 10;
z = randn(N,1);
for i=1:n_lip
    %Multiply with A^TL^TLA
    z1 = reshape(convb(reshape(z,m,n),r),N,1);
    z1 = reshape(convb(reshape(z1,m,n),r),N,1);
    z = z1/norm(z1);
end

%Compute Lipschitz constant
z1 = reshape(convb(reshape(z,m,n),r),N,1);
z1 = reshape(convb(reshape(z1,m,n),r),N,1);

L = z'*z1 + 8*lambda/rho/h^2;

%Compute Finite difference matrix with Neumann b.c.
Dfd_N = spdiags([-ones(n-1,1),ones(n-1,1);0 1]/h,[0,1],n,n);
Dfd_M = spdiags([-ones(m-1,1),ones(m-1,1);0 1]/h,[0,1],m,m);

D = vertcat(kron(speye(n), Dfd_M), kron(Dfd_N, speye(m)));
% Compute "smooth" TV and its gradient
phi = @(y) sqrt(sum(reshape(y,[],2).^2,2)+rho^2);
Dx = D*x;

if nargout>1
    smooth_tv = sum(phi(Dx));
    f_vec = 1/2*norm(reshape(convb(reshape(x,m,n),r),N,1)-b,2)^2+lambda*smooth_tv;
end


h = waitbar(0,'Deblurring image');

while ~converged && k<maxiters
    k = k+1;
    waitbar(k / maxiters)
    %Compute TV-gradient
    Dy = D*y;
    grad_tv = D'*(Dy.*repmat(1./phi(Dy),2,1));    
    
    %Compute full forward mapping
    Ay = reshape(convb(reshape(y,m,n),r),N,1);
    
    %Compute full gradient
    res = Ay - b;
    grad = reshape(convb(reshape(res,m,n),r),N,1) + lambda*grad_tv;
    
    %Take gradient step
    x = y - (1/L)*grad;
    
    %Projection
    x = min(max(0,x),1);
    
    %Momentum
    t = (1 + sqrt(1+4*told^2))/2; 
    y = x + (told-1)/t*(x-xold);
    
    if norm(x-xold)/norm(xold)<epsilon
        disp(['The algorithm stopped at iteration ' num2str(k) ' because the relative change in objective is below the set tolerance'])
        converged = 1;
    end
    
    xold = x;
    told = t;
    
    if nargout>1
        Dx = D*x;
        smooth_tv = sum(phi(Dx));
        f_val = 1/2*norm(reshape(convb(reshape(x,m,n),r),N,1)-b,2)^2+lambda*smooth_tv;
        f_vec = [f_vec f_val];
    end
end
close(h)
%Reshape to get image output
x = reshape(x,m,n);

if k==maxiters
    disp('The algorithm stopped because the number of iterations reached the set maximum')
end
end