function [x_new,f_vec] = x_update(x, mu_r, delta_r, b, sigma_e, Sr, lambda, usechol)
%This function updates the current image x by solving a TV-regularized
%least squares problem taking into account the model error

%INPUT:
%x: Current deblurred image (R^(m x n))
%mu_r: Current estimate of radius (R)
%delta_r: Standard deviation of radius (R)
%b: Original noisy blurred image (R^(m x n))
%Sr: Number of model error samples (Positive Integer)
%lambda: Regularization parameter (Positive Real Number)
%usechol: Boolean variable for chol usage (1 if yes 0 if no)

%OUTPUT:
%x_new: Updated deblurred image (R^(m x n))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialize
[n, m] = size(x);
N = n * m;
eta_samples = zeros(N, Sr);

%Precompute forward operation once for efficiency
A_mu_r_x = convb(x, mu_r);
    
for i = 1:Sr
    %Sample radius
    ri = normrnd(mu_r, delta_r);
    while ri <= 0
        ri = normrnd(mu_r, delta_r);
    end
    
    %Compute model error sample
    eta_i = convb(x, ri) - A_mu_r_x;
    eta_samples(:,i) = eta_i(:);
end

%Compute model error mean
mu_eta = mean(eta_samples, 2);
    
%NAIVE IMPLEMENTATION OF CHOLESKY FACTORIZATION
%C_eta = cov(eta_samples');
%C_total = C_eta + sigma_e^2 * eye(N);
%C_inv = C_total\eye(N);
%L = chol(C_inv);
%b_tilde = L * (b(:) - mu_eta);
%b_tilde = reshape(b_tilde, size(b));
%x_new = tv_weighted_deblurring_chol(b_tilde,L,mu_r,lambda);   
    
%Update x
%x_new = tv_weighted_deblurring_chol_generator(b_tilde,Yt,Zt,c,mu_r,lambda);
%x_new = tv_weighted_deblurring_chol(b_tilde, L, mu_r, lambda);

if usechol == 1 && Sr ~= 0

    %Generator representation of inverse cholesky factorization
    U = (eta_samples - mu_eta)/sqrt(Sr-1);
    Ut = U';
    [Wt,c] = egrss_potrf(Ut,Ut,sigma_e^2);
    [Yt,Zt] = egrss_trtri(Ut,Wt,c);
    b_tilde = egrss_trmv(Yt,Zt,1./c,b(:) - mu_eta);
    b_tilde = reshape(b_tilde,size(b));
    if nargout>1
        [x_new,f_vec] = FISTA_TVsmooth_CholGenerator(mu_r,Yt,Zt,c,b_tilde,lambda,x);
    else
        x_new = FISTA_TVsmooth_CholGenerator(mu_r,Yt,Zt,c,b_tilde,lambda,x);
    end
else
    %x_new = tv_weighted_deblurring(b, mu_r, lambda);
    if nargout>1
        [x_new,f_vec] = FISTA_TVsmooth(mu_r,b,lambda,x);
    else
        x_new = FISTA_TVsmooth(mu_r,b,lambda,x);
    end
end
end