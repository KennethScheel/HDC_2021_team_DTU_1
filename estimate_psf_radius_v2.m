%% load competition data
clear, clc, close all;
addpath('egrssMatlab')
addpath('PSF')

steps = 0:19;

% load exact PSF
step = '0';
psffile = ['focusStep_' step '_PSF.tif']; % psf
psf_0 = im2double(imread(psffile));

% cut out the center region
mid = floor(size(psf_0)/2);
width = 500;
Cpsf_0 = psf_0(mid(1)-width:mid(1)+width, mid(2)-width:mid(2)+width);
re = 1;
Cpsf_0 = imresize(Cpsf_0,re);
% using the noise-free PSF gives radius estimates close to each other, so i
% think we need to use the noisy one
%Cpsf_0 = imcomplement(Cpsf_0);    % noise-free is inverted apparently,
figure(1);
imagesc(Cpsf_0); colormap gray; axis image; 
title(['PSF step: ' num2str(0)],'FontSize',18,'interpret','latex')
drawnow

% parameters
K = 10;
mu_r = 10;
delta_r = 0.3;
alpha = 0.5;
Sr = 200;

mu_r_hist = zeros(length(steps),K+1); mu_r_hist(:,1) = mu_r*ones(length(steps),1);
delta_r_hist = zeros(length(steps),K+1); delta_r_hist(:,1) = delta_r*ones(length(steps),1);

for i = 6
    step = num2str(steps(i));
    psffile = ['focusStep_' step '_PSF.tif']; % psf
    psf = im2double(imread(psffile));
    
    sigma_e = std2(psf(1:20,1:20)); % estimate noise std from small corner patch

    Cpsf = psf(mid(1)-width:mid(1)+width, mid(2)-width:mid(2)+width);
    Cpsf = imresize(Cpsf,re);
    
    % initial guess
    mu_r = 5;
    delta_r = 0.3;

    disp(['step: ' num2str(i-1) ', it: ', num2str(0)])
    disp(['  mu_r: ', num2str(mu_r)])
    disp(['  delta_r: ', num2str(delta_r)])
    
    % estimate radius
    for k = 1:K
        % Update r
        [mu_r, delta_r] = r_update(Cpsf_0, Cpsf, mu_r, delta_r, sigma_e, Sr, alpha);

        disp(['step: ' num2str(i-1) ', it: ', num2str(k)])
        disp(['  mu_r: ', num2str(mu_r)])
        disp(['  delta_r: ', num2str(delta_r)])

        mu_r_hist(i,k+1) = mu_r;
        delta_r_hist(i,k+1) = delta_r;
    end
    
    figure(i);
    imagesc(Cpsf); colormap gray; axis image; 
    title(['PSF step: ' num2str(i-1)],'FontSize',18,'interpret','latex')
    drawnow
end

% PSF radius estimate is final iteration for each step
psf_r_est = mu_r_hist(:,end);

%%
load all_r_est
load all_r_est_v2

r1 = pixel_r_est;
r2 = psf_r_est;

plot([r1,r2])

steps = (0:19)';

r1n = interp1(steps(1:14), r1(1:14), steps, 'linear', 'extrap');
r2n = interp1(steps(1:9), r2(1:9), steps, 'linear', 'extrap');

plot(steps,[r1n,r2n],'linewidth',2)
hold on
scatter(steps,r1,30,'b','filled')
scatter(steps,r2,30,'r','filled')
legend('Extrapolated v1','Extrapolated v2','Estimates v1','Estimates v2','location','best')
title('PSF radius estimate at each step','fontsize',14)
xlabel('blurring step')
ylabel('PSF radius')

