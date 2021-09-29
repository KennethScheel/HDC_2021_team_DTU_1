%% load competition data
clear, clc, close all;
addpath('egrssMatlab')

steps = 0:19;
n = 14;
pixel_r_est = zeros(n,1);

for i = 1:length(steps)
    step = num2str(steps(i));
    sample = '001'; % '001' to '100'
    font = 'times';   
    folder = ['C:\Users\Rainb\OneDrive\Documents\Skole\DTU\HDC 2021 - Image ' ...
        'deblurring project course\HDC 2021 Data\HDC2021_step' step '\step' step '\' font '\'];
%     if strcmp(font,'verdana')
%         bfile = ['CAM02\focusStep_' step '_' font 'Ref_size_30_sample_0' sample '.tif']; % blurred 
%         xfile = ['CAM01\focusStep_' step '_' font 'Ref_size_30_sample_0' sample '.tif']; % exact
%     elseif strcmp(font,'times')
%         bfile = ['CAM02\focusStep_' step '_' font 'R_size_30_sample_0' sample '.tif']; % blurred 
%         xfile = ['CAM01\focusStep_' step '_' font 'R_size_30_sample_0' sample '.tif']; % exact
%     end
%     lsfXfile = ['CAM02\focusStep_' step '_LSF_X.tif']; % lsf x
%     lsfYfile = ['CAM02\focusStep_' step '_LSF_Y.tif']; % lsf y
    psffile = ['CAM02\focusStep_' step '_PSF.tif']; % psf

%     b = im2double(imread([folder bfile]));
%     x = im2double(imread([folder xfile]));
%     lsfx = im2double(imread([folder lsfXfile]));
%     lsfy = im2double(imread([folder lsfYfile]));
    psf = im2double(imread([folder psffile]));

    % We first cut out the center region
    mid = floor(size(psf)/2);
    width = 150;
    Cpsf = psf(mid(1)-width:mid(1)+width, mid(2)-width:mid(2)+width);

    % pixel ratio from original - 200x300
    ratio = mean(size(psf)./[200,300]);

    % differences
    psfE = imbilatfilt(Cpsf,1,15);
    psfT = imbinarize(psfE)<1;
    pixel_r_est(i) = floor(min([max(sum(psfT,1)),max(sum(psfT,2))])/2);
    fprintf('step %.0f: \nestimated radius: %.0f',steps(i),pixel_r_est(i))
    fprintf('\n\n')
    
    figure(i);
    subplot(131)
    imagesc(Cpsf); colormap gray; axis image; 
    title('Original PSF','FontSize',18,'interpret','latex')
    subplot(132)
    imagesc(psfE); colormap gray; axis image; 
    title('Edge enhanced PSF','FontSize',18,'interpret','latex')
    subplot(133)
    imagesc(psfT); colormap gray; axis image; 
    title('Thresholded PSF','FontSize',18,'interpret','latex')
end
%%
save('all_r_est.mat','pixel_r_est')
%%
clc, clear
load all_r_est
figure;
steps = 0:19;
plot(steps,all_r_est,'o','linewidth',2)

