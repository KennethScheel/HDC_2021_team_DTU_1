%function main(input_folder,output_folder,step)
clear; close all force; clc;
step = 2;
input_folder = ['competition_data_single_sample/step' num2str(step)];    % function input
output_folder = ['competition_data_single_output/step' num2str(step)];    % function input

% Add package
addpath('egrssMatlab')

% Options
save_deblur = 1;    %Save output deblurred image?
save_x = 1;         %Save estimate of x
use_egrss = 1;      % Use egrss package for r_update? If 0 only works on small-scale.
use_gpu = 1;        %Use gpu for faster computations?

% pre-estimated r (we use these as initial guesses)
% r0 = [1.0000, 5.8661, 15.2660, 26.2374, 40.1595, 52.2180, 60.2950, ...
%      69.2086, 78.9692, 83.3661, 88.4403, 91.3648, 96.6185, 96.1995, ...
%      99.4622, 112.4390, 116.3366, 115.5400, 132.7100, 133.6777];
% r0 = round(r0);
r0 = [1,4,15,19,26,35,42,50,58,66,74,84,90,96,101,107,114,120,126,130];
dr0 = repmat(0.3,20,1);

mu_r0 = r0(step+1);         %Initial radius
delta_r0 = dr0(step+1);     %Initial variance

% Other Parameters
lambda_deblur = 0.01;    %Regularization parameter for final deblur

% lambda_patch will need to be scaled if we use Sx>0. For Sx>0 the following was
% found: lambda = 1 for step > 5, smaller lambda for smaller steps
% (maybe 0.1, 0.01?) 

speed = 'fast';
% speed of algorithm dependent on patch size
switch speed 
    case 'fast'        
        patch_width  = 300; % Width of patch for radius estimation
        patch_height = 150; % Heightof patch for radius estimation
        lambda_patch = 1;     %Regularization parameter for r update
    case 'medium'
        patch_width  = 600; % Width of patch for radius estimation
        patch_height = 280; % Heightof patch for radius estimation
        lambda_patch = 1;     %Regularization parameter for r update
    case 'slow'
        patch_width  = 1800; % Width of patch for radius estimation
        patch_height = 400; % Heightof patch for radius estimation
        lambda_patch = 1;     %Regularization parameter for r update
end

estimate_r_at_step = 5; % step # where we start to use r_update algorithm
n_iter = 10;      % Number of iterations for r_update
Sr = 100;         % Number of samples for r_update
Sx = 100;         % Number of samples for x_update
alpha = 0.1;      % Relaxation parameter in varience est for r_update
mid_shift = 40;   %Shift center of image to better align letters?


% =============== Algorithm start ===============

% Get list of all .tif files in the directory
imagefiles = dir([input_folder '/*.tif']);
nfiles = length(imagefiles);    % Number of files found

% Loop over all images in folder and deblur
for i = 1:nfiles
    
    %Load current image
    currentfilename = imagefiles(i).name;
    b = im2double(imread([imagefiles(i).folder '\' currentfilename]));
    
    % Estimate noise standard deviation
    sigma_e = std2(b(1:50,1:50)); % estimate noise std from small corner patch
        
    % Initial guess
    mu_r = mu_r0;
    delta_r = delta_r0;
    [mu_r,delta_r,0]
    
    % ==== Iteration for r estimation =====
    if step >= estimate_r_at_step       % only do r_estimation after certain step
        % ==== Prepare patches =====
        mid = floor(size(b)/2)+mid_shift;
        hpatch_width = patch_width/2;
        hpatch_height = patch_height/2;
        b_patch     = b(mid(1)-hpatch_height:mid(1)+hpatch_height, mid(2)-hpatch_width:mid(2)+hpatch_width);
        x = zeros(size(b_patch));
        figure(1); imshow(b_patch); title('Blurred image (patch)'); drawnow;
    
        for k = 1:n_iter
            % Update x
            x = x_update(x, mu_r, delta_r, b_patch, sigma_e, Sx, lambda_patch, 1);
            figure(2); imshow(x); title('Current deblurred patch'); drawnow;

            x_old = x;
            % Filter x
            x = medfilt2(x, [5,5]);                % median filter to remove noise from regularization
            x = imbilatfilt(x,2*sigma_e,4); % edge enhancement
                       
            % Update r
            [mu_r, delta_r] = r_update(x, b_patch, mu_r, delta_r, sigma_e, Sr, alpha, use_egrss);

            % Show result
            [mu_r,delta_r,k]
        end
    end
    
    if use_gpu == 1
        x = gpuArray(zeros(size(b)));
        b = gpuArray(b);
    end
      
    % ==== Deblur with estimate ====
    x_estimate = x_update(x, mu_r, delta_r, b, sigma_e, 0, lambda_deblur, 0);
    if use_gpu == 1
        x_estimate = gather(x_estimate);
    end
    
    figure(5); imshow(x_estimate); title("Deblur with estimate"); drawnow;
   
    % Save to file
    if save_deblur == 1
        % saves image to output folder
        output_file = [output_folder '/' currentfilename(1:end-4) '.png'];
        imwrite(x_estimate,output_file)
    end
    if save_x == 1
        save([output_folder '/' currentfilename(1:end-4) '.mat'],'x_estimate')
    end
end