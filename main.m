function main(input_folder,output_folder,step)

% used for testing
%clear; close all force; clc;
%step = 15;
%input_folder = ['competition_data_single_sample/step' num2str(step)];    % function input
%output_folder = ['competition_data_single_output/step' num2str(step)];    % function input

% Add package
addpath('egrssMatlab')

% Options
save_deblur = 1;    %Save output deblurred image?
save_x = 0;         %Save output deblurred image as .mat file (used in testing)
use_egrss = 1;      %Use egrss package for r_update? If 0 only works on small-scale.
use_gpu = 1;        %Use gpu for faster computations? Requires Parallel computing toolbox

% pre-estimated r (we use these as initial guesses)
r0 = [1, 4,12,18,25,37,43,49,61,66,74,85,90,96,101,107,114,120,126,130];
dr0 = repmat(0.3,20,1);

% binary flag for which steps we wish to use the iterative radius estimation
% with x_update and r_update 
estimate_r = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1] ;  

% Other Parameters
%Regularization parameter for final deblur at each step
lambdas_deblur = ...
    [1e-3,1e-2,1e-2,1e-2,1e-2,1e-2,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3, ...
    1e-3,1e-3,1e-3,1e-3,1e-3,1e-3,1e-3]; 
%Regularization parameters for r_update patch at each step
lambdas_patch = [1,1,1,1,1,1,1,1,1,1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1]; 

lambda_deblur = lambdas_deblur(step+1);    %Regularization parameter for final deblur
lambda_patch = lambdas_patch(step+1);      %Regularization parameter for r_update patch

mu_r0 = r0(step+1);         %Initial radius 
delta_r0 = dr0(step+1);     %Initial variance

% CHANGE THIS IS IF THE ITERATIVE RADIUS ESTIMATION IS TOO SLOW
speed = 'medium';
% speed of algorithm dependent on patch size
switch speed 
    case 'fast'        
        patch_width  = 300; % Width of patch for radius estimation
        patch_height = 150; % Heightof patch for radius estimation
    case 'medium'
        patch_width  = 600; % Width of patch for radius estimation
        patch_height = 280; % Heightof patch for radius estimation
    case 'slow'
        patch_width  = 1800; % Width of patch for radius estimation
        patch_height = 400; % Heightof patch for radius estimation
end

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
    if estimate_r(step+1)       % only do r_estimation at certain steps
        % ==== Prepare patches =====
        mid = floor(size(b)/2)+mid_shift;
        hpatch_width = patch_width/2;
        hpatch_height = patch_height/2;
        b_patch     = b(mid(1)-hpatch_height:mid(1)+hpatch_height, mid(2)-hpatch_width:mid(2)+hpatch_width);
        x = zeros(size(b_patch));
        %figure(1); imshow(b_patch); title('Blurred image (patch)'); drawnow;
    
        for k = 1:n_iter
            % Update x
            x = x_update(x, mu_r, delta_r, b_patch, sigma_e, Sx, lambda_patch, 1);
            %figure(2); imshow(x); title('Current deblurred patch'); drawnow;
            
            % Filter x
            x = medfilt2(x, [5,5]);         % median filter to remove noise from regularization
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
      
    % ==== Deblur with radius estimate ====
    x_estimate = x_update(x, mu_r, delta_r, b, sigma_e, 0, lambda_deblur, 0);
    if use_gpu == 1
        x_estimate = gather(x_estimate);
    end
    
    % use imsharpen to improve contrast in image and filter out some noise
    % parameters have been chosen to get the best OCR score
    x_estimate = imsharpen(x_estimate,'Radius',60,'Amount',5);
    
    %figure(5); imshow(x_estimate); title("Deblur with estimate"); drawnow;
   
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
