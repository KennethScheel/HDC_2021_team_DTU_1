%function main_pipeline(input_folder,output_folder,step)
input_folder = 'pipeline_test_data_medium';    % function input
output_folder = 'pipeline_output_data_medium'; % function input

addpath('egrssMatlab')

profile on

% Set parameters.
mu_r = 10:10:200;        % initial mean radius estimates (true ones here)
delta_r = 0.3;
Sr = 100;
lambda_tv = 0.01;
alpha = 0.5;
use_chol = 0;
K = 2;

% Get list of all .tif files in the directory
imagefiles = dir([input_folder '/*.tif']);      
nfiles = length(imagefiles);    % Number of files found

for i = 8
    currentfilename = imagefiles(i).name;
    re = 1;
    b = imresize(im2double(imread([imagefiles(i).folder '\' currentfilename])),re);
    x = zeros(size(b));   % initial guess just zeros
        
    figure;
    imagesc(b); 
    title('Blurred with noise'); 
    h = colorbar; 
    h.Limits = [0 1];
    colormap('gray');
    drawnow
    
    sigma_e = std2(b(1:20,1:20)); % estimate noise std from small corner patch
    % radius
    if length(currentfilename) == 19
        mu_r = str2double(currentfilename(14:15))+1;
    else
        mu_r = str2double(currentfilename(14:16))+1;
    end
    
    disp(['it: ', num2str(0)])
    disp(['  mu_r: ', num2str(mu_r)])
    disp(['  delta_r: ', num2str(delta_r)])
    tic
    % main deblurring loop
    for k = 1:K
        % Update x
        tic
        [x,~] = x_update(x, mu_r, delta_r, b, sigma_e, Sr, lambda_tv, use_chol);
        toc
        
        x = medfilt2(x);                % median filter to remove noise from regularization
        x = imbilatfilt(x,2*sigma_e,4); % edge enhancement

        % Update psf radius estimate
        tic
        [mu_r, delta_r] = r_update(x, b, mu_r, delta_r, sigma_e, Sr, alpha);
        toc
        
        disp(['it: ', num2str(k)])
        disp(['  mu_r: ', num2str(mu_r)])
        disp(['  delta_r: ', num2str(delta_r)])
    end
    toc

    figure;
    imagesc(x); 
    title('Deblurred image'); 
    h = colorbar; 
    h.Limits = [0 1];
    colormap('gray');
    drawnow
    
    % saves image to output folder
%    output_file = [output_folder '/' currentfilename(1:end-4) '.png'];
%    imwrite(x,output_file)
end

%profile viewer

%end

