% Generate test data to a folder
clear
clc
close all

% Set parameters.
true_radius = 5:5:100;
true_noise_std = 0.03;

% Read in test image.
im = im2double(imread('data/test4.jpg'));
im = im(:,:,1);

for i = true_radius
    % Blur with different levels and add noise
    im_blurred =  convb(im, i);
    b = im_blurred + randn(size(im_blurred))*true_noise_std;

    figure; 
    imagesc(b); 
    title('Blurred with noise'); 
    h = colorbar; 
    h.Limits = [0 1];
    colormap('gray');
    drawnow
    
    % save image to folder
    imwrite(b,['pipeline_test_data_medium/blurred_test_' num2str(i) '.tif'])
end