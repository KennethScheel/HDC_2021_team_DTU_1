# Helsinki Deblur Challenge 2021 
## Team DTU - Group 1 authors
- Maria Knudsen
- Frederik Listov-Saabye Pedersen
- Kenneth Scheel

Supervisors

- Nicolai André Brogaard Riis
- Yiqiu Dong
- Per Christian Hansen

Department of Applied Mathematics and Computer Science, Technical University of Denmark, Kgs. Lyngby, Denmark

E-mail: `s171246@student.dtu.dk, s183984@student.dtu.dk, s174488@student.dtu.dk, nabr@dtu.dk, yido@dtu.dk, pcha@dtu.dk`

## About

This GitHub repository has been created as the submission for group 1 (out of 2) of Team DTU for the Helsinki Deblur Challenge 2021. The Matlab functions and scripts in the repository, implements an image deblurring algorithm with Point-Spread-Function (PSF) radius estimation. The algorithm is based on a paper written by our supervisors where they develop a CT-reconstruction method, which uses uncertainty quantification to estimate view angles using only measured CT data [Riis, N. A. B., Dong, Y. and Hansen, P. C.: Computed Tomography with View Angle
Estimation using Uncertainty Quantification. Inverse Problems, 37, pp. 065007, (2021)].

## Description of algorithm

The algorithms foundation lies in the assumption that a blurred image b consists of two things:
- the exact image x convolved (blurred) with a PSF of radius r
- a Gaussian noise term e, with mean 0 and some standard deviation sigma 

If we denote the "convolution with a given PSF of radius r" operator as A(r), we can write for any blurred image b

b = A(r)x + e 

The PSF radius r is the parameter which controls the level of blurring. The deblurring algorithm is based on two steps:

1) Estimate the PSF radius r from the blurred image b
2) Reconstruct x by image deconvolution with the estimated PSF radius


```bash
Full description of algorithm, pseudocode maybe?
```


## Installation instructions

The main deblurring algorithm is implemented in the `main_pipeline.m` function, which takes three input arguments:

1. (string) Folder where the input image files are located
2. (string) Folder where the output images must be stored
3. (int) Blur category number. Values between 0 and 19

Example of calling function: 
> `main_pipeline('path/to/input/files', 'path/to/output/files', 3)`

The `main_pipeline.m` function opens each image in the specified input folder, performs the deblurring using our algorithm, and then saves the image as a PNG file in the output folder with the same filename. 

The main algorithm needs the following additional Matlab functions to run
- `convb.m`  
- `r_update.m`
- `x_update.m`
- `FISTA_TVsmooth.m`
- `FISTA_TVsmooth_Chol_generator.m`
- egrssMatlab library containing 15 functions

Function descriptions:
- For performing image convolution with a given PSF radius

- Function for estimating PSF radius for current level of blur, given a blurred image b

- Function for performing image deblurring of a blurred image b, given the mean PSF radius estimate r

- Solves a smooth approximation to a total variation regularized least squares problem  

   min_x  1/2*||A(r)x-b||_2^2+ lambda*TV(x)

- Same as above, but with an added Cholesky term for including uncertainty in PSF radius estimate (denote the PSF radius estimate by r and the uncertainty term by mu)

   min_x  1/2*||L{A(r)x-b-mu}||_2^2+ lambda*TV(x)

- The egrssMatlab library is a linear algebra library for smarter matrix calculations, written by Associate Professor Martin S. Andersen at DTU Compute. 

## Usage instructions

For this competition we have been given both blurred and exact data with 20 different levels of blurring. This in principle allows us to estimate the PSF radius r for each blurring level, since we can use the first part of the algorithm (`r_update.m`) on a given pair of images (b,x) from each level. 

## Examples

