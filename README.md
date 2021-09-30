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

This GitHub repository has been created as the submission for group 1 (out of 2) of Team DTU for the Helsinki Deblur Challenge 2021.
The Matlab functions and scripts in the repository implements an image deblurring algorithm with Point-Spread-Function (PSF) radius estimation.
The algorithm is based on a paper written by our supervisors, where they develop a CT-reconstruction method, which uses uncertainty quantification to estimate view angles using only measured CT data [1].

## Description of algorithm

According to the work introduced in [1], we apply it for the image deblurring problem with an uncertain point spread function (PSF). Since the competition dataset as well as PSFs were taken by camera, PSFs are not accurate, which would influence the deblurring results significantly. Our focus is to test the idea from [1] on image deblurring problem and to see how to estimate PSF from degraded images. We assume that the blurring process is the out-of-focus blur, and PSF is a normalized disk with the radius r ∈ R. We apply the method proposed in [1] to estimate the radius r. The method combines the idea of the approximation error approach in Bayesian framework with the variational methods, and quantifies the uncertainty in r via a model-discrepancy term.

The foundation of the algorithm lies in the assumption that a blurred image b consists of two things:
- the exact image x convolved (out-of-focus blurred) with a PSF of radius r > 0
- a Gaussian noise term e, with mean 0 and some standard deviation sigma > 0

We denote the "convolution with a PSF of radius r" operator as A(r), and write the problem as

b = A(r)x + e 


First the algorithm is provided some pre-calculated initial PSF radius estimates, tailored to each step of the competition. This initial radius estimate can be further refined by using the method from [1] iteratively. Due to the high resolutions, in order to reduce the computational complexity, we apply the method in [1] to a small 2D patch of the degraded image. A small patch is extracted from the noisy blurred image. The patch is deblurred by solving a minimization problem with Total Variation (TV) regularization. In principle we can apply any variational methods to deblur the image. Here, we simply utilize the L2-TV method followed by some image enhancement techniques.

   min_x 1/2 * ||A(r)x - b||_2^2 + lambda * TV(x)

The small size of the patch makes it feasible to do uncertainty quantification of the radius estimate and the noise level. This is used to refine the radius estimate, which is then used to improve the deblurring of the patch. 

The algorithm alternates between refining the radius estimate and refining the patch for 10 iterations, after which the radius is assumed to have converged.
Finally, the full image is deblurred by solving a TV regularized minimization problem as above where we use the estimated PSF radius.


## Installation instructions
The main deconvolution algorithm is implemented in the Matlab function `main.m`. It needs the following additional Matlab routines to run correctly.

- `convb.m`  For performing image convolution with a given PSF radius

- `r_update.m`
Function for estimating PSF radius for current level of blur, given a blurred image b

- `x_update.m`
Function for performing image deconvolution of a blurred image b, given the mean PSF radius estimate r

- `FISTA_TVsmooth.m`
Solves a smooth approximation to a total variation regularized least squares problem  

   min_x  1/2*||A(r)x-b||_2^2+ lambda*TV(x)

- `FISTA_TVsmooth_Chol_generator.m`
Same as above, but with an added Cholesky term for including uncertainty in PSF radius estimate (denote the PSF radius estimate by r and the uncertainty term by mu)

   min_x  1/2*||L{A(r)x-b-mu}||_2^2+ lambda*TV(x)

- `egrssMatlab` library containing 15 functions.
The egrssMatlab library is a linear algebra library for smarter matrix calculations, written by Associate Professor Martin S. Andersen at DTU Compute. 


The script `estimate_psf_radius.m` performs some image analysis to estimate a crude initial PSF radius for each step, and the scripts `estimate_psf_radius_v2.m` uses the iterative deblurring/PSF radius estimation method adapted from [1] to do the same. These were run using the blurred PSF images for each steps provided by the competition organizers. These methods along with a grid-search for the optimal PSF radius and regularization parameter combination, have yielded some parameters for each step which we consider optimal by visual inspection. These scripts are not necessary for running the deblurring algorithm, but gives an insight into how we approached the problem of estimating the radius of a real-life point spread function. 

## Usage instructions
The main deblurring algorithm is implemented in the `main.m` function, which takes three input arguments:

1. (string) `inputfolder`: Folder where the input blurred image files are located (e.g 'path/to/input/files')
2. (string) `outputfolder`: Folder where the output deblurred images must be stored (e.g 'path/to/output/files')
3. (int) `step`: Blur category number. Values between 0 and 19

Example of calling function: 
> `main('path/to/input/files', 'path/to/output/files', 3)`

The `main.m` function opens each of the *.tif files in the specified input directory, performs the deblurring using our algorithm, and then saves the deblurred image as a *.PNG file in the output folder with the same filename (apart from the extension). 

## Examples

![examples](assets/steps_examples.png)


## References
[1] Nicolai Andre Brogaard Riis, Yiqiu Dong and Per Christian Hansen, Computed tomography with view angle estimation using uncertainty quantification, Inverse Problems, Vol. 37, pp. 065007, 2021.

