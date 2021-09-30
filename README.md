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

The main deblurring algorithm is implemented in the `main_pipeline.m` function, which takes three input arguments:

1. (string) Folder where the input image files are located
2. (string) Folder where the output images must be stored
3. (int) Blur category number. Values between 0 and 19

Example of calling function: 
> `main_pipeline('path/to/input/files', 'path/to/output/files', 3)`

The `main_pipeline.m` function opens each image in the specified input folder, performs the deblurring using our algorithm, and then saves the image as a PNG file in the output folder with the same filename. 

The main algorithm needs the following additional Matlab functions to run
- `convb.m`  For performing image convolution with a given PSF radius

- `r_update.m`
Function for estimating PSF radius for current level of blur, given a blurred image b

- `x_update.m`
Function for performing image deblurring of a blurred image b, given the mean PSF radius estimate r

- `FISTA_TVsmooth.m`
Solves a smooth approximation to a total variation regularized least squares problem  

   min_x  1/2*||A(r)x-b||_2^2+ lambda*TV(x)

- `FISTA_TVsmooth_Chol_generator.m`
Same as above, but with an added Cholesky term for including uncertainty in PSF radius estimate (denote the PSF radius estimate by r and the uncertainty term by mu)

   min_x  1/2*||L{A(r)x-b-mu}||_2^2+ lambda*TV(x)

- `egrssMatlab` library containing 15 functions
The egrssMatlab library is a linear algebra library for smarter matrix calculations, written by Associate Professor Martin S. Andersen at DTU Compute. 

Function descriptions:





## Usage instructions

For this competition we have been given both blurred and exact data with 20 different levels of blurring. This in principle allows us to estimate the PSF radius r for each blurring level, since we can use the first part of the algorithm (`r_update.m`) on a given pair of images (b,x) from each level. 

## Examples

![examples](assets/steps_examples.png)


## References
[1] Nicolai Andre Brogaard Riis, Yiqiu Dong and Per Christian Hansen, Computed tomography with view angle estimation using uncertainty quantification, Inverse Problems, Vol. 37, pp. 065007, 2021.

