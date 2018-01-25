Hyperspectral ReadMe

The repository contains the Matlab implementations of spatial compositional model (SCM), Gaussian mixture model (GMM), which are used in hyperspectral unmixing. If you find the code helpful, please cite the following papers: 

Zhou, Y., Rangarajan, A. & Gader, P. D. A spatial compositional model for linear unmixing and endmember uncertainty estimation. IEEE Transactions on Image Processing, vol. 25, no. 12, pp. 5987-6002, 2016  
(http://ieeexplore.ieee.org/document/7592431/)

Zhou, Y., Rangarajan, A. & Gader, P. D. A Gaussian mixture model representation of endmember variability in hyperspectral unmixing. IEEE Transactions on Image Processing, To Appear, 2018  
(http://ieeexplore.ieee.org/document/8264812/)

-------------------------------- Summary ----------------------------------

The folder "common" contains common functions that will be used by all the algorithms.  
The folder "SCM" contains functions that implement the spatial compositional model.  
The folder "GMM" contains functions that implement the Gaussian mixture model for spectral unmixing.  

--------------------------------- GMM -------------------------------------

The "GMM" folder contains the files for the algorithm. There are two demo files. 

"test_gmm.m" - demo file that runs the unsupervised GMM algorithm, which is implemented in "gmm_hu.m", which first segments the image, uses the interior pixels of the segmented regions to build distribution parameters, and finally updates the abundances. It also calls "gmm_hu_endmember.m" which estimates pixelwise endmembers.

"test_gmm_ex.m" - demo file that runs the supervised GMM algorithm (implemented in "gmm_hu_ex.m"), which takes a library of spectra as input and outputs the abundances. It also calls "gmm_hu_endmember.m" to estimate pixelwise endmembers.

--------------------------------- SCM -------------------------------------

The "SCM" folder contains the main files for the algorithm. The common folder contains auxiliary files used for displaying and processing. To see the demo, run "test_scm.m" and you can change the "dataset" in it.

The SCM function is used as follows: 

function [A,R,mu,sigma,var_dirs,var_amts] = scm(I,M,options)  
Input:  
  I: row*col*B image data,  
  M: number of endmembers,  
  options: structure for additional parameters,  
Output:  
  A: abundances (N by M),  
  R: endmembers (M by B),  
  mu: noise standard deviations,  
  sigma: covariance matrices of endmember uncertainty,  
  var_dirs: uncertainty directions,  
  var_amts: uncertainty amounts.  

Note: options is a structure containing eta,beta1,beta2,rho1,show_figure,init_mode, etc.

-------------------------------- Contact ----------------------------------

If you have any questions, please contact:

Yuan Zhou  
Department of Computer Information Science and Engineering  
University of Florida  

zhouyuanzxcv@gmail.com

