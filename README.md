Hyperspectral ReadMe

The repository contains the Matlab implementations of several algorithms in hyperspectral image analysis, including:  
	1. Spatial compositional model (SCM) for unmixing with a fixed endmember set (TIP16)  
	2. Gaussian mixture model (GMM) for unmixing with endmember variability (TIP18)  
	3. A registration and fusion algorithm for combining a hyperspectral image and a multispectral image (TGRS19)  
If you find some of the code helpful, please cite the corresponding papers:  


-------------------------------- Summary ----------------------------------

The folder "common" contains common functions that will be used by all the algorithms.  

The folder "SCM" contains functions that implement the spatial compositional model.  

Zhou, Y., Rangarajan, A. & Gader, P. D. A spatial compositional model for linear unmixing and endmember uncertainty estimation. IEEE Transactions on Image Processing, vol. 25, no. 12, pp. 5987-6002, 2016  
(http://ieeexplore.ieee.org/document/7592431/)  

The folder "GMM" contains functions that implement the Gaussian mixture model for spectral unmixing.  

Zhou, Y., Rangarajan, A. & Gader, P. D. A Gaussian mixture model representation of endmember variability in hyperspectral unmixing. IEEE Transactions on Image Processing, vol. 27, no. 5, pp. 2242-2256, 2018  
(http://ieeexplore.ieee.org/document/8264812/)  

The folders "REG" and "Fusion" contains functions that implement the registration algorithm and the fusion algorithm for hyperspectral and multispectral images.  

Zhou, Y., Rangarajan, A. & Gader, P. D. An integrated approach to registration and fusion of hyperspectral and multispectral images, IEEE Transactions on Geoscience and Remote Sensing, To Appear, 2019  
(https://ieeexplore.ieee.org/document/8897135)  


--------------------------------- SCM -------------------------------------

The "SCM" folder contains the main files for the unmixing algorithm with uncertainty estimation. The common folder contains auxiliary files used for displaying and processing. To see the demo, run "test_scm.m" and you can change the "dataset" in it.

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

--------------------------------- GMM -------------------------------------

The "GMM" folder contains the files for the unmxing algorithm with endmember variability modeled by GMM. There are two demo files.  

"test_gmm.m" - demo file that runs the unsupervised GMM algorithm, which is implemented in "gmm_hu.m". The algorithm first segments the image, uses the interior pixels of the segmented regions to build distribution parameters, and finally updates the abundances. It also calls "gmm_hu_endmember.m" which estimates pixelwise endmembers.  

"test_gmm_ex.m" - demo file that runs the supervised GMM algorithm (implemented in "gmm_hu_ex.m"). It takes a library of spectra as input and outputs the abundances. It also calls "gmm_hu_endmember.m" to estimate pixelwise endmembers.  

--------------------------------- REG -------------------------------------

The "REG" folder contains the files for registration of hyperspectral and multispectral images.  

"test_reg.m" - demo file that runs the registration algorithm, which can handle images with different spatial and spectral resolutions.  

"reg_hyper_rgb.m" - main file that implements the registration algorithm.  

--------------------------------- Fusion -------------------------------------

The "Fusion" folder contains the files for fusion of hyperspectral and multispectral images.  

"run_fusion_algo.m" - demo file that runs the fusion algorithm.  

"im_fusion.m" - main file that implements the fusion algorithm.  


-------------------------------- Contact ----------------------------------

If you have any questions, please contact:

Yuan Zhou  
Department of Computer and Information Science and Engineering  
University of Florida  

zhouyuanzxcv@gmail.com

