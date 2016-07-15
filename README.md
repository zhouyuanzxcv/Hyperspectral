Hyperspectral ReadMe

The repository contains the Matlab implementation of spatial compositional model(SCM), which is used in hyperspectral unmixing. If you find the code helpful, please cite the following paper: 

Zhou, Y.; Rangarajan, A. & Gader, P. D. A spatial compositional model (SCM) for linear unmixing and endmember uncertainty estimation. arXiv preprint arXiv:1509.09243, 2015

The SCM folder contains the main files for the algorithm. The common folder contains auxiliary files used for displaying and processing. To see the demo, run "test_scm.m" and you can change the "dataset" in it.

The SCM function is used as follows: 

function [A,R,mu,sigma,var_dirs,var_amts] = scm(I,M,options)
Input:
  I: row*col*B image data
  M: number of endmembers
  options: structure for additional parameters
Output:
  A: abundances (N by M)
  R: endmembers (M by B)
  mu: noise std
  sigma: covariance matrices of endmember uncertainty
  var_dirs: uncertainty directions
  var_amts: uncertainty amounts

Note: options is a structure containing eta,beta1,beta2,rho1,...
  show_figure,init_mode, etc.

If you have any questions, please contact:

Yuan Zhou 
Department of Computer Information Science and Engineering
University of Florida

zhouyuanzxcv@gmail.com

