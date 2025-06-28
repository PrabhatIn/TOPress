# TOPress
## TOPress: a MATLAB implementation for topology optimization of structures subjected to design‑dependent pressure loads
`TOPress.m` provides a MATLAB implementation for topology optimization of structures subjected to design‑dependent pressure loads. Typically, a design-dependent load changes its direction, locationa and/or magnitude as topology optimization advances, and, thus, poses several unique challenges.   
## About
Author: Prabhat Kumar, Department of Mechanical and Aerospace Engineering, Indian Institute of Technology Hyderabad, India. Please send your comments and suggestions to  pkumar@mae.iith.ac.in or prabhatkumar.rns@gmail.com
## How to use
1. Please see the paper: [P. Kumar (2023) TOPress: a MATLAB implementation for topology optimization of structures subjected to design‑dependent pressure loads, Structural and Multidisciplinary Optimization, 66(4), 2023](https://link.springer.com/article/10.1007/s00158-023-03533-9)
2. MMA setting:
   
  
(i) TOPress uses the MMA written in 1999 and updated in the 2002 version. The mmasub function has the following form
## [xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2,f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c,d);
 TOPress code calls it  on line 91 as
## [xmma,~,~,~,~,~,~,~,~,low,upp]=mmasub(mMMA,nMMA,loop,xval,xminvec,xmaxvec,xold1,xold2,obj*normf,objsens(act),objsens(act)*0,Vol,dVol(act)',dVol(act)'*0,low,upp,a0,aMMA,cMMA,dMMA);
(ii) 
With the 2006 version of MMA, one can modify MMA call (line 91) to (may have to reduce 'change' on line 60):
## [xmma,~,~,~,~,~,~,~,~,low,upp]=mmasub(mMMA,nMMA,loop,xval,xminvec,xmaxvec,xold1,xold2,obj*normf,objsens(act),Vol,dVol(act)',low,upp,a0,aMMA,cMMA,dMMA);
## Citation
For citing the paper, please use the following bibtex format:
```
@article{kumar2023TOPress,
  title={{TOPress}: a {MATLAB} implementation for topology optimization of structures subjected to design‑dependent pressure loads},
  author={Kumar, Prabhat},
  journal={Structural and Multidisciplinary Optimization},
  volume={66},
  number={4},
  year={2023},
  publisher={Springer}
}
```
