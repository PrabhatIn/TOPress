# TOPress
## TOPress: a MATLAB implementation for topology optimization of structures subjected to design‑dependent pressure loads
`TOPress.m` provides a MATLAB implementation for topology optimization of structures subjected to design‑dependent pressure loads. Typically, a design-dependent load changes its direction, locationa and/or magnitude as topology optimization advances, and, thus, poses several unique challanges.   
## About
Author: Prabhat Kumar, Department of Mechanical and Aerospace Engineering, Indian Institute of Technology Hyderabad, India. Please send your comments and suggestions to  pkumar@mae.iith.ac.in or prabhatkumar.rns@gmail.com
## How to use
1. Please see the paper: P. Kumar (2022) TOPress: a MATLAB implementation for topology optimization of structures subjected to design‑dependent pressure loads, Structural and Multidisciplinary Optimization, 66(4), 2023
2. MMA setting:
TOPress code uses the MMA written in 1999 and updated in 2002 version. The mmasub function has the followinig form
[xmma,ymma,zmma,lam,xsi,eta,mu,zet,s,low,upp] = mmasub(m,n,iter,xval,xmin,xmax,xold1,xold2,f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0,a,c,d);
in TOPress code calls it  on line 91 as
[xmma, ~,~,~,~,~,~,~,~,low,upp] = mmasub(mMMA,nMMA,loop,xval,xminvec,xmaxvec,xold1,xold2, ...
        obj*normf,objsens(act),objsens(act)*0,Vol,dVol(act),dVol(act)*0,low,upp,a0,aMMA,cMMA,dMMA); % calling MMA
## Citation
For citing the paper, one may use the following bibtex format:
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
