# ta2s2_codes
Scripts for the Transitional Annealed Adaptive Slice Sampling method

### Overview
This repository contains an extension to the Asymptotically Independent Markov Sampling (AIMS) algorithm for stochastic optimisation purposes and has been currently developed for the estimation of the hyper-parameters of a Gaussian process emulator. This work is part of my PhD research project at the [University of Liverpool](https://www.liv.ac.uk/risk-and-uncertainty/). Currently the only version available is on `Matlab` but it is expected to be written in `R` and `Python` in the near future. 

### Contents
* A *startup* file to load all directories needed for the sampler.  
* The `matlab-files` directory contains all files needed for the examples to run. Being  
```Matlab
    [ ... ] = parallel_aims_opt( ... ) 
```
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;the main code for the sampler.  
* The `examples` directory contains several examples used in the paper submitted with the results obtained from such extension.  
* To run the examples first run the **startup** script and run one of the executables in the example directory.  
```Matlab
   startup; execute;        % To run the Franke's emulator example.
```

### Citing
We would be grateful if any results based on this transitional extension of the Slice Sampling algorithm are acknowledged by citing our paper. It is going to be available on [arxiv](http://arxiv.org)

```TeX
@Article{Garbuno2015,
  Title                    = {{Transitional annealed adaptive slice sampling for Gaussian process hyper-parameter 
estimation}},
  Author                   = {Garbuno-Inigo, Alfredo, and DiazDelaO, F.~A., and Zuev, K.~M.},
  Journal                  = {Arxiv},
  Year                     = {2015},

  Archiveprefix            = {arXiv},
  Arxivid                  = {1509.00349},
  Url                      = {http://arxiv.org/abs/1509.00349}
}
```
