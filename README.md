# SLPMF
Fast Non-negative Matrix Factorization under KL-Divergence: A Case for Link Prediction

This is the implementation of the SLPMF framework presented in [1].

## Installation & Usage
* For GSCD:
```python
function [U,V,TraceU,TraceV]=GSCD(incompleteNetwork, unknownsA, latentDim, coordGradMaxIter, gradMaxIter, gradeps,cyclic)
%Inputs:
%   incompleteNetwork - incomplete matrix
%   unknownsA - list of unknown entries
%   latentDim - latent dimension
%   coordGradMaxIter - coordinate descent iterations
%   gradMaxIter - gradient descent iterations
%   gradeps - gradient descent epsilon
%   cyclic - 1=cyclic optimization, 1=greedy optimization
%
%Outputs:
%   U,V - latent factors
%   TraceU - Objective value Traces (U Iterations)
%   TraceV - Objective value Traces (V Iterations)
```

* For PGSCD:
```python
function [U,V,TraceU,TraceV]=PGSCD(incompleteNetwork, unknownsA, latentDim, coordGradMaxIter, gradMaxIter, gradeps, numberOfThreads)
%Inputs:
%   incompleteNetwork - incomplete matrix
%   unknownsA - list of unknown entries
%   latentDim - latent dimension
%   coordGradMaxIter - coordinate descent iterations
%   gradMaxIter - gradient descent iterations
%   gradeps - gradient descent epsilon
%   numberOfThreads - number of parallel threads
%
%Outputs:
%   U,V - latent factors
%   TraceU - Objective value Traces (U Iterations)
%   TraceV - Objective value Traces (V Iterations)

```

## References
* [1] P. Siyari, and H. R. Rabiee, “Fast Non-negative Matrix Factorization under KL-Divergence: A Case for Link Prediction”, Sharif University of Technology, Technical Report, 2013.