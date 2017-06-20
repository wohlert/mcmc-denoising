# Bayesian image denoising

This project uses Markov-chain Monte Carlo (MCMC) methods to sample from a prior distribution of images in order to perform denoising on images. The algorithms perform well on images with *speckle* noise as well as *salt & pepper* noise.

This project is built up of two C++ classes, `Ising` and `Potts` for binary and grayscale images respectively. Additionally, we provide a Cython bridge that enables use directly from Python.

## Examples

**Original image**:

![](data/lars-noisy-gray.jpg)

**Denoised image using Metropolis-Hastings based on the Potts model**.

![](data/lars-denoised.jpg)


## How to use

Compile the Cython library by running.

```bash
python setup.py build_ext --inplace
```

This will generate the library which you can import from Python.

```python
import numpy as np
import matplotlib.pyplot as plt

from denoising import IsingMH, PottsMH

lars = plt.imread("data/lars-noisy-gray.jpg")[:, :, 0]/255

potts = denoising.PottsMH(lars, beta=4, sigma=np.std(lars), bins=10)
denoised = potts.solve(iterations=10)
```
