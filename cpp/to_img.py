import numpy as np
from scipy.misc import imsave

A = np.loadtxt("denoised.txt")
imsave("denoised.jpg", A.reshape(-1, 400))
