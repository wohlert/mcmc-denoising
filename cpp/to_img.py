import numpy as np
from scipy.misc import imsave

A = np.loadtxt("denoised.txt")
imsave("denoised.jpg", A.reshape(-1, 400))

import glob

files = glob.glob("history/*.txt")

for i, file in enumerate(files):
    A = np.loadtxt(file)
    imsave("history/history_{}.jpg".format(i), A.reshape(-1, 400))
