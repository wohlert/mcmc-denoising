import numpy as np
from scipy.misc import imsave

#shape = (320, 480)
shape = (-1, 400)

A = np.loadtxt("denoised.txt")
imsave("denoised.jpg", A.reshape(shape))

import glob

files = glob.glob("history/*.txt")

for i, file in enumerate(files):
    A = np.loadtxt(file)
    imsave("history/history_{}.jpg".format(i), A.reshape(shape))
