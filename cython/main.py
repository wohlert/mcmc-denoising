import potts
from skimage import data, util

X = data.horse()
Y = util.random_noise(X, "gaussian", var=0.3)

U = potts.mh_sample(Y, 10)

from IPython import embed

embed()
