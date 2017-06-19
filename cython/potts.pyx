import numpy as np
cimport numpy as np

DTYPE = np.float
ctypedef np.float_t DTYPE_t

def phi(float u):
    return np.log(1 + u**2)

def neighbourhood_8(np.ndarray X, float colour, np.ndarray coords, int n=4, int order=1):
    assert X.dtype == DTYPE and coords.dtype == np.int
    cdef np.ndarray colours = np.zeros((8,), dtype=DTYPE)
    cdef int i = coords[0]
    cdef int j = coords[1]

    cdef int height = X.shape[0]
    cdef int width  = X.shape[1]

    cdef float s = 0
    cdef float w1 = 1
    cdef float w2 = 0.6

    if i-1 > -1:
        s += w1 * phi(X[i-1,j] - colour)/100
    if i+1 <  height:
        s += w1 * phi(X[i+1,j] - colour)/100
    if j+1 <  width:
        s += w1 * phi(X[i,j+1] - colour)/100
    if j-1 > -1:
        s += w1 * phi(X[i,j-1] - colour)/100

    if n == 8:
        if i-1 > -1 and j-1 > -1:
            s += w2 * phi(X[i-1,j-1] - colour)/100
        if i+1 <  height and j-1 > -1:
            s += w2 * phi(X[i+1,j-1] - colour)/100
        if i+1 <  height and j+1 <  width:
            s += w2 * phi(X[i+1,j+1] - colour)/100
        if i-1 > -1 and j+1 <  width:
            s += w2 * phi(X[i-1,j+1] - colour)/100

    return s

def mh_sample(np.ndarray Y, int n_iter, float noise_var=0.3, float beta=1, int neighbours=4):
    assert Y.dtype == DTYPE

    cdef n = Y.shape[0]
    cdef m = Y.shape[1]
    cdef np.ndarray X = np.random.uniform(0, 1, (n,m))

    cdef DTYPE_t x_candidate
    cdef DTYPE_t h
    cdef DTYPE_t h_candidate
    cdef DTYPE_t d
    cdef DTYPE_t d_candidate

    cdef DTYPE_t p
    cdef DTYPE_t U

    # Generate a chain of length n_iter
    for k in range(n_iter):
        # Iterate over each pixel coordinate (i,j)
        for i in range(n):
            for j in range(m):
                x_candidate = max(min(X[i,j] + np.random.normal(0, 1), 1), 0) #np.random.choice(colors)

                h             = - 1.0/(2*noise_var) * (Y[i,j] - X[i,j])**2
                h_candidate   = - 1.0/(2*noise_var) * (Y[i,j] - x_candidate)**2

                d             = beta * neighbourhood_8(X, X[i,j], np.array([i, j]), n=neighbours) + h
                d_candidate   = beta * neighbourhood_8(X, x_candidate, np.array([i, j]), n=neighbours) + h_candidate

                # Acceptance probability
                p = np.exp(min(0, d_candidate - d))
                U = np.random.uniform()
                if U < p:
                    X[i,j] = x_candidate

    return X
