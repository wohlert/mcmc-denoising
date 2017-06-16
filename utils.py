def mode(ndarray, axis=0):
    # Check inputs
    ndarray = np.asarray(ndarray)
    ndim = ndarray.ndim
    if ndarray.size == 1:
        return (ndarray[0], 1)
    elif ndarray.size == 0:
        raise Exception('Cannot compute mode on empty array')
    try:
        axis = range(ndarray.ndim)[axis]
    except:
        raise Exception('Axis "{}" incompatible with the {}-dimension array'.format(axis, ndim))

    # If array is 1-D and numpy version is > 1.9 numpy.unique will suffice
    if all([ndim == 1,
            int(np.__version__.split('.')[0]) >= 1,
            int(np.__version__.split('.')[1]) >= 9]):
        modals, counts = np.unique(ndarray, return_counts=True)
        index = np.argmax(counts)
        return modals[index], counts[index]

    # Sort array
    sort = np.sort(ndarray, axis=axis)
    # Create array to transpose along the axis and get padding shape
    transpose = np.roll(np.arange(ndim)[::-1], axis)
    shape = list(sort.shape)
    shape[axis] = 1
    # Create a boolean array along strides of unique values
    strides = np.concatenate([np.zeros(shape=shape, dtype='bool'),
                                 np.diff(sort, axis=axis) == 0,
                                 np.zeros(shape=shape, dtype='bool')],
                                axis=axis).transpose(transpose).ravel()
    # Count the stride lengths
    counts = np.cumsum(strides)
    counts[~strides] = np.concatenate([[0], np.diff(counts[~strides])])
    counts[strides] = 0
    # Get shape of padded counts and slice to return to the original shape
    shape = np.array(sort.shape)
    shape[axis] += 1
    shape = shape[transpose]
    slices = [slice(None)] * ndim
    slices[axis] = slice(1, None)
    # Reshape and compute final counts
    counts = counts.reshape(shape).transpose(transpose)[slices] + 1

    # Find maximum counts and return modals/counts
    slices = [slice(None, i) for i in sort.shape]
    del slices[axis]
    index = np.ogrid[slices]
    index.insert(axis, np.argmax(counts, axis=axis))
    return sort[index]

def neighbourhood_second_order(X, colour, coords, phi, n=4, order=1):
    colours = []
    i, j = coords
    H, W = X.shape

    if i-1 > -1:
        colours.append(X[i-1,j])
    if i+1 <  H:
        colours.append(X[i+1,j])
    if j+1 <  W:
        colours.append(X[i,j+1])
    if j-1 > -1:
        colours.append(X[i,j-1])

    if n == 8:
        if i-1 > -1 and j-1 > -1:
            colours.append(X[i-1, j-1])
        if i+1 <  H and j-1 > -1:
            colours.append(X[i+1, j-1])
        if i+1 <  H and j+1 <  W:
            colours.append(X[i+1, j+1])
        if i-1 > -1 and j+1 <  W:
            colours.append(X[i-1, j+1])

    first_order = [phi(next_colour - colour)/10 for next_colour in colours]

    return 1 * sum(first_order)
