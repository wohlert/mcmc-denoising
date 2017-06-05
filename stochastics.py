import itertools

def lcg(x_i, a, c, M):
    while True:
        x_i = (a * x_i + c) % M
        yield x_i/M


lcg_generator = lcg(0, 135412, 5121, 1712289)
take = lambda f, n: list(itertools.islice(f, n))
colors = ["#1A2F4B", "#28475C", "#2F8886", "#84C69B"]


def up_down(x):
    n = len(x)

    def count_occurences(x):
        n = len(x)
        inf = 99999999
        y = x.copy() + [-inf]

        count = 1
        for i in range(n):
            if y[i+1] < y[i]:
                yield count
                count = 0
            count += 1


    A = np.array([
              [ 4529.4,  9044.9, 13568,  18091,  22615,  27892],
              [ 9044.9, 18097,   27139,  36187,  45234,  55789],
              [13568,   27139,   40721,  54281,  67852,  83685],
              [18091,   36187,   54281,  72414,  90470, 111580],
              [22615,   45234,   67852,  90470, 113262, 139476],
              [27892,   55789,   83685, 111580, 139476, 172860]])

    B = np.array([1/6, 5/24, 11/120, 19/720, 29/5040, 1/840])

    occurences = list(count_occurences(x))

    R =  [sum(filter(lambda x: x == i, occurences)) for i in range(1, 6)]
    R += [sum(filter(lambda x: x >= 6, occurences))]

    Z = 1/(n - 6) * (R - n*B).T @ A @ (R - n*B)

    return R, Z


def head(l):
    return l[0]


def tail(l):
    return l[1:]