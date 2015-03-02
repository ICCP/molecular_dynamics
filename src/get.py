import math
import numpy as np


# Returns the standard deviation in the last 100 components of given vector x
def standard_deviation(x):
    y = np.zeros(100)
    for i in xrange(100):
        y[i] = x[len(x)-i-1]
    return np.std(y)


# Returns the mean in the last 100 components of given vector x
def mean(x):
    m=0
    for i in xrange(100):
        m += x[len(x)-i-1]/100.0
    return m

# Returns the variance in the last 100 components of given vector x
def var(x):
    y = np.zeros(100)
    for i in xrange(100):
        y[i] = x[len(x)-i-1]
    return np.var(y)
