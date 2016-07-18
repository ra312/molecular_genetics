from numpy import array
from scipy.cluster.vq import vq, kmeans, whiten
features  = array([[ 1.9,2.3],
                    [ 1.5,2.5],
                    [ 0.8,0.6],
                    [ 0.4,1.8],
                    [ 0.1,0.1],
                    [ 0.2,1.8],
                    [ 2.0,0.5],
                    [ 0.3,1.5],
                    [ 1.0,1.0]])

print "Features"
print features
whitened = whiten(features)
print "Features after whitening"
print features
book = array((whitened[0],whitened[2]))
print "book"
print book
kmeans(whitened,book)
(array([[ 2.3110306 ,  2.86287398],    # random
       [ 0.93218041,  1.24398691]]), 0.85684700941625547)

from numpy import random

random.seed((1000,2000))
codes = 3
print kmeans(whitened,codes)

(array([[ 2.3110306 ,  2.86287398],    # random
       [ 1.32544402,  0.65607529],
       [ 0.40782893,  2.02786907]]), 0.5196582527686241)
