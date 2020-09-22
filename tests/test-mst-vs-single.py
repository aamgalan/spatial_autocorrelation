# testing to see if 2D implementation of minimum spanning tree clustering 
# avoids the ~n^2/2 memory issue of the scipy implementation 
# scipy implementation is more general than just the 
# 2D space our implementation works in 

import sys, os
import numpy as np

# adding the directory where S_A is located 
path_S_A = os.path.realpath(__file__).split('/tests')[0]
sys.path.append(path_S_A)

from S_A import *

for size in [10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000]:
    
    x = np.random.randn(size, 2)
    z = np.random.randn(size)
    
    print('handwritten MST clustering at n =', size)
    merge_order = get_mst_merge_order(x)
    skienaA = SkienaA(z, merge_order)
    print(skienaA)

    print('scipy single-linkage clustering at n =', size)
    try:
        merge_order = get_merge_order(x, 'single')
        skienaA = SkienaA(z, merge_order)
        print(skienaA)

    except MemoryError:
        print('memory issue')
    
