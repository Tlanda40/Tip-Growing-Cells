import numpy as np
import gdist

'''
This is a python program that calculates geodesic distances from
the pole to any point on the surface. It uses the gdist python 
package implementation of the MMP algorithm.
'''

distances = gdist.compute_gdist(
    numpy.ndarray[numpy.float64_t, ndim=2] vertices, 
    numpy.ndarray[numpy.int32_t, ndim=2] triangles, 
    numpy.ndarray[numpy.int32_t, ndim=1] source_indices = None, 
    numpy.ndarray[numpy.int32_t, ndim=1] target_indices = None, 
    double max_distance = GEODESIC_INF, #change GEODESIC_INF to a real value to restrict recalculated geodesics
    bool is_one_indexed = True #matlab starts indexing at 1
)