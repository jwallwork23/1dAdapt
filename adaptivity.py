from firedrake import *

import numpy as np


def intervalAdapt(mesh, f, tolCoarsen, tolRefine):
    """
    :param mesh: interval mesh.
    :param f: function defined on interval mesh.
    :param tolCoarsen: tolerance below which mesh is to be coarsened.
    :param tolRefine: tolerance above which mesh is to be refined.
    :return: adapted mesh.
    """
    x1 = mesh.coordinates.dat.data      # Get mesh coordinates
    d = x1[-1] - x1[0]                   # Get length of interval
    # x1.append(x1[0])                    # Append coordinate list for periodic boundary
    x2 = []                             # List for new coordinate set
    fData = f.dat.data
    # fData.append(fData[0])              # Append function value list for periodic boundary

    for i in range(len(x1)-1):
        fAverage = 0.5 * (fData[i+1] + fData[i])    # Get average function value over cell
        if (i == 0) | (fData[i] > tolCoarsen):
            x2.append(x1[i])                        # Add in any nodes to be kept (implicitly deleting otherwise)
        if fAverage > tolRefine:
            x2.append(0.5 * (x1[i+1] + x1[i+1]))    # Add in any extra nodes for refinement
    x2.append(x1[-1])
        
    mesh2 = IntervalMesh(len(x2)-1, d)              # Establish a new mesh object
    mesh2.coordinates.dat.data[:] = x2              # Set the new coordinates

    return mesh2

# TODO: above does not work for PeriodicIntervalMesh, as the .dat.data structure is different

# def interpolateVariable(mesh, f):
#     # TODO

if __name__ == '__main__':
    mesh = IntervalMesh(8, 1)
    print('Mesh coords pre-adapt: ', mesh.coordinates.dat.data)
    X = SpatialCoordinate(mesh)
    f = Function(FunctionSpace(mesh, 'CG', 1)).interpolate(X**2)

    mesh2 = intervalAdapt(mesh, f, 0.1, 0.5)
    print('Mesh coords post-adapt: ', mesh2.coordinates.dat.data)

