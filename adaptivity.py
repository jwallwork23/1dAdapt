from firedrake import *

import numpy as np
import matplotlib.pyplot as plt


def intervalAdapt(mesh, f, tolCoarsen, tolRefine):
    """
    :param mesh: interval mesh.
    :param f: function defined on interval mesh.
    :param tolCoarsen: tolerance below which mesh is to be coarsened.
    :param tolRefine: tolerance above which mesh is to be refined.
    :return: adapted mesh.
    """
    # TODO: only currently works for P1 fields

    x1 = mesh.coordinates.dat.data      # Get mesh coordinates
    d = x1[-1] - x1[0]                  # Get length of interval
    x2 = []                             # List for new coordinate set
    fData = f.dat.data

    for i in range(len(x1)-1):
        if (i == 0):
            x2.append(x1[i])
        elif (fData[i-1] > tolCoarsen) | (fData[i] > tolCoarsen) | (fData[i+1] > tolCoarsen):
            x2.append(x1[i])                        # Add in any nodes to be kept (implicitly deleting otherwise)
        if (fData[i] > tolRefine) & (fData[i+1] > tolRefine):
            x2.append(0.5 * (x1[i] + x1[i+1]))    # Add in any extra nodes for refinement
    x2.append(x1[-1])
        
    mesh2 = IntervalMesh(len(x2)-1, d)              # Establish a new mesh object
    mesh2.coordinates.dat.data[:] = x2              # Set the new coordinates

    return mesh2

# TODO: above does not work for PeriodicIntervalMesh, as the .dat.data structure is different

def interpolateVariable(mesh, f):
    """
    :param mesh: new mesh to interpolate onto.
    :param f: function to interpolate.
    :return: function f interpolated onto new mesh.
    """
    V = f.function_space()
    family = V.ufl_element().family()
    try:
        assert family in ('Lagrange', 'Discontinuous Lagrage')
    except:
        raise NotImplementedError('Function space type %s not yet implemented' % family)
    g = Function(FunctionSpace(mesh, family, V.ufl_element().degree()))
    for i in range(len(mesh.coordinates.dat.data)):
        g.dat.data[i] = f.at(mesh.coordinates.dat.data[i])
    return g


if __name__ == '__main__':
    mesh = IntervalMesh(8, np.pi)
    X, = SpatialCoordinate(mesh)
    f = Function(FunctionSpace(mesh, 'CG', 1)).interpolate(sin(X))
    pre_adapt_mesh = mesh.coordinates.dat.data
    plt.plot(pre_adapt_mesh, f.dat.data)
    plt.plot(pre_adapt_mesh, np.zeros(len(pre_adapt_mesh)), 'o')

    mesh2 = intervalAdapt(mesh, f, 0.1, 0.5)
    f2 = interpolateVariable(mesh2, f)
    post_adapt_mesh = mesh2.coordinates.dat.data
    plt.plot(post_adapt_mesh, f2.dat.data)
    plt.plot(post_adapt_mesh, np.zeros(len(post_adapt_mesh)) - 0.02, 'x')
    plt.show()
