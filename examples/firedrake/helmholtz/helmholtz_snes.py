# Simple Helmholtz equation but cast as a nonlinear problem
# =========================
#
# Let's start by considering the modified Helmholtz equation on a unit square,
# :math:`\Omega`, with boundary :math:`\Gamma`:
#
# .. math::
#
#    -\nabla^2 u + u &= f
#
#    \nabla u \cdot \vec{n} &= 0 \quad \textrm{on}\ \Gamma
#
# for some known function :math:`f`. The solution to this equation will
# be some function :math:`u\in V`, for some suitable function space
# :math:`V`, that satisfies these equations. Note that this is the
# Helmholtz equation that appears in meteorology, rather than the
# indefinite Helmholtz equation :math:`\nabla^2 u + u = f` that arises
# in wave problems.
#
# We transform the equation into weak form by multiplying by an arbitrary
# test function in :math:`V`, integrating over the domain and then
# integrating by parts. The variational problem so derived reads: find
# :math:`u \in V` such that:
#
# .. math::
#
#    \require{cancel}
#    \int_\Omega \nabla u\cdot\nabla v  + uv\ \mathrm{d}x = \int_\Omega
#    vf\ \mathrm{d}x + \cancel{\int_\Gamma v \nabla u \cdot \vec{n} \mathrm{d}s}
#
# Note that the boundary condition has been enforced weakly by removing
# the surface term resulting from the integration by parts.
#
# We can choose the function :math:`f`, so we take:
#
# .. math::
#
#    f = (1.0 + 8.0\pi^2)\cos(2\pi x)\cos(2\pi y)
#
# which conveniently yields the analytic solution:
#
# .. math::
#
#    u = \cos(2\pi x)\cos(2\pi y)
#
# However we wish to employ this as an example for the finite element
# method, so lets go ahead and produce a numerical solution.
#
# First, we always need a mesh. Let's have a :math:`10\times10` element unit square::

from firedrake import *

N = 10  
mesh = UnitSquareMesh(N, N, quadrilateral=True)

# We need to decide on the function space in which we'd like to solve the
# problem. Let's use piecewise linear functions continuous between
# elements::

V = FunctionSpace(mesh, "CG", 2)

# We'll also need the test and trial functions corresponding to this
# function space::

u = Function(V)
v = TestFunction(V)

# We declare a function over our function space and give it the
# value of our right hand side function::

f = Function(V)
x, y = SpatialCoordinate(mesh)
f.interpolate((1+8*pi*pi)*cos(x*pi*2)*cos(y*pi*2))

# We can now define the weak form of the residual function F(u)::

F = (inner(grad(u), grad(v)) + inner(u, v)) * dx - inner(f, v) * dx

# set the petsc solver parameters to use the nonlinear solver SNES with a direct LU factorisation preconditioner. We also set some additional parameters to monitor the convergence of the nonlinear and linear solvers and to set the relative tolerance for convergence. For more details on how to specify solver parameters, see the section of the manual on :doc:`solving PDEs <../solving-interface>`.

parameters = {'ksp_type': 'preonly', 
              'pc_type': 'lu',
              'pc_factor_mat_solver_type': 'mumps',
              'snes_monitor': None,
              'snes_converged_reason': None,
              'snes_rtol': 1e-6,
              'snes_view': None,
              'ksp_monitor': None,
              'ksp_converged_reason': None,
              'ksp_rtol': 1e-6}

solve(F == 0, u, solver_parameters=parameters)


# For more details on how to specify solver parameters, see the section
# of the manual on :doc:`solving PDEs <../solving-interface>`.
#
# Next, we might want to look at the result, so we output our solution
# to a file::

VTKFile("helmholtz.pvd").write(u)

# This file can be visualised using `paraview <http://www.paraview.org/>`__.
#
# We could use the built-in plotting functions of firedrake by calling
# :func:`tripcolor <firedrake.pyplot.tripcolor>` to make a pseudo-color plot.
# Before that, matplotlib.pyplot should be installed and imported::

try:
  import matplotlib.pyplot as plt
except:
  warning("Matplotlib not imported")

try:
  from firedrake.pyplot import tripcolor, tricontour
  fig, axes = plt.subplots()
  colors = tripcolor(u, axes=axes)
  fig.colorbar(colors)
except Exception as e:
  warning("Cannot plot figure. Error msg: '%s'" % e)

# The plotting functions in Firedrake mimic those of matplotlib; to produce a
# contour plot instead of a pseudocolor plot, we can call
# :func:`tricontour <firedrake.pyplot.tricontour>` instead::

try:
  fig, axes = plt.subplots()
  contours = tricontour(u, axes=axes)
  fig.colorbar(contours)
except Exception as e:
  warning("Cannot plot figure. Error msg: '%s'" % e)

# Don't forget to show the image::

try:
  plt.show()
except Exception as e:
  warning("Cannot show figure. Error msg: '%s'" % e)

# Alternatively, since we have an analytic solution, we can check the
# :math:`L_2` norm of the error in the solution::

f.interpolate(cos(x*pi*2)*cos(y*pi*2))
print(sqrt(assemble(dot(u - f, u - f) * dx)))

# A python script version of this demo can be found :demo:`here <helmholtz.py>`.
