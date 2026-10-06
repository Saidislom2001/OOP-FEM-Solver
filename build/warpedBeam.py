import pyvista as pv
import os

# Automatically find the directory where this Python script is located
script_dir = os.path.dirname(os.path.abspath(__file__))
vtk_path = os.path.join(script_dir, 'cantilever_2d.vtk')

# Load and plot
mesh = pv.read(vtk_path)
warped = mesh.warp_by_vector('Displacement', factor=500)

# Removed the 'title' argument
warped.plot(scalars='Displacement', cmap='jet', show_edges=True)