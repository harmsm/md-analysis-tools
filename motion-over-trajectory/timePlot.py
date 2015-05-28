__description__ = \
"""
Takes a pdb file with multiple snapshots from a trajectory and colors by time
point using user-specified color gradients. Can take multiple trajectories 
simultaneously.
"""
__author__ = "Michael J. Harms"
__date__ = "110223"
__usage__ = "pymol -a *.pdb timePlot.py [pdb files *must* have .pdb extension]"

import cmd, sys

# This is a set of in-line functions used to specify color gradients.  fx is 
# the fractional time of a given step within the trajectories (ranging from 0
# to 1). "color_graidents" is a list of references to these functions.  
# If multiple trajectories are loaded, the gradients are applied in the order
# specified in color_gradients.  The trajectories are ordered by where they 
# occur in the command line.
def a(fx): return [fx,0,(1.0-fx)]   #fx,fx,1-fx]
def b(fx): return [1,fx,1]   #1,fx,fx]
def c(fx): return [1,1,fx]   #1-fx,1-fx,fx]

color_gradients = [a,b,c]

# Default representations
cmd.hide("everything","all")
cmd.show("ribbon","all")
cmd.show("spheres","resn N06 or resn N15 or resn E40 or resn EST")

# Create a list of unique trajectories loaded by parsing the command line
unique = [arg[:-4] for arg in sys.argv if arg.endswith(".pdb")]

# Go through each trajecotry
models = cmd.get_names("all")
for i, u in enumerate(unique):

    traj = [m for m in models if "_".join(m.split("_")[:-1]) == u]

    # Go through each step in the trajectory
    num_steps = float(len(traj))
    for j, t in enumerate(traj):

        # Determine rgb of new color
        fx = j/num_steps
        color_list = color_gradients[i](fx) 
    
        # Create a new color
        color_name = "col%s%i" % (u,j)
        cmd.set_color(color_name,color_list)

        # Apply the color to step j of trajectory u
        cmd.color(color_name,t)
    
