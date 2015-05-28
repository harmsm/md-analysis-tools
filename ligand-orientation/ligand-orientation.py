#!/usr/bin/env python
__description__ =\
"""
"""

__author__ = "Michael J. Harms"
__date__ = "101122"
__usage__ = "calcEuler.py vmd_output_file"

import sys
from math import acos, sqrt, pi, asin, atan2, cos, sin
from numpy import matrix, array, dot, zeros, cross

class CalcEulerError(Exception):
    """
    General error class for this module.
    """
    
    pass

# -----------------------------------------------------------------------------
# Basic geometery functions 
# -----------------------------------------------------------------------------

def scalarAngle(v1,v2):
    """
    Use a dot product to calculate the angle between two vectors.
    """

    l_v1 = sqrt(sum([v1[0,i]**2 for i in range(3)]))
    l_v2 = sqrt(sum([v2[0,i]**2 for i in range(3)]))

    d = dot([v1[0,i] for i in range(3)],
            [v2[0,i] for i in range(3)])

    # The "round" call is a hack to deal with a float round error that
    # accumulates, causing acos to throw a range error with values of
    # 1.00000001 etc.

    return acos(round(d/(l_v1*l_v2),10))



def genRotMatrix(axis,theta):
    """
    Generate a rotation matrix for rotation of theta about axis.
    """

    M = matrix(zeros((3,3)))

    axis_length = sqrt((axis[0,0]**2 + axis[0,1]**2 + axis[0,2]**2))
    xNorm = axis[0,0]/axis_length
    yNorm = axis[0,1]/axis_length
    zNorm = axis[0,2]/axis_length

    sin_theta = sin(theta)
    cos_theta = cos(theta)
    one_costheta = 1.0 - cos_theta

    M[0,0] = cos_theta                + xNorm*xNorm*one_costheta
    M[0,1] = xNorm*yNorm*one_costheta - zNorm*sin_theta
    M[0,2] = xNorm*zNorm*one_costheta + yNorm*sin_theta
    M[1,0] = xNorm*yNorm*one_costheta + zNorm*sin_theta
    M[1,1] = cos_theta                + yNorm*yNorm*one_costheta
    M[1,2] = yNorm*zNorm*one_costheta - xNorm*sin_theta
    M[2,0] = xNorm*zNorm*one_costheta - yNorm*sin_theta
    M[2,1] = yNorm*zNorm*one_costheta + xNorm*sin_theta
    M[2,2] = cos_theta                + zNorm*zNorm*one_costheta

    return M


def arbRotCoord(coord,axis,theta):
    """
    Rotate all vectors in coord about an arbitray axis by theta.
    """

    matrix = genRotMatrix(axis,theta)
  
    return coord*matrix
 


#------------------------------------------------------------------------------
# High level geometry functions
#------------------------------------------------------------------------------

def assignAxes(principal_axes,reference_axes):
    """
    The principal axes are "natural" axes assigned to describe the steroid 
    orientation by vmd.  The reference axes are ~orthoganol vectors connecting
    atoms in the steroid that can be used to assign the arbitrarily assigned
    principal axes to consistent orientations along the steroid.  Equivalent
    reference and principal axes are determined by finding which axes are 
    closest to parallel to one another.  The principal axes are then rotated
    such that principal x aligns to reference x, principal y to reference y,
    and principal z to reference z.  The signs of the axes are then determined
    to make sure the axes not only align, but point in the correct direction.  

    To assign each principal axis to its appropriate reference axis:

    1. Calculate the scalar angle between a principal axis and each reference
       axis.
    2. Subtract pi/2 from each angle.  This will rotate a perfectly 
       parallel axis to -pi/2 and a perfectly anti-parallel axis to pi/2.  All
       other angles will fall somewhere between -pi/2 and pi/2.
    3. Take the absolute value of these rotated angles.
    4. Select reference axis with the maximum value of theta.

    If two axes are assigned the same reference axis, the function dies 
    noisily.
    """
    
    # Figure out the reference axis that is closest to parallel for each 
    # principal axis.  
    assignments = []
    for q in principal_axes:
        theta_list = [] 
        for r in reference_axes:
            theta_list.append(abs(scalarAngle(q,r) - pi/2))

        assignments.append(theta_list.index(max(theta_list)))       

    # Make sure that each axis was assigned to a unique reference axis
    tmp_dict = dict([(a,()) for a in assignments])
    if len(tmp_dict.keys()) != len(principal_axes):
        err = "Could not uniquely assign axes to reference axes!\n"
        raise CalcEulerError(err)

    # Transform the principal axes into the correct orientation relative to the
    # reference axes
    num_rotations = 0
    for i, q in enumerate(principal_axes):

        # Two rotations are sufficient to get us into the correct orientation 
        if num_rotations == 2:
            break

        # If an axis is already lined up, it counts as a rotation
        if i == assignments[i]:
            num_rotations += 1
            continue

        # Figure out which axis we are going to do a rotation about.  (Take the 
        # axis that is not the principal axis or its new assignment).
        to_rotate = [0,1,2]
        to_rotate.remove(i)
        to_rotate.remove(assignments[i])

        # Perform the rotation on the principal axes
        rotation_axis = principal_axes[to_rotate[0]]
        principal_axes = arbRotCoord(principal_axes,rotation_axis,pi/2)

        num_rotations += 1
   
    # Check for handedness of principal_axes.  If the handedness is different
    # than the reference axes, 1 or 3 axes will face in the incorrect direction.
    # If the handedness is correct, 0 or 2 axes will face in the incorrect 
    # direction.  For an incorrect handedness, arbitrarily flip the x-axis.
    theta_list = [scalarAngle(principal_axes[i],reference_axes[i])
                  for i in range(3)]

    num_large_rot = 0
    for t in theta_list:
        if t > pi/2:
            num_large_rot += 1

    #if num_large_rot == 1 or num_large_rot == 3:
    
   #     print "# flipping handedness of principal axes", theta_list

     #   principal_axes[0,0] *= -1
     #   principal_axes[0,1] *= -1
     #   principal_axes[0,2] *= -1

    # Now check to see if we need to do a final rotation to get the signs of the 
    # axes correct.  Because we are only rotating and have a defined handed-ness,
    # the only possible problem is that two axes are pointing the wrong way, one
    # the right way.  If we need to rotate, rotate 180 degrees around the axis
    # that is already facing the correct direction to get everyone lined up.
    theta_list = [scalarAngle(principal_axes[i],reference_axes[i])
                  for i in range(3)]
    num_large_rot = len([t for t in theta_list if t > pi/2])
    if num_large_rot == 2:

        rotation_axis = principal_axes[theta_list.index(min(theta_list))]
        principal_axes = arbRotCoord(principal_axes,rotation_axis,pi) 
        
 
    # Make sure that axes have all been aligned properly and die if this is
    # not true
    check_assignments = []
    for q in principal_axes:
        theta_list = [] 
        for r in reference_axes:
            theta_list.append(abs(scalarAngle(q,r) - pi/2))

        check_assignments.append(theta_list.index(max(theta_list)))       

    # If the axes have been aligned, check_assignments will be [0,1,2]
    if check_assignments != [0,1,2]:
        err = "Problem aligning principal and reference axes!\n"
        raise CalcEulerError(err)

    # Now verify that the axes are all facing in the correct direction
    theta_list = [scalarAngle(principal_axes[i],reference_axes[i])
                  for i in range(3)]

    #if len([t for t in theta_list if t > pi/2]) > 0:
    #    print "#", theta_list
    #   
    #    err = "Axis has the incorrect sign.  Check handedness of reference!\n"
    #    raise CalcEulerError(err)

    return principal_axes


def determineTranslation(translation,reference_orientation):
    """
    Convert a translation in the world reference frame to a translation in
    the steroid reference frame.
    """

    return translation*reference_orientation.I
   
 
def printMatrix(M):
    """
    Print a numpy matrix.
    """

    out = []
    for i in range(3):
        out.append("[")
        for j in range(3):
            out.append("%10.3f" % M[i,j])
        out.append("]\n")
    print "".join(out)    

def eulerFromMatrix(R):
    """
    Determine the euler angles given by a rotation matrix.  By personal 
    convention, I use the euler angles given by "beta", but have left the
    definitions of beta2, alpha2, and gamma2 in for reference.  This function
    assumes that the rotation matrix R is given by Rz(gamma)*Ry(beta)*Rx(alpha)
    (i.e. rotate about x, then y, then z).

    Greg Slabaugh posted psuedocode for this function in a "Computing Euler
    angles from a rotation matrix" (1999) on his website:
    http://www.gregslabaugh.name/publicationsAll.html
    """


    if R[2,0] != 1 and R[2,0] != -1:

        beta = -asin(R[2,0])
        alpha = atan2(R[2,1]/cos(beta),R[2,2]/cos(beta))
        gamma = atan2(R[1,0]/cos(beta),R[0,0]/cos(beta))
        
        beta2 = pi - beta
        alpha2 = atan2(R[2,1]/cos(beta2),R[2,2]/cos(beta2))
        gamma2 = atan2(R[1,0]/cos(beta2),R[0,0]/cos(beta2))
   
    else:

        alpha = 0

        if R[2,0] == -1:
            beta = pi/2
            gamma = alpha + atan2(R[0,1],R[0,2])
        else:
            beta = -pi/2
            gamma = -alpha + atan2(-R[0,1],-R[0,2])

        alpha2 = alpha
        beta2 = beta
        gamma2 = gamma
  
    return (alpha, beta, gamma), (alpha2, beta2, gamma2)


def calcEquivalentEuler(euler):
    """
    Given a set of euler angles, calculate the equivalent (alternate)
    transformation.  
    """

    alpha = euler[0]
    beta = euler[1]
    gamma = euler[2]

    if cos(beta) == 0: 
        return alpha, beta, gamma

    R11 = cos(beta)*cos(gamma)
    R21 = cos(beta)*sin(gamma)
    R32 = sin(alpha)*cos(beta)
    R33 = cos(alpha)*cos(beta)

    beta2 = pi - beta
    alpha2 = atan2(R32/cos(beta2),R33/cos(beta2))
    gamma2 = atan2(R21/cos(beta2),R11/cos(beta2))
  
    return alpha2, beta2, gamma2


def determineRotation(steroid_orientation,reference_orientation,cutoff=1e-2):
    """
    Only does rotation if angle differs significantly from 0 or pi.  This 
    minimizes the number of transformations that are done.  At the end, we 
    then decide the correct way to describe the signs of each flip.  
    """

    # Rotate the steroid coordinates so that they're in the local reference
    # frame (e.g. i = {1,0,0}, j = {0,1,0}, k = {0,0,1}) by multiplying by
    # the inverse of the reference orientation.
    M = steroid_orientation*reference_orientation.I

    # Determine the rotation angles that would have lead to this rotation
    # matrix 
    v1, v2 = eulerFromMatrix(M.I)

    if sum([abs(v) for v in v1]) < sum([abs(v) for v in v2]):
        output = list(v1[:])
    else:
        output = list(v2[:]) 

    # Collapse the range of roll, pitch, and yaw
    for i in range(3):
        output[i] = output[i] % (2*pi)
    
        if output[i] > pi:
            output[i] = output[i] - (2*pi)
        elif output[i] < -pi:
            output[i] = output[i] + (2*pi)

    # Deal with gimbal lock
    real_values = {( 2,-2,-2):pi/2,
                   ( 2, 2, 2):pi/2,
                   (-2,-2,-2):-pi/2,
                   ( 2,-2, 2):-pi/2}              
     
    gimbal_lock = [True for o in output if abs(abs(o)-pi/2) < cutoff]
    if len(gimbal_lock) == 3:
        key = tuple([int(round(o,0)) for o in output])

        output[0] = 0
        output[1] = real_values[key]
        output[2] = 0
    



    return tuple(output)




# -----------------------------------------------------------------------------
# Interface functions
# -----------------------------------------------------------------------------
 
def parseVMDLine(line,reference_orientation=None,min_euler_axis=1):
    """
    Parse a line of output from ligand-orientation.tcl and return either
    the reference axes relating the steroid reference frame to the world
    reference frame (if no reference_orientation is supplied) or the 
    translation and rotations specified in an output line relative to 
    the reference frame.
    """

    # Parse the line
    split_line = line.split("|")
  
    # Frame number 
    frame = int(split_line[1])

    # Pull the principal and internal axes specified on a vmd output line
    # and align the principal axes to the internal axes
    principal_axes = []
    reference_axes = []
    for i in range(3):
        principal_axes.append([float(v) for v in split_line[3+i].split()])
        reference_axes.append([float(v) for v in split_line[6+i].split()])

    principal_axes = matrix(principal_axes)
    reference_axes = matrix(reference_axes)
    
    steroid_orientation = assignAxes(principal_axes,reference_axes)

    # If no reference orientation was supplied, just return the
    # steroid_orientation (This is useful for creating the
    # reference_orientation in the first place). 
    if reference_orientation == None:
        return steroid_orientation

    # Grab the translation vector and use the steroid_orientation to determine
    # how the translation in the world reference frame maps to the tranlsation
    # in the steroid reference frame
    translation = matrix([float(v) for v in split_line[2].split()])
    translation = determineTranslation(translation,reference_orientation)

    # Now calculate the rotation angles describing the position of the
    # steroid_orientation relative to the reference_orientation.  
    rotation = determineRotation(steroid_orientation,reference_orientation)
 
    # Make a vector description of output
    out = [frame]
    out.extend([translation[0,i] for i in range(3)])
    out.extend([e*180/pi for e in rotation])

    return tuple(out)

def main(argv=None):
    """
    Go through the output of ligand-orientation.tcl and convert the raw vectors
    output by VMD into translations and rotation angles within the steroid
    reference frame.  Spit out in R-compatible format.
    """

    # Parse command line
    if argv == None:
        argv = sys.argv[1:]

    try:
        input_file = argv[0]
    except IndexError:
        err = __usage__
        raise CalcEulerError(err)

    # Read the input file
    f = open(input_file,'r')
    lines = f.readlines()
    f.close()

    lines = [l for l in lines if l.startswith("->")]

    # Grab the reference axis from the first line of output
    reference_axes = parseVMDLine(lines[0]) 

    # Parse each line in the file
    out = []
    for l in lines:
        out.append("%10i%10.3f%10.3f%10.3f%10.3f%10.3f%10.3f\n" % 
                   parseVMDLine(l,reference_axes))

    # Add line numbers and header
    out = ["%10i%s" % (i,x) for i, x in enumerate(out)]
    out.insert(0,"%10s%10s%10s%10s%10s%10s%10s%10s\n" % 
               (" ","frame","x","y","z","roll","pitch","yaw"))

    print "".join(out)

 
if __name__ == "__main__":
    
    main()     
