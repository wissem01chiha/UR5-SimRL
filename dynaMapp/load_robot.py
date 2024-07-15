import os
from urdfpy import URDF

def check_urdf_file(robot_urdf_file):
    """
    Check if the given file path points to a valid URDF file.
    """
    if os.path.isfile(robot_urdf_file):
        try:
            robot = URDF.load(robot_urdf_file)
            return True
        except Exception as e:
            print(f"Error parsing URDF file: {e}")
            return False
    else:
        print("Invalid URDF file path.")

def robot_joint_num(robot_urdf_file):
    """
    Compute the number of joints of a robot from URDF file.
    """
    if check_urdf_file(robot_urdf_file):
        robot = URDF.load(robot_urdf_file)
        return len(robot.joints)
    return 0

def robot_link_num(robot_urdf_file):
    """ 
    Compute the number of links of a robot from URDF file.
    """
    if check_urdf_file(robot_urdf_file):
        robot = URDF.load(robot_urdf_file)
        return len(robot.links)
    return 0


def link_rotation_matrix(robot, link_number):
    """
    Compute the 3*3 orthoognal roation matrix of 
    """
    return 0
 
 
 
 
 
 
robot = URDF.load("C:/Users/chiha/OneDrive/Bureau/Dynamium/DynaMapp/robots/kinovaGen3.urdf")

print(type(robot.links[1].inertial.inertia))