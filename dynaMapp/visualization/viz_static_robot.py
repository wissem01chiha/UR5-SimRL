"""

"""
import pybullet as p
import time
import random


def main():
    # Connect to the PyBullet physics server
    physicsClient = p.connect(p.GUI)  # Use p.DIRECT for non-graphical version

    # Set gravity
    p.setGravity(0, 0, -9.81)

    # Load the plane (ground)
    #planeId = p.loadURDF("plane.urdf")

    # Load the robot
    robotStartPos = [0, 0, 0]  # Initial position of the robot
    robotStartOrientation = p.getQuaternionFromEuler([0, 1, 0])  # Initial orientation of the robot (no rotation)

    robotId = p.loadURDF("robots/kinovaGen3.urdf", robotStartPos, robotStartOrientation)

    # Set colors for each link
    num_joints = p.getNumJoints(robotId)
    for joint_index in range(num_joints):
        random_color = [random.uniform(0, 1), random.uniform(0, 1), random.uniform(0, 1), 1]
        p.changeVisualShape(robotId, joint_index, rgbaColor=random_color)  # Change color to red (R,G,B,A)
    # Enable moving the cursor
    p.configureDebugVisualizer(p.COV_ENABLE_MOUSE_PICKING, 1)

    # Simulation loop
    for _ in range(5000):  # Simulate for 5000 steps
        p.stepSimulation()  # Advance the simulation by one timestep
        time.sleep(1./300.)  # Control the simulation speed (240 Hz)
        

    # Disconnect from the PyBullet physics server
    p.disconnect()

if __name__ == "__main__":
    main()
