import numpy as np


def recursive_newton_euler(robot, q, qdot, qddot, gravity):
    """
    Recursive Newton-Euler Implementation.

    Parameters:
        robot   - Robot structure model
        q       - Joints position vector (1 * n)
        qdot    - Joints velocities vector (1 * n)
        qddot   - Joints accelerations vector (1 * n)
        gravity - Gravity acceleration value

    Returns:
        Q       - The computed torque
    """
    ndof = robot.degrees_of_freedom
    Q = np.zeros(ndof)
    w = np.zeros((3, ndof))
    wdot = np.zeros((3, ndof))
    vdot = np.zeros((3, ndof))
    z0 = np.array([0, 0, 1])

    for i in range(ndof):
        R = robot.rotation_matrix(i)
        p = get_robot_property(robot, "linkFramePosition", i)
        
        if i > 0:
            w[:, i] = np.dot(R.T, w[:, i-1] + z0 * qdot[i])
            wdot[:, i] = np.dot(R.T, wdot[:, i-1] + z0 * qddot[i] + np.cross(w[:, i-1], z0 * qdot[i]))
            vdot[:, i] = np.dot(R.T, vdot[:, i-1]) + np.cross(wdot[:, i], p) + ...
            np.cross(w[:, i], np.cross(w[:, i], p))
        else:
            w[:, i] = np.dot(R.T, z0 * qdot[i])
            wdot[:, i] = np.dot(R.T, z0 * qddot[i])
            vdot[:, i] = np.dot(R.T, gravity) + np.cross(wdot[:, i], p) + ...
            np.cross(w[:, i], np.cross(w[:, i], p))

    n = np.zeros((3, ndof))
    f = np.zeros((3, ndof))

    for i in range(ndof-1, -1, -1):
        p = get_robot_property(robot, "linkFramePosition", i)
        I = get_robot_property(robot, "linkInertiaTensor", i)
        
        vcdot = vdot[:, i] + np.cross(wdot[:, i], robot.link_com_vect(i)) + ...
        np.cross(w[:, i], np.cross(w[:, i], robot.link_com_vect(i)))
        
        F = robot.links[i].mass * vcdot
        N = np.dot(I, wdot[:, i]) + np.cross(w[:, i], np.dot(I, w[:, i]))
        
        if i < ndof - 1:
            R = robot.rotation_matrix(i+1)
            n[:, i] = np.dot(R, n[:, i+1] + np.cross(np.dot(R.T, p), f[:, i+1])) + ...
            np.cross(robot.link_com_vect(i) + p, F) + N
            f[:, i] = np.dot(R, f[:, i+1]) + F
        else:
            n[:, i] = np.cross(robot.link_com_vect(i) + p, F) + N
            f[:, i] = F

        R = robot.rotation_matrix(i)
        if robot.links[i].type == "prismatic":
            Q[i] = np.dot(f[:, i], np.dot(R.T, z0))
        elif robot.links[i].type == "revolute":
            Q[i] = np.dot(n[:, i], np.dot(R.T, z0))

    return Q