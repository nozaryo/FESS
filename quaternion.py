import numpy as np
import numpy.linalg as nplin

#########################################################################
def int_quat(rail_ang, rail_dir):
    """
    Generate initial Quaternion

    """

    ### Compute euler angle
    theta = np.deg2rad(-1 * rail_ang)
    psi = np.deg2rad(rail_dir)

    ### Compute quaternion rotation angle theta
    q_theta = np.arccos(-1*(1-np.cos(theta)-np.cos(psi)-np.cos(theta)*np.cos(psi))/2)

    q_lambda = np.empty(3)
    q_lambda[0] = (-1*np.sin(theta)*np.sin(psi))/(2*np.sin(q_theta))
    q_lambda[1] = (np.sin(theta)*np.cos(psi)+np.sin(theta))/(2*np.sin(q_theta))
    q_lambda[2] = (np.cos(theta)*np.sin(psi)+np.sin(psi))/(2*np.sin(q_theta))

    ### Compute initial quaternion Q
    quat0 = np.empty(4)
    quat0[0:3] = q_lambda * np.sin(q_theta/2)
    quat0[3] = np.cos(q_theta/2)

    return quat0


#########################################################################
def omega2dQuat(omega, quat):

    L_omega = np.array([
                    [0, omega[2], -omega[1], omega[0]],
                    [-omega[2], 0, omega[0], omega[1]],
                    [omega[1], -omega[0], 0, omega[2]],
                    [-omega[0], -omega[1], -omega[2], 0]
                    ])

    dQuat = 0.5 * L_omega @ quat

    return dQuat


#########################################################################
def quat2dcm(quat):
    q1 = quat[0]
    q2 = quat[1]
    q3 = quat[2]
    q4 = quat[3]

    dcm = np.empty([3,3])

    dcm[0,0] = (q1**2) - (q2**2) - (q3**2) + (q4**2)
    dcm[0,1] = 2 * (q1*q2 + q3*q4)
    dcm[0,2] = 2 * (q3*q1 - q2*q4)
    dcm[1,0] = 2 * (q1*q2 - q3*q4)
    dcm[1,1] = (q2**2) - (q3**2) - (q1**2) + (q4**2)
    dcm[1,2] = 2 * (q2*q3 + q1*q4)
    dcm[2,0] = 2 * (q3*q1 + q2*q4)
    dcm[2,1] = 2 * (q2*q3 - q1*q4)
    dcm[2,2] = (q3**2) - (q1**2) - (q2**2) + (q4**2)

    return dcm


#########################################################################
def euler2dcm(phi, theta, psi):
    """
    Use 1-2-3 type (k-j-i type)

    """
    dcm = np.zeros([3,3])
    dcm[0,0] = np.cos(theta) * np.cos(psi)
    dcm[0,1] = np.cos(theta) * np.sin(psi)
    dcm[0,2] = -1.0 * np.sin(theta)

    dcm[1,0] = np.sin(phi) * np.sin(theta) * np.cos(psi) \
               - np.cos(phi) * np.sin(psi)
    dcm[1,1] = np.sin(phi) * np.sin(theta) * np.sin(psi) \
               + np.cos(phi) * np.cos(psi)
    dcm[1,2] = np.sin(phi) * np.cos(theta)

    dcm[2,0] = np.cos(phi) * np.sin(theta) * np.cos(psi) \
               + np.sin(phi) *  np.sin(psi)
    dcm[2,1] = np.cos(phi) * np.sin(theta) * np.sin(psi) \
               - np.sin(phi) * np.cos(psi)
    dcm[2,2] = np.cos(phi) * np.cos(theta)

    return dcm


#########################################################################
def dcm2quat(dcm):
    """ """

    temp = np.zeros(4)
    quat = np.empty(4)

    temp[0] = 0.5 * np.sqrt(1 + dcm[0,0]**2 - dcm[1,1]**2 - dcm[2,2]**2)
    temp[1] = 0.5 * np.sqrt(1 - dcm[0,0]**2 + dcm[1,1]**2 - dcm[2,2]**2)
    temp[2] = 0.5 * np.sqrt(1 - dcm[0,0]**2 - dcm[1,1]**2 + dcm[2,2]**2)
    temp[3] = 0.5 * np.sqrt(1 + dcm[0,0]**2 + dcm[1,1]**2 + dcm[2,2]**2)

    qmax = np.max(temp)

    if (temp[0] == qmax):
        quat[0] = qmax
        quat[1] = 0.25 * (1/qmax) * (dcm[0,1] + dcm[1,0])
        quat[2] = 0.25 * (1/qmax) * (dcm[0,2] + dcm[2,0])
        quat[3] = 0.25 * (1/qmax) * (dcm[1,2] - dcm[2,1])

    elif (temp[1] == qmax):
        quat[0] = 0.25 * (1/qmax) * (dcm[0,1] + dcm[1,0])
        quat[1] = qmax
        quat[2] = 0.25 * (1/qmax) * (dcm[2,1] + dcm[1,2])
        quat[3] = 0.25 * (1/qmax) * (dcm[2,0] - dcm[0,2])

    elif (temp[2] == qmax):
        quat[0] = 0.25 * (1/qmax) * (dcm[2,0] + dcm[0,2])
        quat[1] = 0.25 * (1/qmax) * (dcm[2,1] + dcm[1,2])
        quat[2] = qmax
        quat[3] = 0.25 * (1/qmax) * (dcm[0,1] - dcm[1,0])

    elif (temp[3] == qmax):
        quat[0] = 0.25 * (1/qmax) * (dcm[1,2] - dcm[2,1])
        quat[1] = 0.25 * (1/qmax) * (dcm[2,0] - dcm[0,2])
        quat[2] = 0.25 * (1/qmax) * (dcm[0,1] - dcm[1,0])
        quat[3] = qmax

    return quat


#########################################################################
def quat_normalize(quat):

    quat_norm = nplin.norm(quat)
    quat_n = quat/quat_norm

    return quat_n


#########################################################################
def int_quat_other(rail_ang, rail_dir):
    """
    Generate initial Quaternion

    """

    # Define euler angle
    phi = 0.0
    theta = np.deg2rad(-1.0 * rail_ang)
    psi = np.deg2rad(rail_dir)

    # Convert from Euler to DCM
    dcm = euler2dcm(phi, theta, psi)

    # Convert from DCM to quaternion
    quat = dcm2quat(dcm)
    quat0 = quat_normalize(quat)

    return quat0
