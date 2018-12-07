import numpy as np
import math as m

def SigmaAcceleration(a,t,d,d_rail,D_ball,theta,Delta_theta,sigma_t,sigma_d,sigma_rail,sigma_ball,sigma_theta):
    """
    Computes the error on gravitational acceleration g for a
    'Ball Incline Experiment'

    INPUT parameters
    -----------
    a : acceleration array of shape [n_samples, 1]
    t : time array of shape [n_samples, 1]
    d : distance between gates array of shape [n_samples, 1]
    d_rail : diameter of rail as float
    D_ball : diameter of ball as float
    theta : angle in radians as float
    Delta_theta : difference angle in radians as float
    sigma_rail : error on rail as float
    sigma_ball : error on the ball measurement as float
    sigma_theta : error on angle as float
    sigma_t : error on time as float
    sigma_d : error on distance between gates as float

    OUTPUT returns
    ---------
    sigma_g : error on 'g' as a numpy array
    """

    # make sure that we have numpy arrays; also
    # reshape the array X to ensure that we have
    # a multidimensional numpy array (ndarray) np.array(a).reshape((a.shape[0], -1))
    a = np.array(a).reshape((len(a),1))
    d = np.array(d).reshape((len(d),1))
    t = np.array(t).reshape((len(t),1))

    #Rewriting parameters to ease computation
    delta = D_ball**2 / (D_ball**2 - d_rail**2)
    phi = Delta_theta

    #computes the error on acceleration
    sigma_a = (2.0 / t**2)**2 * sigma_d**2 + (-4.0*d/t**3)**2 * sigma_t**2

    #derivatives of g w.r.t variable
    diff_theta = 1/5*a*(2*delta+5)*(1./np.tan(theta+phi))*(1./np.sin(theta+phi))
    diff_ball = (-4*a*d_rail**2*D_ball*1./np.sin(theta))/(5*(D_ball**2-d_rail**2)**2)
    diff_rail = 4*a*d_rail*D_ball**2*1./np.sin(theta)/(5*(D_ball**2-d_rail**2)**2)
    diff_a = (2*delta/5+1)*(1./np.sin(theta+phi))

    #calculates error on g
    sigma_g = np.sqrt(2.0*diff_theta**2*sigma_theta**2 + diff_ball**2*sigma_ball**2 \
                     + diff_rail**2*sigma_rail**2 + diff_a**2*sigma_a**2)
    return sigma_g
