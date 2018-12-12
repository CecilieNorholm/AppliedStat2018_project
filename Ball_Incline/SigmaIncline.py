import numpy as np
import math as m

def SigmaIncline(a,d_rail,D_ball,theta,sigma_a,sigma_rail,sigma_ball,sigma_theta):
    """
    Computes the error on gravitational acceleration g for a
    'Ball Incline Experiment'

    INPUT parameters
    -----------
    a : acceleration array as float
    d_rail : diameter of rail as float
    D_ball : diameter of ball as float
    theta : angle in radians as float
    sigma_a : error on acc as float
    sigma_rail : error on rail as float
    sigma_ball : error on the ball measurement as float
    sigma_theta : error on angle as float


    OUTPUT returns
    ---------
    sigma_g : error on 'g' as a numpy array
    """

    #   rewrtie to radians
    sigma_theta = np.radians(sigma_theta)
    theta = np.radians(theta)


    #Rewriting parameters to ease computation
    delta = D_ball**2 / (D_ball**2 - d_rail**2)

    #derivatives of g w.r.t variable
    diff_theta = 1/5*a*(2*delta+5)*(1./np.tan(theta))*(1./np.sin(theta))
    diff_ball = (-4*a*d_rail**2*D_ball*1./np.sin(theta))/(5*(D_ball**2-d_rail**2)**2)
    diff_rail = 4*a*d_rail*D_ball**2*1./np.sin(theta)/(5*(D_ball**2-d_rail**2)**2)
    diff_a = (2*delta/5+1)*(1./np.sin(theta))

    #calculates error on g
    sigma_g = np.sqrt(2.0*diff_theta**2*sigma_theta**2 + diff_ball**2*sigma_ball**2 \
                     + diff_rail**2*sigma_rail**2 + diff_a**2*sigma_a**2)
    return sigma_g
