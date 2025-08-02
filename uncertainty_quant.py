from sgp4.api import Satrec
import numpy as n
from sgp4.api import WGS72OLD, WGS72, WGS84
from sgp4.api import jday
import scipy.optimize as sio
import scipy.stats as sstat
import matplotlib.pyplot as plt

import numpy as np
from scipy.stats import multivariate_normal
from scipy.integrate import nquad

def collision_probability_3d(miss_vector, sigma_vector, R):
    """
    probability of collision. assumes that at closes approach the two satellites are predicted to be 
    within miss_vector. The uncertainty in the position (combined for both objects) is given by the standard deviation vector
    R is the sum of the radii of the two satellites. If the realized miss distance is less than R, a collision occurs. 
    Returns:
    - Pc: probability of collision (unitless, 0 to 1)
    """
    mean = miss_vector
    cov = np.diag([sigma_vector[0]**2, sigma_vector[1]**2, sigma_vector[2]**2])
    rv = multivariate_normal(mean=mean, cov=cov)

    # Integration over the sphere of radius R centered at origin
    def integrand(z, y, x):
        return rv.pdf([x, y, z])

    def bounds_x():
        return [-R, R]
    
    def bounds_y(x):
        return [-np.sqrt(R**2 - x**2), np.sqrt(R**2 - x**2)]
    
    def bounds_z(x, y):
        return [-np.sqrt(R**2 - x**2 - y**2), np.sqrt(R**2 - x**2 - y**2)]

    Pc, _ = nquad(integrand, [bounds_z, bounds_y, bounds_x])
    return Pc

def get_satellite(jdepoch,
                  bstar=3.8792e-05, 
                  e=0.0007417, 
                  aop=0.3083420829620822,
                  inc=0.9013560935706996, 
                  mu=1.4946964807494398,
                  no_kozai=0.06763602333248933,
                  nodeo=3.686137125541276):
    # TLE epoch
    jd0, fr0 = jday(1949, 12, 31, 0, 0, 0)

    satellite2 = Satrec()
    satellite2.sgp4init(
        WGS72,                # gravity model
        'i',                  # 'a' = old AFSPC mode, 'i' = improved mode
        25544,                # satnum: Satellite number
        jdepoch-jd0,       # epoch: days since 1949 December 31 00:00 UT
        bstar,           # bstar: drag coefficient (1/earth radii)
        0.0,                  # ndot: ballistic coefficient (radians/minute^2)
        0.0,                  # nddot: mean motion 2nd derivative (radians/minute^3)
        e,            # ecco: eccentricity
        aop,   # argpo: argument of perigee (radians 0..2pi)
        inc,   # inclo: inclination (radians 0..pi)
        mu,   # mo: mean anomaly (radians 0..2pi)
        no_kozai,  # no_kozai: mean motion (radians/minute)
        nodeo,    # nodeo: R.A. of ascending node (radians 0..2pi)
    )
    return(satellite2)

def drag_uncertainty_sweep():
    # JD epoch
    epoch=2458826.86525

    # evaluate these epochs (2 days into future)
    jdtimes = n.linspace(epoch,epoch+2,num=1000)
    jdfr = n.zeros(1000)#linspace(epoch,epoch+2,num=1000)

    # evaluate these uncertainties in atmospheric drag
    uncertainties=[0,10,20,30,40,50,60,70,80,90,100]
    for u in uncertainties:
        satellite=get_satellite(jdepoch=epoch,bstar=3.8792e-05)
        satellite2=get_satellite(jdepoch=epoch,bstar=3.8792e-05*(u/100 + 1.0))

        pos_diff=[]
        for i in range(len(jdtimes)):
            e, r, v = satellite.sgp4(jdtimes[i], 0)
            e, r2, v = satellite2.sgp4(jdtimes[i], 0)
            pos_diff.append(n.linalg.norm(n.array(r)-n.array(r2))*1e3)

        plt.plot( (jdtimes-epoch)*24, pos_diff,label="Drag uncertainty=%1.0f %%"%(u))
    plt.xlabel("Time since epoch (hours)")
    plt.ylabel("Position uncertainty (m)")
    plt.legend()
    plt.show()

def cartesian2tle(pos,vel, jdepoch, bstar=3.8792e-05,
                  guess=n.array([0.0007417,0.3083420829620822,0.9013560935706996,1.4946964807494398,0.06763602333248933,3.686137125541276])):
    """ 
    go from position and velocity in TEME at epoch jdepoch to sgp4 format keplerian elements
    we need to fit.
    note that this fit would need a number of positions along an orbit to be stable. now it doesn't really work
    """
    def ss(x):
        sat=get_satellite(jdepoch,
                          bstar=bstar, 
                          e=x[0], aop=x[1],inc=x[2], mu=x[3], no_kozai=x[4], nodeo=x[5])
        e, r, v = sat.sgp4(jdepoch, 0)
        return(n.sum(n.abs(n.array(r)-pos)**2) + n.sum(n.abs(n.array(v)-vel)**2))
    x=sio.fmin(ss,guess)
    #print(x)
    sat=get_satellite(jdepoch,bstar=bstar, 
                      e=x[0], aop=x[1],inc=x[2], mu=x[3], no_kozai=x[4], nodeo=x[5])
    return(sat)

def collision_probability(r0=10, r1=1, pstd0=200, pstd1=200, miss_distance=10):
    """
    collision probability 
    radiuses of two spacecraft and standard deviations of positions
    assume symmetry and normal distributions
    """
    return(collision_probability_3d(miss_vector=n.array([0,0,miss_distance]), 
                                    sigma_vector=n.array([n.sqrt(pstd0**2.0+pstd1**2.0),n.sqrt(pstd0**2.0+pstd1**2.0),n.sqrt(pstd0**2.0+pstd1**2.0)]), 
                                    R=r0+r1))


def drag_uncertainty_pc_sweep(miss_distance=200, rad0=10, rad1=10, initial_uncertainty=200.0):
    # JD epoch
    epoch=2458826.86525

    # evaluate these epochs (2 days into future)
    jdtimes = n.linspace(epoch,epoch+4,num=10)
    jdfr = n.zeros(1000)#linspace(epoch,epoch+2,num=1000)

    # evaluate these uncertainties in atmospheric drag
    uncertainties=[10,50,100]#,20,30,40,50,60,70,80,90,100]
    for u in uncertainties:
        satellite=get_satellite(jdepoch=epoch,bstar=3.8792e-05)
        satellite2=get_satellite(jdepoch=epoch,bstar=3.8792e-05*(u/100 + 1.0))

        pos_diff=[]
        pcs=[]

        for i in range(len(jdtimes)):
            print(i)
            e, r, v = satellite.sgp4(jdtimes[i], 0)
            e, r2, v = satellite2.sgp4(jdtimes[i], 0)
            pos_err=n.sqrt( (n.linalg.norm(n.array(r)-n.array(r2))*1e3)**2.0 + initial_uncertainty**2.0)
            pos_diff.append(pos_err)#n.linalg.norm(n.array(r)-n.array(r2))*1e3)
            print(pos_err)
            # TBD. we should evaluate the position errros of both objects!
            pcs.append(collision_probability(pstd0=pos_err,pstd1=pos_err,miss_distance=miss_distance,r0=rad0,r1=rad1))
        pcs=n.array(pcs)
        plt.semilogy( (jdtimes-epoch)*24, pcs*1e6,label="Drag uncertainty=%1.0f %%"%(u))
    plt.xlabel("Time since epoch (hours)")
    plt.ylabel("Collision probability (ppm)")
    plt.legend()
    plt.show()

if __name__ == "__main__":
    drag_uncertainty_sweep()#miss_distance=200, r0=10, r1=10)
    drag_uncertainty_pc_sweep(miss_distance=200, rad0=10, rad1=10)
#    print(collision_probability(pstd0=200,pstd1=200))
 #   print(collision_probability(pstd0=300,pstd1=300))
  #  print(collision_probability(pstd0=400,pstd1=400))
   ##print(collision_probability(pstd0=600,pstd1=600))

    #drag_uncertainty_sweep()
# satellite position
#print(r)

#satellite.bstar=satellite.bstar*0.9
# satellite position
#print(r)