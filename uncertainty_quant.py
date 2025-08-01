from sgp4.api import Satrec
import numpy as n
from sgp4.api import WGS72OLD, WGS72, WGS84
from sgp4.api import jday
import matplotlib.pyplot as plt

def get_satellite(jdepoch,bstar=3.8792e-05):
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
        0.0007417,            # ecco: eccentricity
        0.3083420829620822,   # argpo: argument of perigee (radians 0..2pi)
        0.9013560935706996,   # inclo: inclination (radians 0..pi)
        1.4946964807494398,   # mo: mean anomaly (radians 0..2pi)
        0.06763602333248933,  # no_kozai: mean motion (radians/minute)
        3.686137125541276,    # nodeo: R.A. of ascending node (radians 0..2pi)
    )
    return(satellite2)

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
# satellite position
#print(r)

#satellite.bstar=satellite.bstar*0.9
# satellite position
#print(r)