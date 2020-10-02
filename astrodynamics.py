import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#find an ode45 python script

'''
numpy vectors recap
[ 1 2 3
 4 5 6] is np.array([1,2,3],[4,5,6])
'''
pi = np.pi
mu = 398600 #km^3/s^2
gc = 6.674e-11 #m^3 / kg s^2

def orbit_eq(angular_momentum, eccentricity, theta, mu=398600):
    """
    Given angular momentum (h), eccentricity of the orbit (apogee - perigee ) / (apogee + perigee), and the angle, we get the radial distance from the planet
    """
    r = (angular_momentum**2/mu)/(1+eccentricity*np.cos(theta))
    return r

def plot_orbit(h,e,phi,mu=398600,planet_radius=6371):
    """
    Make a three dimensional plot of the orbit trajectory given angular momentum h, eccentricity e, inclination phi in degrees, gravitational parameter mu, and planet radius
    
    """
    plt.style.use('dark_background')
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')

    #make_planet(planet_radius,ax)
    x,y,z = generate_orbit(h,e,mu)
    x_i,y_i,z_i = inclination(x,y,z,(-np.pi/180 * phi))
    equal_axes(x_i,y_i,z_i,ax)

    ax.plot(x_i,y_i,z_i,label = "orbit",color = 'r')
    make_planet(planet_radius,ax)
    plt.xlabel('kilometers')       
    ax.legend()
    plt.show()



def make_planet(planet_radius,ax):
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    x = planet_radius * np.outer(np.cos(u), np.sin(v))
    y = planet_radius * np.outer(np.sin(u), np.sin(v))
    z = planet_radius * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, color='b',alpha = .5)

def generate_orbit(h,e,mu):
    x = []
    y = []
    z = []
    thetas = np.linspace(0,2*np.pi,1000)
    for i,theta in enumerate(thetas):
        r = orbit_eq(h,e,theta,mu)
        x.append(r*np.cos(theta))
        y.append(r*np.sin(theta))
        z.append(0)

    coordinate = np.array([[x],[y],[z]])
    return x,y,z

def equal_axes(x,y,z,ax):
    X = np.array(x)
    Y = np.array(y)
    Z = np.array(z)

    # Create cubic bounding box to simulate equal aspect ratio
    max_range = np.array([X.max()-X.min(), Y.max()-Y.min(), Z.max()-Z.min()]).max()
    Xb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][0].flatten() + 0.5*(X.max()+X.min())
    Yb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][1].flatten() + 0.5*(Y.max()+Y.min())
    Zb = 0.5*max_range*np.mgrid[-1:2:2,-1:2:2,-1:2:2][2].flatten() + 0.5*(Z.max()+Z.min())
    # Comment or uncomment following both lines to test the fake bounding box:
    for xb, yb, zb in zip(Xb, Yb, Zb):
        ax.plot([xb], [yb], [zb], 'w')
    plt.grid()

def inclination(x,y,z,phi):
    x_i = []
    y_i = []
    z_i = []

    for a,b,c in zip(x,y,z):
        x_i.append( a*np.cos(phi) + 0*b + c*np.sin(phi))
        y_i.append( a*0 + b*1 + c*0)
        z_i.append( -a*np.sin(phi) + 0*b + c*np.cos(phi))
    
    return x_i,y_i,z_i