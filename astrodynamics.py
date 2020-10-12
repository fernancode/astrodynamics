import numpy as np
from scipy.optimize import root_scalar
from scipy.optimize import fsolve
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


def radius_from_dtheta(theta=[], r0=[], v0=[], mu=mu):
    h = np.linalg.norm(np.cross(r0,v0))
    r1 = np.linalg.norm(r0)
    vr0 = np.dot(r0, v0/r1)
    r2 = (h**2)/mu * 1/(1 + (h**2 / (mu*r1)-1)*np.cos(theta) - (h*vr0)/mu * np.sin(theta))
    return r1, r2


def Period(apogee,perigee,mu=mu):
    '''
    Given apogee, perigee, orbital parameter, and radius of body, calculate period of orbit in seconds
    
    **apogee and perigee are measured from center of earth
    '''
    a = (apogee+perigee)/2
    T = (2*np.pi) / mu**.5 * ((a))**1.5
    return T


def to_radian(deg):
    r = deg * np.pi/180
    return r


def to_degrees(rad):
    d = rad * 180/np.pi
    return d


def get_eccentricity(h=[], r0=[], theta=[], mu=mu):
    e = (h**2 / mu - r0)/ (r0*np.cos(theta))
    return e


def lagrange(delta_theta=[], r=[], r0=[], h=[], mu=mu):
    '''
    Calculate the Lagrange coefficients given change in theta, radius, initial radius, angular momentum, and gravitational constant
    '''
    def lagrange_f(delta_theta=[], r=[], h=[], mu=mu):
        f = 1- (r*mu)/h**2 * (1-np.cos(delta_theta))
        return f

    def lagrange_g(delta_theta=[], r=[], r0=[], h=[], mu=mu):
        g = (r*r0) / h * np.sin(delta_theta)
        return g

    def lagrange_fdot(delta_theta=[], r=[], r0=[], h=[], mu=mu):
        fdot = mu/h * (1-np.cos(delta_theta))/np.sin(delta_theta) * (mu/h**2 *(1-np.cos(delta_theta)) -1/r0 - 1/r)
        return fdot

    def lagrange_gdot(delta_theta=[], r0=[], h=[], mu=mu):
        gdot = 1 - (mu*r0)/h**2 * (1-np.cos(delta_theta))
        return gdot
    
    f = lagrange_f(delta_theta, r, h)
    g = lagrange_g(delta_theta, r, r0, h)
    fdot = lagrange_fdot(delta_theta, r, r0, h)
    gdot = lagrange_gdot(delta_theta, r0, h)
    return f,g,fdot,gdot


def esc_velocity(r=[], mu=mu):
    v_esc = (2*mu / r)**.5
    return v_esc


def chobotov_approx(r0=[], dt=[], mu=mu):
    '''
    Chobotov approximation for the universal variable approach. First guess at x0
    '''
    r = np.linalg.norm(r0)
    x0 = (mu**.5)/np.linalg.norm(r) * dt
    return x0


def stumpf_c(z=[]):
    c = (np.cosh((-z)**.5) -1) / -z
    return c


def stumpf_s(z=[]):
    s = (np.sinh((-z)**.5) - (-z)**.5) / ((-z)**.5)**3
    return np.abs(s)


def universal_parameters(r0=[], v0=[], mu=mu):
    r = np.linalg.norm(r0)
    v = np.linalg.norm(v0)
    alpha = 2/r - (v**2)/mu
    a = 1/alpha
    h = np.linalg.norm(np.cross(r0,v0))
    theta = np.arctan(r0[1]/r0[0])
    e = (np.linalg.norm(h)**2 / mu - r) / (np.linalg.norm(r0) * np.cos(theta))
    return alpha, a, e, h, theta


def universal_variable(r0=[], v0=[], alpha=[], dt=[], mu=mu):
    '''

    '''
    
    #get the first approximation
    x0 = chobotov_approx(r0, dt)
    
    #define the function to be solved
    def zero(x, r0=[], v0=[], alpha=[], dt=[], mu=mu):
        r = np.linalg.norm(r0)
        vr0 = np.dot(r0,v0)/r
        z = alpha*x**2
    
        t1 = r*vr0 /(mu**.5) * x**2 * stumpf_c(z)
        t2 = (1-alpha*r)*x**3 * stumpf_s(z)
        t3 = r*x
        t4 = mu**.5 * dt
        return (t1+t2+t3 - t4)
    
    #iterate for  the solution
    sol = fsolve(zero, args=(r0, v0, alpha, dt), x0=x0)
    x = sol[0]
    return x


def universal_lagrange(x=[], r0=[], v0=[], alpha=[], dt=[], mu=mu):
    z = alpha*x**2
    r = np.linalg.norm(r0)
    #v = np.linalg.norm(v0)
    z = alpha*x**2

    f = 1 - x**2/r * stumpf_c(z)
    g = dt - 1/mu**.5 * x**3 * stumpf_s(z)
    r_new = f*r0 + g*v0

    fdot = mu**.5 / (r*np.linalg.norm(r_new)) * (alpha*x**3 * stumpf_s(z) - x)
    gdot = 1 - x**2/np.linalg.norm(r_new) * stumpf_c(z)
    v_new = fdot*r0 + gdot*v0
    return r_new, v_new


def hyperbolic_anomaly(r0=[], v0=[], e=[] , dt=[], mu=mu):
    '''

    '''
    #create the function to solve
    def zero(F, r0=[], v0=[], e=[] , dt=[], mu=mu):
        h = np.linalg.norm(np.cross(r0,v0))
        r = np.linalg.norm(r0)
        a = r / (1-e)

        Mh = (mu/(-a**3))**.5 * dt
        return Mh - (e*np.sinh(F) -F)
    
    #and get the first guess
    x0 = zero(0,r0,v0,e,dt)

    #solve hyperbolic anomaly for F
    sol = fsolve(zero, args=(r0,v0,e,dt), x0=x0)
    F = sol[0]
    #get anomaly
    theta = 2*np.arctan( ((e+1) / (e-1))**.5 * np.tanh(F/2))
    return theta


#TODO: maybe change some of this later.
#Stuff i may change later
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
