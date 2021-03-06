import astrodynamics as ad
import numpy as np

###########################problem 2################################
radius_earth = 6371
apogee = 6000 + radius_earth
perigee = 740 + radius_earth
T = ad.Period(apogee,perigee) / 60**2
e = (apogee - perigee) / (apogee + perigee)
p = (apogee + perigee)/2 * (1-e)
h = (p * ad.mu) **.5
v_perigee = ((perigee*(1+e)*ad.mu)**.5)/perigee
v_apogee =  ((apogee*(1-e)*ad.mu)**.5)/apogee
max_angle, max_angle_anomaly =  ad.max_flight_angle(h,e)
max_angle = ad.to_degrees(max_angle)
max_angle_anomaly = ad.to_degrees(max_angle_anomaly)

print('Problem 2).')
print('Period is %0.2f hours.' % T)
print('Velocity at perigee is %0.3f' % v_perigee)
print('Velocity at apogee is %0.3f' % v_apogee)
print('Max flight angle : %0.2f' % max_angle) 
print('Anomaly at max flight angle : %0.2f' % max_angle_anomaly )
print('\n')

###################problem 3###################################
r0 = np.array((12300, 8000))
v0 = np.array((-1.7, 7.0))
alpha,a,e,h,theta = ad.universal_parameters(r0=r0, v0=v0)
theta = ad.to_radian(73)
r1, r2 = ad.radius_from_dtheta(theta=theta, r0=r0, v0=v0)
f,g,fdot,gdot = ad.lagrange(theta,r2,r1,h)
r_new = f*r0 + g*v0
v_new = fdot*r0 + gdot*v0

print('Problem 3). ')
print('e: %0.3f ' % e)
print('New radial vector is [ %0.2f'%r_new[0],', %0.2f]'%r_new[1])
print('New veloity vector is [ %0.2f'%v_new[0],', %0.2f]'%v_new[1])
print('\n')

################################Problem 4#########################################

r0 = np.array((6500,0,0))
v0 = np.array((0, 1.2 * ad.esc_velocity(r0),0))
dt = 18 * 60**2 #hours to seconds
alpha, a, e, h, theta = ad.universal_parameters(r0=r0, v0=v0)
theta = ad.hyperbolic_anomaly(r0,v0,e,dt)
theta_deg = ad.to_degrees(theta)
r = ad.orbit_eq(h, e, theta)

print('Problem 4).')
print('True Anomaly: %0.2f ' % theta_deg)
print('Distance from earth: %0.0f km' % r)
print('\n')



###############################Problem 5##################################################
r0 = ad.np.array((10000,6000,800))
v0 = ad.np.array((-2.0, 8.0, 3))
dt = 20*60

alpha, a, e, h, theta = ad.universal_parameters(r0=r0, v0=v0)
x = ad.universal_anomaly(r0=r0, v0=v0, alpha=alpha, dt=dt)
r_new, v_new = ad.universal_lagrange(x, r0, v0, alpha, dt)
print('Problem 5). ')
print('New radial vector is [ %0.2f'%r_new[0],', %0.2f]'%r_new[1])
print('New veloity vector is [ %0.2f'%v_new[0],', %0.2f]'%v_new[1])

##############################Problem 6####################################################
r0 = np.array((4500, 5500, 3000))
v0 = np.array((-7, .96, 5.0))

sol = ad.numerical_solve(r0=r0, v0=v0, dt=100000, plot=True)
#loop through sol.y's and get magnitude to search for closest approach
