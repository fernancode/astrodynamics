import astrodynamics as astro

#problem 1
radius_earth = 6378
apogee = 6000 + radius_earth
perigee = 740 + radius_earth
T = astro.Period(apogee,perigee) / 60**2

print('1).')
print('Period is ', T, 'hours.')

e = (apogee - perigee) / (apogee + apogee)
v_perigee = ((perigee*(1+e)*astro.mu)**.5)/perigee
v_apogee =  ((apogee*(1-e)*astro.mu)**.5)/apogee

print('Velocity at perigee is ', v_perigee)
print('Velocity at apogee is ', v_apogee)

a = (apogee + perigee)  /2
h = (a*astro.mu*(1-e**2))**.5

thetas = astro.np.linspace(0,astro.np.pi,1000)

gamma=[]
for theta in thetas:
    r =  h**2 / astro.mu / (1 + e*astro.np.cos(theta))
    v = h/r
    print(v)
    #gamma = astro.np.arctan( (e*astro.np.sin(theta))/(1+ e* astro.np.cos(theta)))
#    print(gamma * 180 / astro.np.pi)


#    r =  h**2 / astro.mu / (1 + e*astro.np.cos(theta))
#    v = h/r
#    v_theta = astro.mu/h * (1+e*astro.np.cos(theta))
#    print(v)
#    print(v_theta)
#    print('\n')
    
#max_angle = gamma.index(max(gamma))
#print(thetas[max_angle])