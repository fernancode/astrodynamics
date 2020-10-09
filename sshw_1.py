import astrodynamics as astro

#problem 1
radius_earth = 6378
apogee = 6000 + radius_earth
perigee = 740 + radius_earth
T = astro.Period(apogee,perigee) / 60**2

print('1).')
print('Period is ', T, 'hours.')

e = (apogee - perigee) / (apogee + apogee)
p = (apogee + perigee)/2 * (1-e)
h = (p * astro.mu) **.5

v_perigee = ((perigee*(1+e)*astro.mu)**.5)/perigee
v_apogee =  ((apogee*(1-e)*astro.mu)**.5)/apogee

print('Velocity at perigee is ', v_perigee)
print('Velocity at apogee is ', v_apogee)

a = (apogee + perigee)  /2
#h = (astro.mu*a*(1-e**2))**.5

 
thetas = astro.np.linspace(0,2*astro.np.pi,10000)
gammas=[]
for theta in thetas:
    r =  h**2 / astro.mu / (1 + e*astro.np.cos(theta))

    sin_gamma = (astro.mu / h * (e*astro.np.sin(theta)))
    cos_gamma = (astro.mu / h * (1+e*astro.np.cos(theta)))
    tan_gamma = sin_gamma/cos_gamma

    gamma = (astro.np.arctan(tan_gamma))
    gammas.append(gamma)
    
max_angle = gammas.index(max(gammas))
print(astro.to_degrees(max(gammas)))
print(astro.to_degrees(thetas[max_angle]))


#problem 3
r0 = astro.np.array([12300, 8000])
v0 = astro.np.array([-1.7, 7.0])
h = astro.np.cross(r0, v0)
theta = astro.np.arctan(r0[1]/r0[0])
v_not = astro.np.linalg.norm(v0)
r_not = astro.np.linalg.norm(r0)

e = (h**2 / astro.mu - r_not) / (r_not*astro.np.cos(theta))
print(e)
dtheta = astro.to_radian(73)
v_r_not = astro.np.dot(r0, v0/r_not)

r = (h**2)/astro.mu * 1/(1 + (h**2 / (astro.mu*r_not)-1)*astro.np.cos(dtheta) - (h*v_r_not)/astro.mu * astro.np.sin(dtheta))
f,g,fdot,gdot = astro.lagrange(dtheta,r,r_not,h,astro.mu)

r_new = f*r0 + g*v0
v_new = fdot*r0 + gdot*v0
print('New radial vector is ', r_new)
print('New veloity vector is ', v_new)