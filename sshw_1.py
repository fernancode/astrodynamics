import astrodynamics as astro

###########################problem 2################################
radius_earth = 6378
apogee = 6000 + radius_earth
perigee = 740 + radius_earth
T = astro.Period(apogee,perigee) / 60**2

e = (apogee - perigee) / (apogee + apogee)
p = (apogee + perigee)/2 * (1-e)
h = (p * astro.mu) **.5

v_perigee = ((perigee*(1+e)*astro.mu)**.5)/perigee
v_apogee =  ((apogee*(1-e)*astro.mu)**.5)/apogee

a = (apogee + perigee)  /2
 
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

print('2).')
print('Period is ', T, 'hours.')
print('Velocity at perigee is ', v_perigee)
print('Velocity at apogee is ', v_apogee)

print('Max flight angle :', astro.to_degrees(max(gammas)))
print('Anomaly at max flight angle :', astro.to_degrees(thetas[max_angle]))
print('\n')

###################problem 3###################################
r0 = astro.np.array([12300, 8000])
v0 = astro.np.array([-1.7, 7.0])
h = astro.np.cross(r0, v0)
theta = astro.np.arctan(r0[1]/r0[0])
v_not = astro.np.linalg.norm(v0)
r_not = astro.np.linalg.norm(r0)
e = (h**2 / astro.mu - r_not) / (r_not*astro.np.cos(theta))

dtheta = astro.to_radian(73)
v_r_not = astro.np.dot(r0, v0/r_not)
r = (h**2)/astro.mu * 1/(1 + (h**2 / (astro.mu*r_not)-1)*astro.np.cos(dtheta) - (h*v_r_not)/astro.mu * astro.np.sin(dtheta))
f,g,fdot,gdot = astro.lagrange(dtheta,r,r_not,h,astro.mu)
r_new = f*r0 + g*v0
v_new = fdot*r0 + gdot*v0
print('Problem 3). ')
print(e)
print('New radial vector is ', r_new)
print('New veloity vector is ', v_new)

print('\n')

################################Problem 4#########################################
r0 = astro.np.array([6500,0,0])
v0 = astro.np.array([0, 1.2 * astro.esc_velocity(6500),0])
dt = 18 * 60**2 #hours to seconds

alpha, a = astro.universal_parameters(r0=r0, v0=v0)
x0 = astro.chobotov_approx(r0, dt)

sol = astro.root_scalar(astro.universal_variable, args=(r0, v0, alpha, dt), x0=x0, x1=2*x0)
x = sol.root
r_new, v_new = astro.universal_lagrange(x, r0, v0, alpha, dt)

r = astro.np.linalg.norm(r_new)

print('Problem 4')
print(r)
print('\n')


###############################Problem 5##################################################
r0 = astro.np.array([10000,6000,800])
v0 = astro.np.array([-2.0, 8.0, 3])
dt = 20*60#minutes to seconds

alpha, a = astro.universal_parameters(r0=r0, v0=v0)
x0 = astro.chobotov_approx(r0, dt)

sol = astro.root_scalar(astro.universal_variable, args=(r0, v0, alpha, dt), x0=x0, x1=2*x0)
x = sol.root
r_new, v_new = astro.universal_lagrange(x, r0, v0, alpha, dt)

print('Problem 5). ')
print(r_new)
print(v_new)


##############################Problem 6####################################################
