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