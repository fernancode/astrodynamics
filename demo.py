import astrodynamics as astro

#angular_momentum = 80000
#eccentricity = .5

#determine the apogee
#apogee = astro.orbit_eq(angular_momentum,eccentricity,3.14159)
#print(apogee)

#make a graphic of the orbit
#astro.plot_orbit(angular_momentum,eccentricity,22.5)

r0 = astro.np.array((4500, 5500, 3000))
v0 = astro.np.array((-7, .96, 5.0))

h,i,right_asc, e, omega, theta = astro.orbital_elements(r0,v0)
print(' h: %0.2f \n' % h , 'i: %0.2f \n' % i, 'r.a.: %0.2f \n' % right_asc, 'e: %0.4f \n' % e, 'omega: %0.2f \n' % omega, 'theta: %0.2f' % theta)