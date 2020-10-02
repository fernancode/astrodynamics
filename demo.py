import astrodynamics as astro

angular_momentum = 80000
eccentricity = .5

#determine the apogee
apogee = astro.orbit_eq(angular_momentum,eccentricity,3.14159)
print(apogee)

#make a graphic of the orbit
astro.plot_orbit(angular_momentum,eccentricity)
astro.plot_orbit()

