# astrodynamics
Fundamentals of Astrodynamics equations, plots, and visualizations. Code was made from content learned from AE412 Space Systems class at O.S.U. and Fundamentals of Astrodynamics textbook.

Running plot_orbit(h,e,i)
```
plot_orbit(h,e,22.5)
```
where h is angular momentum, e is orbit eccentricity, and i is 22.5 degrees inclination returns a graphic, representing that satelling orbit. Optional arguments are the gravitational parameter mu and the planet radius for non earth systems. 
ex:

```
plot_orbit(h,e,22.5,mu=42830,planet_radius = 3397)
```

![plot_orbit](https://github.com/fernancode/astrodynamics/blob/master/plot_orbit.png)


Running
```
orbit_eq(h,e,theta)
```
where angular momentum, eccentricity, and theta are known, returns radius of the satellite from the body. 
