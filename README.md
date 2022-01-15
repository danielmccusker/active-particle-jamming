# active particle simulations
Simulations of Vicsek-like dynamics with repulsive particles, associated with the publication

**Active particle dynamics beyond the jamming density**, Daniel R. McCusker, Ruben van Drongelen, and Timon Idema 2019 [EPL 125 36001](https://dx.doi.org/10.1209/0295-5075/125/36001)

NOTE: This code is © D. McCusker, 2019, and it is made available under the MIT license enclosed with the software. Over and above the legal restrictions imposed by this license, if you use this software for an academic publication then you are obliged to provide proper attribution. 

The code in this repository implements molecular dynamics in two dimensions for overdamped, self-propelled particles using Euler integration and a combined Verlet list/cell list algorithm for tracking neighbors. The code also computes physical quantities of interest: effective diffusion constants, velocity correlations, density fluctuations.

Simulation snapshots are displayed below along with their corresponding locations in phase space.

![alt text](https://github.com/danielmccusker/active-particle-jamming/blob/master/images/snapshots.png?raw=true)

![alt text](https://github.com/danielmccusker/active-particle-jamming/blob/master/images/phase-diagram.png?raw=true)
