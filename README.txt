particle.py - Contains the Particles class and methods to update the accelerations, velocities, positions, linear momenta, angular momenta, potential energies and kinetic energies of the particles.
Requires numpy, astropy.

main.py - Creates instances of the Particles class, advances the system in time and saves results to a designated .npy file.
Requires numpy, astropy, poliastro, spiceypy.

analysis.py - Reads the .npy file obtained by running main.py and uses the data to plot graphs as required.
Requires numpy, matplotlib.

main_twobody.py - Copy of main.py edited to work more efficently for the simulation of a two body system.
Requires numpy, astropy, poliastro, spiceypy.

analysis_twobody.py - Copy of analysis.py edited to work more efficently for the analysis of a two body system.
Requires numpy, matplotlib.

In order to run the simulation, main.py should be run within the same directory as particle.py. analysis.py can then be run to import the generated .npy file.