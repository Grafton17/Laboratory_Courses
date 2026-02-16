# Laboratory_Courses
Here I put some of my (janky) hand written code developed during LCs.

The purpose of this repository is to keep at hand some functions that I have developed personally in the context of the laboratoy courses in my Master degree. Some of these routines already exist in Python, therefore they are just a rough rewriting of more sophisticated and optimized code, while some of them are specifically designed for my lab experiences. 

I have to say, I do not put __ALL__ my code here, but just what would be cumbersome and tedious to rewrite in the possible future. As mentioned in my [README.md document](https://github.com/Grafton17), the purpose of my repositories is just that of note keeping. I do not intend to be extremely systematic and tidy.

# Scope of laboratory courses
A short and not exhaustive list of the objectives of the attended laboratory courses:
- PLASMA PHYSICS LAB I: investigate the formation of travelling linear waves in a simulated plasma system;
- PLASMA PHYSICS LAB I: analyze the turbulent plasma state inside the [Thorello device](https://fusenet.eu/node/517);
- PLASMA PHYSICS LAB II: learn how to use a LaBr3:Ce detector and analyze properly its data;
- PLASMA PHYSICS II: technically not a laboratory course, but the read of Ex 1.7 and 1.8 of Bellan's "Fundamentals of Plasma Physics" sparked the idea of writing a numerical integrator for the analysis of the Rutherford scattering problem;
- COMPUTATIONAL MATERIALS SCIENCE: learn how to write basic Classical Molecular Dynamics (CMD), Kinetic Monte Carlo (KMC) and Metropolis Monte Carlo (MMC) codes.

# Summary

In the [Libraries.py](Libraries.py) file I list all the modules used for the "💾 files 💾". The other codes were independently written (as I said, the repo is quite messy, beacuse I am stil learning how to use Python properly...). Just add the content of this file at the beginning of your executable!

- 💾 [Bi-spectral analysis](Bi-coherence) 💾: bi-spetrum, bi-coherence and summed bi-coherence, useful tools in the search of wave-wave coupling phenomena.
- The 💾 [Conditional Sampling](Conditional-Sampling) 💾 technique, in my case used for the search of coherent structures inside a Simply Magnetized Plasma (SMP). This code proved to be essential in the study of the transport barrier in the GOLEM tokamak when I attended the SUMTRAIC 2025 summer school.
- Statistical evaluation of the 💾 [Dispersion Relation k(f)](Dispersion-Relation) 💾: in case the estimated linear coherence between two signals is way less than one, this approach permits to establish whether a travelling linear wave between two spots exists or not.
- 💾 [Monte Carlo evaluation](MC-estimator) 💾 of a good error estimator for the rescaling of an energy spectrum: when a LaBr3:Ce scintillator is used, its inherent radioactivity counts sum up with those of external sources. In order to remove it, you must subtract from the energy spectrum of interest a 'reference spectrum', acquired previously (maybe with a different acquisition period!). How do you propagate the error from the reference spectrum onto the new one?
- The essence of collisional plasmas is the presence of small angle [Coulomb collisions](Rutherford_scattering), which permit the transfer of energy and momentum between plasma particles. In this code I follow the wise suggestions given by Paul Bellan in his 1.7 axercise and built a numerical integrator for the motion of a charge particle in the presence of variable EM fields, then I use it to simulate the scattering of two charged particles given that the scattering problem can be recast into a one particle problem in the relative frame of motion (exercise 1.8).
- [Classical Molecular Dynamics](Classical_Molecular_Dynamics): using numerical integration we can study the behaviour of Ag atoms arranged in a finite FCC lattice with or without the periodic boundary conditions. Specifically, the LJ 6-12 potential cut with a polynomial junction and the velocity Verlet algorithms were used, and the thermalization and diffusion processes were studied.
- [Kineti Monte Carlo](Kinetic_Monte_Carlo): CMD becomes unsuitable for phenomenon happening on a way bigger time scale than that of vibrations, therefore it is necessary to move into the land of stochastic dynamics and simulate the behaviour of physical systems knowing the probability of something to happen (like a deposition and/or diffusion event). Using a Time-Continuous Markov Chain I study the layer-by-layer and the multi-layer epitaxial growth of Ag on top of Ag.
- [Metropolis Monte Carlo](Metropolis_Monte_Carlo): the evaluation of macroscopical observables is nothing else than the computation of multi-dimensional integrals across an improper subset of the phase space (like the canonical or micro-canonical ensembles). Due to the unreliable application of the ergodic theorem to a simulated trajectory with CMD or KMC, we resolve to sample with an artifically built stochastic process only those microstates that are important enough. In doing so we use the Metropolis acceptance rule in the canonical ensemble.
