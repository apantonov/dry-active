# Description
Simulation of active Ornstein-Uhlenbeck particles with Coulomb (dry) friction

# Citation
If you use our code for your own research, we ask you to cite the following paper:
A. P. Antonov, H. Löwen, Physical Review E 113, 025407 (2026)

# Dependencies
Our software requires boost, which can be downloaded from the following source:
https://www.boost.org/

Random numbers were generated using the Mersenne Twister algorithm implemented in MersenneTwister.h (Matsumoto & Nishimura, 1998).

# How to use
To compile the code, use the provided makefile:

make dryfric

To see available parameters, run ./dryfric --help
