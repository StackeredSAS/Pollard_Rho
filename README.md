# Pollard_Rho

A GMP implementation of the parallelized Pollard's Rho algorithm.

The python wrapper manages the C workers and computes the private key once a collision is found.

Each worker walks its own pseudo-random sequence until it finds a distinguished point $X_i = a_iG + b_iQ$. Once, one is found, it is printed to STDOUT for the wrapper to collect. The wrapper checks if this point has already been found by another worker and if so computes the private key (only if the coefficients $a_i$ and $b_i$ are different). Otherwise, a new worker is instanciated. The amount of workers in parallel is equal to the CPU count of the machine the wrapper is running on.

A point is considered distinguished when its 24 least significant bits are 0. The distinguisher can be chosen arbitrarily but distinguished points shouldn't be too rare or too common.

## Requirements

You will need to install the GMP developer tools before building. This can be done like this on Ubuntu :

```bash
sudo apt install libgmp-dev
```

## Usage

Edit `pollard_rho_worker.h` with your curve parameters and point coordinates, then :

```
make
./wrapper.py
```

## References

[Parallel Collision Search with Cryptanalytic Applications](https://people.scs.carleton.ca/~paulv/papers/JoC97.pdf)
