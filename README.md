# Pollard_Rho

A GMP implementation of Pollard's Rho algorithm using Brent's path finding algorithm.

A python wrapper is provided to spread the search over multiple processes to speed things up.

## Usage

Edit `pollard_rho.h` with your curve parameters and point coordinates, then :

```
make
./wrapper.py
```
