# lazy-batched-ksw

A Lazy and Batched Key-Switching Framework for Accelerating Computations in Homomorphic Encryption

## Overview

This repository aims to provide a lazy and batched key-switching framework designed to accelerate homomorphic encryption (HE) computations.

> 🚧 Implementation in progress. Core features will be added soon.

## Requirements

- C++ Compiler (g++ ≥ 9.4.0)
- CMake (≥ 3.5)
- Git
- OpenFHE library
## Installing OpenFHE

This project includes a fork of the OpenFHE library as a Git submodule under `openfhe-development/`.  
The lazy-batched key-switching framework is implemented directly within this submodule.

To build and install it:

```bash
cd openfhe-development
mkdir build && cd build
cmake ..
make -j
sudo make install
```

To verify that OpenFHE (including the lazy-batched key-switching extension) was built correctly, you can run:

```bash
make testall
```

You can also try running a sample program:

```bash
bin/examples/pke/simple-integers
```

To clean the build artifacts later:

```bash
make clean
```

> ✅ You do **not** need to clone OpenFHE separately — the required version is already included and contains the key-switching code.



## Clone the Repository

```bash
git clone --recursive https://github.com/oksuman/lazy-batched-ksw.git
cd lazy-batched-ksw
```

## License

This project is licensed under the [MIT License](./LICENSE).
