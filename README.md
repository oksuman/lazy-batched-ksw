# lazy-batched-ksw

A Lazy and Batched Key-Switching Framework for Accelerating Computations in Homomorphic Encryption

## Overview

This repository aims to provide a lazy and batched key-switching framework designed to accelerate homomorphic encryption (HE) computations.


## Clone the Repository

```bash
git clone --recursive -submodules https://github.com/oksuman/lazy-batched-ksw.git
cd lazy-batched-ksw
```


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

> ✅ You do **not** need to clone OpenFHE separately — the required version is already included and contains the key-switching code.


## How to Start 

- The `benchmark` directory contains experiments that compare the runtime of *m* sequential key-switchings with a batched key-switching on an *m*-lazy ciphertext. You can run them as follows:
```
cd \benchamrk
mkdir build
cd build
cmake ..
make -j
./measure-batched-ks
```
or simply run `run_batched_ks_experiments.sh`.

- The `app` directory contains applications built on the lazy batched key-switching framework. You can run them as follows:

```
cd polycircuit
mkdir build
cd build
cmake ..
make -j 
./cifar
```

You can use `OMP_NUM_THREADS=<num_threads>` to set the desired number of threads. Matrix multiplication is similar.

## License

This project is licensed under the [MIT License](./LICENSE).
