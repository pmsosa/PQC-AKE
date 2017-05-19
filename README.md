# PQC-AKE #

## Overview ##
This repo contains a C++ implementation for the Post-Quantum Lattice-Based Authenticated Key Exchange (AKE) based on the work of [Del Pino, Lyubashevsky, Pointcheval](https://eprint.iacr.org/2016/435.pdf).

## PreReqs ##
```
sudo apt-get install libgmp-dev libntl-dev libssl-dev
```

## Running the AKE ##
```
cd AKE
make AKE
./AKE
```

## Code Structure ##
This AKE consists of a Key Encapsulation Mechanism (KEM), a Digital Signature (DS) scheme and huffman encoder (for Message-Recovery)
- [KEM.cc](../master/AKE/KEM.cc) - Key Encapsulation Mechanism 
- [DigitalSignature.cc](../master/AKE/DigitalSignature.cc) - Digital Signature Scheme (with and without Message Recovery):
- [Huffman.cc](../master/AKE/huffman.cc) - Huffman encoding.
- [AKE.cc](../master/AKE/AKE.cc) - Two-Way Authenticated Key Exchange
- [Param.h](../master/AKE/params.h) - This file defines parameters n, q, etc...

## Notes ##
This repo contains the Post-Quantum research work that I am doing for my Master's Thesis for UCSB (2016-2017).
- [AKE Folder](../master/AKE): Implementing a Lattice-Based AKE based on the work of [Del Pino, Lyubashevsky, Pointcheval](https://eprint.iacr.org/2016/435.pdf)
- [MISC Folder](../master/Misc): This is a temporary folder where I keep all my WIP, interesting tid-bits of code, or things that didn't work (but I kept anyways to look fondly at in the future).


