# PQC-AKE #

## Overview ##
This repo contains a C++ implementation for the Post-Quantum Lattice-Based Authenticated Key Exchange (AKE) based on the work of [Del Pino, Lyubashevsky, Pointcheval](https://eprint.iacr.org/2016/435.pdf)

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
- KEM.cc - Key Encapsulation Mechanism 
- DS.cc - Digital Signature Scheme (with and without Message Recovery):
- Huffman.cc - Huffman encoding.

## Notes ##
This repo contains the Post-Quantum research work that I am doing for my Master's Thesis for UCSB (2016-2017).
- AKE Folder: Implementing a Lattice-Based AKE based on the work of [Del Pino, Lyubashevsky, Pointcheval](https://eprint.iacr.org/2016/435.pdf)
- MISC Folder: This is a temporary folder where I keep all my WIP, interesting tid-bits of code, or things that didn't work (but I kept anyways to look fondly at in the future).


