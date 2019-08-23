# tree tensor network strong disorder renormalization group

tSDRG is a powerful and efficient algorithm for disorder system.

This algorithm is based on [PhysRevB.89.214203](https://link.aps.org/doi/10.1103/PhysRevB.89.214203) and [Andrew M. Goldsborough's github](https://github.com/AMGoldsborough/tSDRG) for Matlab, and this code is written in C++ language.

More detail, plz read [tSDRG_algorithm.pdf](./tSDRG_algorithm.pdf)

## Requirements

* C++ compiler (need C++11)
* BLAS and LAPACK libraries and header files
* [Uni10](https://gitlab.com/uni10/uni10) (Note that uni10 is not support function uni10::EigHLazy(), I will upload my uni10 or change EigHLazy() to EigH() or open this issue)

## Installation

```shell
$ git clone https://github.com/ZlinTsai/tSDRG.git
$ cd Main
```

## Compiler command-line by makefile

PLZ chnage makefile UNI10_ROOT and ROOTS to your uni10 directory

UNI10_ROOT    := /usr/local/uni10

ROOTS         := /usr/local/uni10

```shell
make code=yourcode.cpp name=yourcode.exe
```

## Directory

* Operator
    * spin operator

* MPO
    * matrix product operator

* Hamiltonian
    * Hamiltonian model

* tSDRG_tools
    * tSDRG main algorithm and other tools of  tSDRG 
    
* tSDRG_net
    * network file

* Main
    * main.cpp

## Tutorials

```c++
generateTTN(L, chi, Pdis, Jdis, algo, S, Jz, h, Jseed);
```

Coming sooooooooon
