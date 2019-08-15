# tree tensor network strong disorder renormalization group

tSDRG is a powerful and efficient algorithm for disorder system.

This algorithm is based on [PhysRevB.89.214203](https://link.aps.org/doi/10.1103/PhysRevB.89.214203) and [Andrew M. Goldsborough's github](https://github.com/AMGoldsborough/tSDRG) for Matlab, and this code is written in C++ language.

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
## Introduction my code

### Include header file

* Operator.h
    * spin operator

* MPO.h
    * matrix product operator

* Hamiltonian.h
    * Hamiltonian model

* tSDRG_tools.h
    * tSDRG main algorithm and other tools of  tSDRG 


Coming sooooooooon
