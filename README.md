# BayesOmics
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)



BayesOmics is now ready to use. See the [wiki](kk) page for further documentation. For comments and questions, please contact shouye.liu@uq.edu.au.


## Installation

On Linux/Unix, to build and make the test:

    git clone https://github.com/ShouyeLiu/BayesOmics.git
    mkdir BayesOmics/build && cd $_
    cmake ..
    make -j10

By default, the project will be built in RELEASE mode, use

    cmake .. -DCMAKE_BUILD_TYPE=DEBUG

to build in DEBUG mode.

## test
Run all tests

    ctest




## Authors

Shouye Liu (University of Queensland)

