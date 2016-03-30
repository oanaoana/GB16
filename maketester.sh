#!/bin/bash

ftn -free -O3 nekbox_libxsmm_bench.f -Ilibxsmm/include -Llibxsmm/lib/ -lxsmm -lxsmmf -mkl=sequential -Wl,--whole-archive,-ldmapp,--no-whole-archive -o nektester
