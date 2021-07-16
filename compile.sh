#!/bin/bash
#compiler=clang++
#compiler=icpc
compiler=g++
<<<<<<< HEAD
#type=Debug
type=Release
#mkdir build
=======
type=Debug
#type=Release
mkdir build
>>>>>>> dev_WDM
cd build
if [ -n "$1" ]
  then
  if [ $1 = "-n" ] || [ $1 = "--new" ]; then
    make clean
    rm -rf CMakeFiles #force to rerun cmake configuration
    cmake -DCMAKE_BUILD_TYPE=$type -DCMAKE_CXX_COMPILER=$compiler ../src
  fi
fi
make -j
make install
cd -
