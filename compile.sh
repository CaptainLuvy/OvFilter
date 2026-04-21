#!/bin/bash
# ???? ovfilter.cpp
# ????????? g++ ?? zlib-devel (Ubuntu/Debian) ?? zlib-devel (CentOS/RHEL)

echo "Compiling ovfilter.cpp..."
g++ -O3 -std=c++11 -pthread -march=native ovfilter.cpp -o ovfilter_cpp -lz

if [ $? -eq 0 ]; then
    echo "Compilation successful! Executable generated: ./ovfilter_cpp"
else
    echo "Compilation failed. Please make sure g++ and zlib are installed."
    echo "Ubuntu/Debian: sudo apt-get install g++ zlib1g-dev"
    echo "CentOS/RHEL: sudo yum install gcc-c++ zlib-devel"
fi
