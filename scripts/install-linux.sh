#!/bin/sh

if ! grep adx /proc/cpuinfo > /dev/null
then
    git clone https://github.com/mrsmkl/snark-tools
    cd snark-tools
    npm i -g
    cd ..
    sudo apt-get -y install g++ build-essential
    sudo mkdir -p /usr/local/lib && sudo mkdir -p /usr/local/include && sudo mkdir -p /usr/local/share/keytransfer
    g++ -shared -o alt/libkeytransfer.so -fPIC alt/capi.cpp
    sudo cp alt/libkeytransfer.so /usr/local/lib/libkeytransfer.so
    sudo cp keytransfer.zkey keytransfer.dat keytransfer.wasm /usr/local/share/keytransfer
    sudo ldconfig
else
    sudo apt-get -y install g++ nasm libgmp-dev libsodium-dev build-essential
    sudo mkdir -p /usr/local/lib && sudo mkdir -p /usr/local/include && sudo mkdir -p /usr/local/share/keytransfer
    git clone https://github.com/brainhub/SHA3IUF
    cd SHA3IUF && make && sudo make install && cd ..
    rm -rf SHA3IUF
    git submodule update --init --recursive
    npm i 
    npx task createFieldSources
    cd build
    sed -i 's/extern Fq_fail/;extern Fq_fail/g' *.asm
    sed -i 's/extern Fr_fail/;extern Fr_fail/g' *.asm
    sed -i 's/call    Fq_fail/;call    Fq_fail/g' *.asm
    sed -i 's/call    Fr_fail/;call    Fr_fail/g' *.asm
    nasm -felf64 fq.asm
    nasm -felf64 fr.asm
    cd ..
    npx task buildApi
    sudo cp build/libkeytransfer.so /usr/local/lib/libkeytransfer.so
    sudo cp keytransfer.zkey keytransfer.dat keytransfer.wasm /usr/local/share/keytransfer
    sudo ldconfig

fi