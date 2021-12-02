FROM ubuntu:21.04

RUN apt-get update -y

RUN apt-get -y install sudo cron curl git g++ nasm libgmp-dev libsodium-dev build-essential

RUN curl -sL https://deb.nodesource.com/setup_14.x | sudo bash -
RUN sudo apt-get -y install nodejs

RUN sudo mkdir -p /usr/local/lib \
 && sudo mkdir -p /usr/local/include \
 && sudo mkdir -p /usr/local/share/keytransfer

RUN git clone https://github.com/brainhub/SHA3IUF
RUN cd SHA3IUF && make && sudo make install

RUN git clone   https://github.com/nevermined-io/rapidsnark
RUN cd rapidsnark \
 && git submodule update --init --recursive \
 && npm i \
 && npx task createFieldSources \
 && cd build \
 && sed -i 's/extern Fq_fail/;extern Fq_fail/g' *.asm \
 && sed -i 's/extern Fr_fail/;extern Fr_fail/g' *.asm \
 && sed -i 's/call    Fq_fail/;call    Fq_fail/g' *.asm \
 && sed -i 's/call    Fr_fail/;call    Fr_fail/g' *.asm \
 && nasm -felf64 fq.asm \
 && nasm -felf64 fr.asm \
 && cd .. \
 && npx task buildApi \
 && sudo cp build/libkeytransfer.so /usr/local/lib \
 && sudo cp keytransfer.zkey keytransfer.dat keytransfer.wasm /usr/local/share/keytransfer \
 && sudo ldconfig \
 && python3.9 test.py

