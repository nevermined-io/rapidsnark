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

RUN git clone     https://github.com/nevermined-io/rapidsnark
RUN cd rapidsnark \
 && git submodule update --init --recursive \
 && sh ./scripts/install-linux.sh
