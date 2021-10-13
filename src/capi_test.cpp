
#include "capi.hpp"
#include <iostream>

int main() {
    void *ptr = make("keytransfer.zkey", "keytransfer.dat");
    char * res = fullprove(ptr, "tmp.wtns", "input2.json");
    std::cerr << "Proof: " << res << "\n";
    res = fullprove(ptr, "tmp.wtns", "input2.json");
    std::cerr << "Proof: " << res << "\n";
}
