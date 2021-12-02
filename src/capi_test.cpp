
#include "capi.hpp"
#include <iostream>

int main() {
    void *ptr = make("keytransfer.zkey", "keytransfer.dat");
    std::cerr << "Whats brokne\n";
    char * res = fullprove(ptr, "tmp.wtns", "input2.json");
    std::cerr << "Proof: " << res << "\n";
    res = fullprove(ptr, "tmp.wtns", "input2.json");
    std::cerr << "Proof: " << res << "\n";
}
