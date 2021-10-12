
#include "capi.hpp"

int main() {
    void *ptr = make("keytransfer.zkey", "keytransfer.dat");
    fullprove(ptr, "tmp.wtns", "input.json");
}
