
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

struct my_data {
};

extern "C" {
    void *make(const char *_zkey, const char *_dat) {
        my_data *dta = new my_data;
        return (void*)dta;
    }
    const char *fullprove(void *ptr, const char *_wtns, const char *_in) {
        std::cerr << "??? " << _in << "\n";
        std::stringstream s;
        s << "keytransfer-prover " << _in;
        system(s.str().c_str());
        std::ifstream t("/tmp/output.json");
        std::stringstream buffer;
        buffer << t.rdbuf();
        std::string res = buffer.str();
        char *copy = new char[res.size()+1];
        for (int i = 0; i < res.size(); i++) {
            copy[i] = res[i];
        }
        copy[res.size()] = 0;
        return copy;
    }
}
