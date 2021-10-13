
#ifndef CAPI_H
#define CAPI_H

extern "C" {
    void *make(const char *_zkey, const char *_dat);
    char *fullprove(void *ptr, const char *_wtns, const char *_in);
}

#endif // CAPI_H
