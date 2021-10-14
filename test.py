from ctypes import *

cdll.LoadLibrary("libkeytransfer.so")

libkey = CDLL('libkeytransfer.so')

libkey.make.restype = c_void_p
libkey.fullprove.restype = c_char_p
handle = libkey.make(b"keytransfer.zkey", b"keytransfer.dat")
res = libkey.fullprove(handle, b"tmp.wtns", b"input.json")
print(res)
