from ctypes import *

cdll.LoadLibrary("libkeytransfer.so")

libkey = CDLL('libkeytransfer.so')

libkey.make.restype = c_void_p
handle = libkey.make(b"keytransfer.zkey", b"keytransfer.dat")
libkey.fullprove(handle, b"tmp.wtns", b"input.json")

