print('hmmmmm')

from ctypes import *

print('loading library')
cdll.LoadLibrary("libkeytransfer.so")
print('loaded library')

libkey = CDLL('libkeytransfer.so')

libkey.make.restype = c_void_p
libkey.fullprove.restype = c_char_p
print('calling make')
handle = libkey.make(b"/usr/local/share/keytransfer/keytransfer.zkey", b"/usr/local/share/keytransfer/keytransfer.dat")
#print('calling fullprove')
#res = libkey.fullprove(handle, b"/tmp/tmp.wtns", b"input.json")
#print(res)
