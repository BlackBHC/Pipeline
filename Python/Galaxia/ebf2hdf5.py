# python 2
import ebf
import h5py
import sys

argv = sys.argv
argv = argv[1:]

f=h5py.File(argv[0]+".hdf5","w")
data = ebf.read(argv[0], '/')

n = 1
while(n < len(argv)):
	f.create_dataset(name=argv[n], data=data[argv[n]])
	n += 1

f.close()
