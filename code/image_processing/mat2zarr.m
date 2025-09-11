pyenv('Version','/usr/bin/python')
	
mask = uint32(NMseg); 
py.zarr.save('NMseg.zarr', mask)
	
