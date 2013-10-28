# See www.zlib.net
include ../../make.inc
MW=../../mwrap

gzmex:
	$(MW) -mex gzmex -mb gzfile.mw
	$(MW) -mex gzmex -c gzmex.c gzfile.mw
	$(MEX) gzmex.c -lz

clean:
	rm -f gz*.m gzmex.* *.o* eye.gz
