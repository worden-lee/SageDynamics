# to install the python package into sage.  possibly needs to be run sudo.
# you may prefer to use "setup.py develop" if you'll be editing in these files
install:
	sage -python setup.py install

test:
	sage -t ./dynamicalsystems
