# to install the python package into sage.  possibly needs to be run sudo.
install:
	sage -python setup.py install

# you may prefer to use "setup.py develop" if you'll be editing in these files
install-develop:
	sage -python setup.py develop

test:
	sage -t ./dynamicalsystems
