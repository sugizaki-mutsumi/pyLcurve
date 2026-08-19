
all : build

MODULE_NAME = lxdata

build:
	python setup.py build_ext --inplace

clean:
	rm -rf build __pycache__ *.pyc *.so ${MODULE_NAME}_wrap.cpp ${MODULE_NAME}.py


