## Download Prerequisites

You will need to check that the following libraries are available on your system.

- GNU C compiler (gcc/g++)
- GNU scientific library (http://www.gnu.org/software/gsl/)
- Boost C++ libraries (http://www.boost.org)

These are usually pre-installed on modern linux distributions, and if not, a simple command such as 'yum install gsl-devel' and 'yum install boost-devel' (if YUM is the package manager) will suffice.

For boost, if you do not have root access, or if you want to install locally, unpack to `/home/kaiwang/usr/boost/boost_1_58_0/`, then do a `./bootstrap.sh --prefix=/home/kaiwang/usr/boost/boost`, then do a `./b2 install`

For GSL, if you do not have root access, or if you want to install locally, unpack the downloaded package, then do a `./configure --prefix=/home/kaiwang/usr/gsl/gsl`, then do a `make` and `make install`

## Configure installation

Edit the locations.mk file to specify the correct locations for boost and gsl.

```
BOOST_INC_FLAGS = -I/home/kaiwang/usr/boost/boost_1_58_0/
BOOST_LIB_FLAGS = -L/home/kaiwang/usr/boost/boost/lib

GSL_INC_FLAGS = -I/home/kaiwang/usr/gsl/gsl/include
GSL_LIB_FLAGS = -lgsl -lgslcblas -L/home/kaiwang/usr/gsl/gsl/lib
```

## Building

To compile, run:

```
make
```

You will see that a new `analyzer` program is created in the current directory. Installation is done!
