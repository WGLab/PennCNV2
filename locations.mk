# USER NEEDS TO SET THE PATHS HERE

# download Boost from http://www.boost.org
# for root account, just do a 'yum install boost-devel' (if yum is the package manager)
# if you want to install locally, unpack to /home/kaiwang/usr/boost/boost_1_58_0/, then do a './bootstrap.sh --prefix=/home/kaiwang/usr/boost/boost', then do a './b2 install'

BOOST_INC_FLAGS = -I/home/kaiwang/usr/boost/boost_1_58_0/
BOOST_LIB_FLAGS = -L/home/kaiwang/usr/boost/boost/lib

# download GSL from http://www.gnu.org/s/gsl
# for root account, you can just do 'yum install gsl-devel' (if yum is the package manager)
# if you want to install locally, unpack the downloaded package, then do a './configure --prefix=/home/kaiwang/usr/gsl/gsl', then do a 'make' and 'make install'

GSL_INC_FLAGS = -I/home/kaiwang/usr/gsl/gsl/include
GSL_LIB_FLAGS = -lgsl -lgslcblas -L/home/kaiwang/usr/gsl/gsl/lib
