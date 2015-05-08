include locations.mk

COMMONCFLAGS = -Wall  -g
COMMONLINKFLAGS = -lm
COMMON_OBJS=analyzer.o main.o io.o utility.o
# compiler
CC = g++
LINKER = g++

OBJS=$(COMMON_OBJS)

PREPROC=
USE_GPU=0
USE_MPI=0
USE_DB=0

CFLAGS=$(COMMONCFLAGS) $(BOOST_INC_FLAGS) $(GSL_INC_FLAGS)
LINKFLAGS=$(COMMONLINKFLAGS) $(BOOST_LIB_FLAGS) $(GSL_LIB_FLAGS)

# below we can dynamically add in the dependent libraries and include flags

ifneq ($(wildcard gwas/*),)
	OBJS+=stepwise.o univariate.o
	PREPROC+=-Dgwas
	USE_DB=1
endif

ifneq ($(wildcard pimsa/*),)
	OBJS+=db_manager.o pathwaysearch.o ramrepository.o mysqlrepository.o
	PREPROC+=-Dpimsa
	USE_DB=1
endif

ifneq ($(wildcard cnv/*),)
	OBJS+=hmm_cnv.o
	PREPROC+=-Dcnv
	USE_GPU=0
endif
	
ifneq ($(wildcard parallel-lasso/*),)
	OBJS+=cross_validation.o power.o stability.o slave_lasso2.o master_lasso2.o
	USE_MPI=1
	USE_GPU=1
	PREPROC+=-Dlasso
endif

ifneq ($(wildcard gpu-impute/*),)
	OBJS+=hmm_impute.o
	USE_MPI=1
	USE_GPU=1
	PREPROC+=-Dgpu_impute
endif

ifneq ($(wildcard ssvs/*),)
	OBJS+=ssvs.o
	PREPROC+=-Dssvs
	USE_DB=1
endif

# now add in the appropriate include and link flags

ifeq ($(USE_GPU),1)
	PREPROC+=-DUSE_GPU
	CFLAGS+=$(OPENCL_INC_FLAGS)
	LINKFLAGS+=$(OPENCL_LIB_FLAGS)
	OBJS+=clsafe.o
endif

ifeq ($(USE_MPI),1)
	PREPROC+=-DUSE_MPI
	CFLAGS+=$(MPI_INC_FLAGS)
	LINKER = mpicxx
endif

ifeq ($(USE_DB),1)
	PREPROC+=-DUSE_DB
	CFLAGS+=$(MYSQL_INC_FLAGS)
	LINKFLAGS+=$(MYSQL_LIB_FLAGS)
endif

PROGRAM=analyzer

analyzer: $(OBJS)
	$(LINKER) -o $(PROGRAM) $(OBJS) $(LINKFLAGS)

main.o: common/main.cpp common/main.hpp
	$(CC) $(CFLAGS) $(PREPROC) -c $<

analyzer.o: common/analyzer.cpp common/analyzer.hpp
	$(CC) $(CFLAGS) $(PREPROC) -c $<

io.o: common/io.cpp common/io.hpp
	$(CC) $(CFLAGS) $(PREPROC) -c $<

utility.o: common/utility.cpp common/utility.hpp
	$(CC) $(CFLAGS) $(PREPROC) -c $<

clsafe.o: common/clsafe.c common/clsafe.h
	$(CC) $(CFLAGS) $(PREPROC) -c $<

#BEGIN COMPILE OF GWAS

stepwise.o: gwas/stepwise.cpp gwas/stepwise.hpp
	$(CC) $(CFLAGS) $(PREPROC) -c $<

univariate.o: gwas/univariate.cpp gwas/univariate.hpp
	$(CC) $(CFLAGS) $(PREPROC) -c $<

#END COMPILE OF GWAS

#BEGIN COMPILE OF TUMOR_CNV

hmm_cnv.o: cnv/hmm_cnv.cpp cnv/hmm_cnv.hpp
	$(CC) $(CFLAGS) $(PREPROC) -c $<

#END COMPILE OF TUMOR_CNV

#BEGIN COMPILE OF PARALLEL-LASSO

master_lasso2.o: parallel-lasso/master_lasso2.cpp parallel-lasso/lasso_mpi2.hpp
	$(CC) $(CFLAGS) $(PREPROC)  -c $<

slave_lasso2.o:  parallel-lasso/slave_lasso2.cpp parallel-lasso/lasso_mpi2.hpp
	$(CC) $(CFLAGS) $(PREPROC)  -c $<

cross_validation.o:  parallel-lasso/cross_validation.cpp parallel-lasso/cross_validation.hpp 
	$(CC) $(CFLAGS) $(PREPROC)  -c $<

power.o:  parallel-lasso/power.cpp parallel-lasso/power.hpp 
	$(CC) $(CFLAGS) $(PREPROC)  -c $<

stability.o:  parallel-lasso/stability.cpp parallel-lasso/stability.hpp 
	$(CC) $(CFLAGS) $(PREPROC)  -c $<

#END COMPILE OF PARALLEL-LASSO

#BEGIN COMPILE OF GPU-IMPUTE

hmm_impute.o: gpu-impute/hmm_impute.cpp gpu-impute/hmm_impute.hpp gpu-impute/hmm_impute_dimensions.h
	$(CC) $(CFLAGS) $(PREPROC) -c $<

#BEGIN COMPILE OF GPU-IMPUTE

#BEGIN COMPILE OF SSVS

ssvs.o: ssvs/ssvs.cpp ssvs/ssvs.hpp ssvs/ssvs_dimensions.h
	$(CC) $(CFLAGS) $(PREPROC) -c $<

#END COMPILE OF SSVS

#BEGIN COMPILE OF PIMSA

pathwaysearch.o :pimsa/pathwaysearch.cpp pimsa/pathwaysearch.hpp
	$(CC) $(CFLAGS) -c $<

db_manager.o: pimsa/db_manager.cpp pimsa/db_manager.hpp
	$(CC) $(CFLAGS) -c $<

ramrepository.o: pimsa/ramrepository.cpp pimsa/db_manager.hpp pimsa/ramrepository.hpp
	$(CC) $(CFLAGS) -c $<

mysqlrepository.o: pimsa/mysqlrepository.cpp pimsa/mysqlrepository.hpp pimsa/db_manager.hpp
	$(CC) $(CFLAGS) -c $<


#END COMPILE OF PIMSA

clean :
	rm -fr *.o $(PROGRAM)
