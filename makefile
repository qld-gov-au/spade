objects1 = copy.o err.o matrixio.o memory.o vecop.o matop.o pxop.o \
	submat.o init.o otherio.o machine.o matlab.o ivecop.o version.o \
	meminfo.o memstat.o
objects2 = lufactor.o bkpfacto.o chfactor.o qrfactor.o solve.o hsehldr.o \
	givens.o update.o norm.o hessen.o symmeig.o schur.o svd.o fft.o \
	mfunc.o bdfactor.o
HBASE = meschach/err.h meschach/meminfo.h meschach/machine.h meschach/matrix.h
HLIST = $(HBASE) meschach/iter.h meschach/matlab.h meschach/matrix2.h \
	meschach/oldnames.h meschach/sparse.h meschach/sparse2.h \
	meschach/zmatrix.h meschach/zmatrix2.h

spade: spade.o libmeschach.a
	cc spade.o -lmeschach -lm -o spade
spade.o : spade.c
	cc -c spade.c -std=c99
libmeschach.a: $(objects1) $(objects2)
	ar rc libmeschach.a $(objects1) $(objects2)
copy.o : meschach/copy.c $(HBASE)
	cc -c meschach/copy.c
err.o : meschach/err.c $(HBASE)
	cc -c meschach/err.c
matrixio.o : meschach/matrixio.c $(HBASE)
	cc -c meschach/matrixio.c
memory.o : meschach/memory.c $(HBASE)
	cc -c meschach/memory.c
vecop.o : meschach/vecop.c $(HBASE)
	cc -c meschach/vecop.c
matop.o : meschach/matop.c $(HBASE)
	cc -c meschach/matop.c
pxop.o : meschach/pxop.c $(HBASE)
	cc -c meschach/pxop.c
submat.o : meschach/submat.c $(HBASE)
	cc -c meschach/submat.c
init.o : meschach/init.c $(HBASE)
	cc -c meschach/init.c
otherio.o : meschach/otherio.c $(HBASE)
	cc -c meschach/otherio.c
machine.o : meschach/machine.c $(HBASE)
	cc -c meschach/machine.c
matlab.o : meschach/matlab.c $(HBASE)
	cc -c meschach/matlab.c
ivecop.o : meschach/ivecop.c $(HBASE)
	cc -c meschach/ivecop.c
version.o : meschach/version.c $(HBASE)
	cc -c meschach/version.c
meminfo.o : meschach/meminfo.c $(HBASE)
	cc -c meschach/meminfo.c
memstat.o : meschach/memstat.c $(HBASE)
	cc -c meschach/memstat.c
lufactor.o : meschach/lufactor.c $(HLIST)
	cc -c meschach/lufactor.c
bkpfacto.o : meschach/bkpfacto.c $(HLIST)
	cc -c meschach/bkpfacto.c
chfactor.o : meschach/chfactor.c $(HLIST)
	cc -c meschach/chfactor.c
qrfactor.o : meschach/qrfactor.c $(HLIST)
	cc -c meschach/qrfactor.c
solve.o : meschach/solve.c $(HLIST)
	cc -c meschach/solve.c
hsehldr.o : meschach/hsehldr.c $(HLIST)
	cc -c meschach/hsehldr.c
givens.o : meschach/givens.c $(HLIST)
	cc -c meschach/givens.c
update.o : meschach/update.c $(HLIST)
	cc -c meschach/update.c
norm.o : meschach/norm.c $(HLIST)
	cc -c meschach/norm.c
hessen.o : meschach/hessen.c $(HLIST)
	cc -c meschach/hessen.c
symmeig.o : meschach/symmeig.c $(HLIST)
	cc -c meschach/symmeig.c
schur.o : meschach/schur.c $(HLIST)
	cc -c meschach/schur.c
svd.o : meschach/svd.c $(HLIST)
	cc -c meschach/svd.c
fft.o : meschach/fft.c $(HLIST)
	cc -c meschach/fft.c
mfunc.o : meschach/mfunc.c $(HLIST)
	cc -c meschach/mfunc.c
bdfactor.o : meschach/bdfactor.c $(HLIST)
	cc -c meschach/bdfactor.c

.PHONY : clean

clean:
	rm -f $(objects1) $(objects2)
