## Process this file with automake to produce Makefile.in
if ICC_NOLIMITH
AM_CFLAGS = -D_GCC_LIMITS_H_ -I../src
else 
AM_CFLAGS = -I../src/
endif

if WITH_ZLIB
ZLIB_LDFLAG= -lz
else
ZLIB_LDFLAG=
endif

bin_PROGRAMS = simread fetchseq splitreads mixreads \
	splitmates readstats trunkreads basqcol simqual

dist_data_DATA = SAM.py

simread_SOURCES = simread.c rsample.c rsample.h ../src/randef.h
simread_LDADD = ../src/libsgm.a ../src/libseq.a -lm $(ZLIB_LDFLAG)

fetchseq_SOURCES = fetchseq.c 
fetchseq_LDADD = ../src/libseq.a $(ZLIB_LDFLAG)

splitreads_SOURCES = splitreads.c
splitreads_LDADD = ../src/libseq.a $(ZLIB_LDFLAG)

mixreads_SOURCES = mixreads.c
mixreads_LDADD = ../src/libseq.a $(ZLIB_LDFLAG)

readstats_SOURCES = readstats.c
readstats_LDADD = ../src/libseq.a $(ZLIB_LDFLAG)

splitmates_SOURCES = splitmates.c
splitmates_LDADD = ../src/libseq.a $(ZLIB_LDFLAG)

trunkreads_SOURCES = trunkreads.c
trunkreads_LDADD = ../src/libseq.a $(ZLIB_LDFLAG)

basqcol_SOURCES = basqcol.c basqual.c basqual.h
if ICC_NOLIMITH
basqcol_CFLAGS = -D_GCC_LIMITS_H_ -I../src
else
basqcol_CFLAGS = -I../src
endif
basqcol_LDADD = ../src/libseq.a $(ZLIB_LDFLAG)

simqual_SOURCES = simqual.c basqual.c basqual.h ../src/diffstr.c ../src/diffstr.h \
	 	  ../src/randef.h
if ICC_NOLIMITH
simqual_CFLAGS = -D_GCC_LIMITS_H_ -I../src
else
simqual_CFLAGS = -I../src
endif
simqual_LDADD = ../src/libseq.a -lm $(ZLIB_LDFLAG)

