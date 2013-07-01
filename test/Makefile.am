## Process this file with automake to produce Makefile.in
DIST_SUBDIRS = data

INCLUDE_PATTERN = -I%
LIBNAM_PATTERN = -l%

AM_CFLAGS = -I../src/

if ICC_COMPILER
if ICC_NOLIMITH
AM_CFLAGS += -D_GCC_LIMITS_H_
endif
else
AM_CFLAGS += -Wall -Wextra
endif

if APPLE_64BIT
AM_CFLAGS += -m64 -D_FILE_OFFSET_BITS=64
APPLELDFLAGS = -pthread
endif

if SSE2
if ICC_COMPILER
SSE_INTRINS_FLAGS = -xO
else
SSE_INTRINS_FLAGS = -msse2
endif
else
SSE_INTRINS_FLAGS =
endif

SEQ_LIB = ../src/libseq.a
SGM_LIB = ../src/libsgm.a ${SEQ_LIB}

ALI_SRC = ../src/score.c ../src/score.h ../src/diffstr.c ../src/diffstr.h \
	  ../src/alibuffer.c ../src/alibuffer.h ../src/alibuffer_struct.h \
	  ../src/alignment.c ../src/alignment.h
if SSE2
ALI_SRC += ../src/swsimd.c ../src/swsimd.h
endif

HIT_SRC = ../src/sort.c ../src/sort.h \
	../src/interval.c ../src/interval.h \
	../src/hashidx.c ../src/hashidx.h \
	../src/hashhit.c ../src/hashhit.h

SGM_SRC = ${HIT_SRC} ../src/array.c ../src/array.h ../src/vdef.h \
         ../src/segment.c ../src/segment.h 

MAP_SRC = ../src/results.c ../src/results.h \
	../src/resultpairs.c ../src/resultpairs.h \
	../src/insert.c ../src/insert.h \
	../src/report.c ../src/report.h ../src/rmap.c ../src/rmap.h 


SMALT_SRC = ../src/smalt.c ../src/menu.c ../src/menu.h \
	../src/threads.c ../src/threads.h \
	../src/infmt.c ../src/infmt.h \
	${MAP_SRC}

SMALT_TEST_SRC = ${SMALT_SRC} ${SGM_SRC} ${ALI_SRC} tracker.c tracker.h

TEST_DRIVERS = sequence_test.py consistentFR_test.py \
		diffstr_test.py \
		splitReads_test.py xali_test.py \
		results_split_test.py \
		ioform_test.py ouform_cigar_test.py \
		sample_test.py \
		cigar_test.py \
		mthread_test.py
                #regression_test.py #alignment_test.py

if WITH_TESTDRIVER
TESTS_ENVIRONMENT = $(PYTHON)


TESTS = ${TEST_DRIVERS}

if WITH_SAMTOOLS_EXEC
TESTS += bam_cigar_test.py
endif

else
TESTS =
endif

## distribute test drivers
dist_data_DATA = ${TEST_DRIVERS} bam_cigar_test.py formats.py testdata.py

check_PROGRAMS = sequenceReadWrite_test sequenceCompress_test \
		sequenceDeCompress_test \
		sequenceReverseComplement_test \
		diffstr_test \
		smalt_Xali_test

EXTRA_PROGRAMS = seqSetIO_test seqSetCompress_test seqSetSort_test \
                hashhit_test hashhit_track_test hashidx_test \
		segment_test segment_simple_test \
	        array_test diffstr_slice_test bam_test\
		candseg_test \
		alignment_test \
	        rmap_test smalt_test\
		smalt_basq basqual_test \
	 	infmt_test

if WITH_SAMTOOLS
check_PROGRAMS += samtools_api_test
endif

if WITH_BAMBAMC
check_PROGRAMS += bambam_test
endif

if WITH_ZLIB
EXTRA_PROGRAMS += zlib_test gzip_test
ZLIB_LDFLAG= -lz
else
ZLIB_LDFLAG=
endif

if WITH_TESTALL
if SSE2	
check_PROGRAMS += mergesort_test
endif
check_PROGRAMS += ${EXTRA_PROGRAMS}
endif

## test routines for array.c
nodist_array_test_SOURCES = array_test.c
array_test_LDADD = $(SGM_LIB)
array_test_LDFLAGS = -lm $(ZLIB_LDFLAG)

## test routines for diffstr.c
nodist_diffstr_slice_test_SOURCES = diffstr_slice_test.c ../src/diffstr.c ../src/diffstr.h ../src/elib.c ../src/elib.h
diffstr_test_SOURCES = diffstr_test.c ../src/diffstr.c ../src/diffstr.h ../src/elib.c ../src/elib.h

## test routines for BAM file format
nodist_bam_test_SOURCES = bam_test.c ../src/elib.c ../src/elib.h

## merge sort using SIMD
if SSE2
nodist_mergesort_test_SOURCES = mergesort_test.c ../src/elib.c ../src/elib.h ../src/randef.h ../src/sortdef.h \
	../src/sort.c ../src/sort.h
mergesort_test_CFLAGS = $(AM_CLFAGS) $(SSE_INTRINS_FLAGS)
endif

## test routines for sequence.c
sequenceReadWrite_test_SOURCES = sequenceReadWrite_test.c
sequenceReadWrite_test_LDADD = $(SEQ_LIB) $(ZLIB_LDFLAG)

sequenceReverseComplement_test_SOURCES = sequenceReverseComplement_test.c
sequenceReverseComplement_test_LDADD = $(SEQ_LIB) $(ZLIB_LDFLAG)

sequenceCompress_test_SOURCES = sequenceCompress_test.c
sequenceCompress_test_LDADD = $(SEQ_LIB) $(ZLIB_LDFLAG)

sequenceDeCompress_test_SOURCES = sequenceDeCompress_test.c
sequenceDeCompress_test_LDADD = $(SEQ_LIB) $(ZLIB_LDFLAG)

nodist_seqSetIO_test_SOURCES = seqSetIO_test.c
seqSetIO_test_LDADD = $(SEQ_LIB) $(ZLIB_LDFLAG)

nodist_seqSetCompress_test_SOURCES = seqSetCompress_test.c
seqSetCompress_test_LDADD = $(SEQ_LIB) $(ZLIB_LDFLAG)

nodist_seqSetSort_test_SOURCES = seqSetSort_test.c
seqSetSort_test_LDADD = $(SEQ_LIB) $(ZLIB_LDFLAG)

# test routines for hashing
nodist_hashhit_test_SOURCES = hashhit_test.c
hashhit_test_LDADD = $(SGM_LIB) -lm $(ZLIB_LDFLAG)

nodist_hashhit_track_test_SOURCES = hashhit_track_test.c ${HIT_SRC} \
	tracker.c tracker.h \
	../src/diffstr.c ../src/diffstr.h
hashhit_track_test_CFLAGS= -I../test -I../src
hashhit_track_test_CPPFLAGS = -DRESULTS_TRACKER
hashhit_track_test_LDADD = $(SEQ_LIB) -lm $(ZLIB_LDFLAG)

nodist_hashidx_test_SOURCES = hashidx_test.c
hashidx_test_LDADD = $(SGM_LIB) $(ZLIB_LDFLAG)

# tests for candidate segment selection
nodist_segment_test_SOURCES = segment_test.c
segment_test_LDADD = $(SGM_LIB) -lm $(ZLIB_LDFLAG)

nodist_segment_simple_test_SOURCES = segment_simple_test.c
segment_simple_test_LDADD = $(SGM_LIB) -lm $(ZLIB_LDFLAG)

# tests for candidate segment filters
nodist_candseg_test_SOURCES = candseg_test.c ../src/candseg.c ../src/candseg.h ${SGM_SRC}\
	../src/diffstr.c ../src/diffstr.h tracker.c tracker.h
candseg_test_CFLAGS =  $(AM_CFLAGS) -I../test
candseg_test_CPPFLAGS = -DCANDSEG_TRACKER -DRESULTS_TRACKER
candseg_test_LDADD = $(SEQ_LIB) -lm $(ZLIB_LDFLAG)

# tests of alignment routines
nodist_alignment_test_SOURCES = alignment_test.c ${ALI_SRC}
alignment_test_CFLAGS = $(AM_CFLAGS) $(SSE_INTRINS_FLAGS)
alignment_test_CPPFLAGS = -Dalignment_debug -Dalignment_timing
alignment_test_LDADD = $(SEQ_LIB) -lm $(ZLIB_LDFLAG)


# test of mapping routines
nodist_rmap_test_SOURCES = rmap_test.c ${MAP_SRC}
#if WITH_BAMBAMC
# Dummy C++ source to trigger C++ linking.
#nodist_EXTRA_rmap_test_SOURCES = dummy.cxx
#endif
rmap_test_CFLAGS = $(AM_CFLAGS) $(BAMBAMC_INCL)
rmap_test_LDADD = $(SGM_LIB) $(BAMBAMC_DEPS) -lm $(ZLIB_LDFLAG)
rmap_test_LDFLAGS = $(BAMBAMC_LDFLG)

# test splitting of alignments across neighbouring reference sequences
smalt_Xali_test_SOURCES = ${SMALT_SRC}
smalt_Xali_test_CFLAGS = $(AM_CFLAGS) $(SSE_INTRINS_FLAGS) $(BAMBAMC_INCL)
smalt_Xali_test_CPPFLAGS = -DSMALT_DEBUG_XALI
smalt_Xali_test_LDADD = $(SGM_LIB) $(BAMBAMC_DEPS) -lm -lpthread $(ZLIB_LDFLAG)
smalt_Xali_test_LDFLAGS =  $(APPLELDFLAGS) $(BAMBAMC_LDFLG)

# testing smalt with simulated reads (not distributed)
nodist_smalt_test_SOURCES = ${SMALT_TEST_SRC}
smalt_test_CFLAGS = $(AM_CFLAGS) $(SSE_INTRINS_FLAGS) -I../test $(BAMBAMC_INCL)
smalt_test_CPPFLAGS = -DRESULTS_TRACKER

smalt_test_LDADD = $(SEQ_LIB) $(BAMBAMC_DEPS) -lm -lpthread $(ZLIB_LDFLAG)
smalt_test_LDFLAGS = $(APPLELDFLAGS) $(BAMBAMC_LDFLG)

# output bas quality stats with smalt
nodist_smalt_basq_SOURCES = ../src/smalt.c ../src/menu.c ../src/menu.h \
		../src/infmt.c ../src/infmt.h\
		../src/threads.c ../src/threads.h ${MAP_SRC}
smalt_basq_CFLAGS = $(AM_CFLAGS) $(BAMBAMC_INCL)
smalt_basq_CPPFLAGS = -DRESULTS_BASQUAL
smalt_basq_LDADD = $(SGM_LIB) $(BAMBAMC_DEPS) -lm -lpthread $(ZLIB_LDFLAG)
smalt_basq_LDFLAGS = $(APPLELDFLAGS) $(BAMBAMC_LDFLG) 

# extracting quality values for statistical analysis
nodist_basqual_test_SOURCES = basqual_test.c ../misc/basqual.c ../misc/basqual.h
basqual_test_CFLAGS = $(AM_CFLAGS) 
basqual_test_LDADD = $(SEQ_LIB) $(ZLIB_LDFLAG)

# generic interface to sequence input formats
nodist_infmt_test_SOURCES = infmt_test.c ../src/infmt.c ../src/infmt.h
# conditional dummy C++ source to cause C++ linking.
#if WITH_BAMBAMC
#nodist_EXTRA_infmt_test_SOURCES = dummy.cxx
#endif
infmt_test_CFLAGS = $(AM_CFLAGS) $(BAMBAMC_INCL)
infmt_test_LDFLAGS = $(BAMBAMC_LDFLG)
infmt_test_LDADD = $(SEQ_LIB) $(BAMBAMC_DEPS) $(ZLIB_LDFLAG)

# optional samtools api
if WITH_SAMTOOLS
samtools_api_test_SOURCES = samtools_api_test.c ../src/elib.c ../src/elib.h
samtools_api_test_CFLAGS = $(AM_CFLAGS) $(SAMTOOLS_CFLAGS)
samtools_api_test_LDFLAGS = $(SAMTOOLS_LDFLAGS)
samtools_api_test_LDADD =  $(SAMTOOLS_DEPS)
endif

# optional libary for reading/writing BAM/SAM by German Tischler
# has to be linked with g++
if WITH_BAMBAMC
bambam_test_SOURCES = bambam_test.c ../src/elib.c ../src/elib.h
# Dummy C++ source to cause C++ linking.
#nodist_EXTRA_bambam_test_SOURCES = dummy.cxx
bambam_test_CFLAGS = $(AM_CFLAGS) $(BAMBAMC_INCL)
bambam_test_CPPFLAGS = -DRESULTS_BASQUAL
bambam_test_LDFLAGS = $(BAMBAMC_LDFLG)
bambam_test_LDADD = $(BAMBAMC_DEPS)
#
#bambam_test_CFLAGS = -I../src $(filter-out -std=%, $(BAMBAMC_CFLAGS))
#bambam_test_LDFLAGS = $(BAMBAMC_LIBS)
endif

# zlib test
if WITH_ZLIB
nodist_zlib_test_SOURCES = zlib_test.c ../src/elib.c ../src/elib.h
zlib_test_LDADD = -lz

# gzip functions in zlib
nodist_gzip_test_SOURCES = gzip_test.c ../src/elib.c ../src/elib.h
gzip_test_LDADD = -lz
endif
