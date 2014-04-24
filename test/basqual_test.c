/** Test gathering empirical base quality profile from a set of reads
 * and sampling from empirical distribution for simulated reads.
 */

/****************************************************************************
 ****************************************************************************
 *                                                                          *
 *  Copyright (C) 2010 - 2014 Genome Research Ltd.                          * 
 *                                                                          *        
 *  Author: Hannes Ponstingl (hp3@sanger.ac.uk)                             *
 *                                                                          *
 *  This file is part of SMALT.                                             *
 *                                                                          *
 *  SMALT is free software: you can redistribute it and/or modify it under  *
 *  the terms of the GNU General Public License as published by the Free    *
 *  Software Foundation, either version 3 of the License, or (at your       *
 *  option) any later version.                                              *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful, but     *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of              *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU        *
 *  General Public License for more details.                                *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License along *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.         *
 *                                                                          *
 ****************************************************************************
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#include "elib.h"
#include "sequence.h"
#include "basqual.h"

enum {
  DUMMY_COUNT = 8,
};

const double MAX_TOLERANCE = 0.2;

static int cmpBasQualFreq(const BasQualFreq_t *ap, const BasQualFreq_t *bp)
{
  uint8_t i, a_nq, b_nq, a_qmin, b_qmin;
  uint32_t a_rlen, b_rlen;
  const uint32_t *a_q0p, *b_q0p, *a_qtp, *b_qtp;
  uint64_t j, maxj;

  a_rlen = basQualFreqGetData(ap, &a_q0p, &a_nq, &a_qmin, &a_qtp, NULL);
  b_rlen = basQualFreqGetData(bp, &b_q0p, &b_nq, &b_qmin, &b_qtp, NULL);
  
  if (a_rlen != b_rlen ||
      a_nq != b_nq ||
      a_qmin != b_qmin)
    return ERRCODE_FAILURE;
  
  for (i=0; i<a_nq; i++)
    if (a_q0p[i] != b_q0p[i])
      return ERRCODE_FAILURE;

  maxj = a_nq*a_nq*(a_rlen-1);
  for (j=0; j<maxj; j++)
    if (a_qtp[j] != b_qtp[j])
      return ERRCODE_FAILURE;

  return ERRCODE_SUCCESS;
}

static int cmpBasQualFreqWithTolerance(const BasQualFreq_t *ap, const BasQualFreq_t *bp, double tol)
{
  int i, j, nq;
  uint8_t a_nq, b_nq, a_qmin, b_qmin;
  uint32_t a_rlen, b_rlen, rlen;
  uint32_t r, k, maxk, bs, bt;
  const uint32_t *a_q0p, *b_q0p, *a_qtp, *b_qtp;
  const uint64_t *a_qsp, *b_qsp;
  uint64_t max_bad, bad_ctr, a_tot, b_tot;
  double tol_lo, tol_hi;

  if (tol < 0.001 || tol > 0.999)
    return ERRCODE_ASSERT;
  
  tol_lo = ((double) 1) - tol;
  tol_hi = ((double) 1) + tol;

  a_rlen = basQualFreqGetData(ap, &a_q0p, &a_nq, &a_qmin, &a_qtp, &a_qsp);
  b_rlen = basQualFreqGetData(bp, &b_q0p, &b_nq, &b_qmin, &b_qtp, &b_qsp);
  
  if (a_rlen != b_rlen)
    return ERRCODE_FAILURE;

  if (a_nq != b_nq)
    return ERRCODE_FAILURE;
    
  if (a_qmin != b_qmin)
    return ERRCODE_FAILURE;
  
  nq = a_nq;
  rlen = a_rlen;

  for (i=0; i<nq; i++)
    if (a_q0p[i] + DUMMY_COUNT < b_q0p[i]*tol_lo || 
	a_q0p[i] > b_q0p[i]*tol_hi + DUMMY_COUNT +.5)
      return ERRCODE_FAILURE;
      
  maxk = nq*(rlen-1);
  a_tot = b_tot = 0;
  for (k=0; k<maxk; k++) {
    a_tot += a_qsp[k];
    b_tot += b_qsp[k];
  }
  if (a_tot + DUMMY_COUNT < b_tot*tol_lo ||
      a_tot > b_tot*tol_hi + DUMMY_COUNT)
    return ERRCODE_FAILURE;

  max_bad = maxk*tol;
  bad_ctr=0;
  for (r=1; r<rlen; r++) {
    bs = (r-1)*nq;
    for (i=0; i<nq; i++) {
      bt = bs + i;
      if (a_qsp[bt] + DUMMY_COUNT < b_qsp[bt]*tol_lo || 
	  a_qsp[bt] > b_qsp[bt]*tol_hi + DUMMY_COUNT +.5) {
	printf("BAD_SUM[%llu] (r,i) = (%u, %i) a:%llu <-> b:%llu\n", 
	       (long long unsigned int) bad_ctr, 
	       r, i,
	       (long long unsigned int) a_qsp[bt],
	       (long long unsigned int) b_qsp[bt]);
	bad_ctr++;
	if (bad_ctr > max_bad)
	  return ERRCODE_FAILURE;
      }
    }
  }

  bad_ctr = 0;
  maxk = nq*nq*rlen;
  max_bad = maxk*tol*tol;
  for (r=1; r<rlen; r++) {
    bs = (r-1)*nq;
    for (i=0; i<nq; i++) {
      bt = (bs + i)*nq;
      for (j=0; j<nq; j++) {
	k =  bt + j;
	if (a_qtp[k] + DUMMY_COUNT < b_qtp[k]*tol_lo || 
	    a_qtp[k] > b_qtp[k]*tol_hi + DUMMY_COUNT +.5) {
	  printf("BAD_TRANSITION[%llu] (r, i, j) = (%u, %i, %i) a:%u <-> b:%u\n", 
		 (long long unsigned int) bad_ctr, 
		 r, i, j,
		 a_qtp[k], b_qtp[k]);
	  bad_ctr++;
	  if (bad_ctr > max_bad)
	    return ERRCODE_FAILURE;
	}
      }
    }
  }

  return ERRCODE_SUCCESS;
}

int main(int argc, char *argv[])
{
  int errcode;
  uint8_t maxq, minq;
  uint32_t maxlen, minlen, rlen;
  uint64_t nreads;
  char *infilnam, *oufilnam, *qualp;
  SeqIO *sifp, *sofp;
  SeqFastq *sqbufp;
  BasQualFreq_t *bqfp, *bqfBp;
  ErrMsg *errmsgp=0;
  
  ERRMSG_CREATE(errmsgp);

  if (argc < 3) {
    fprintf(stderr, "usage: %s <FASTQ file> <base quality file (output)>\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  infilnam = argv[1];
  oufilnam = argv[2];
  if (!(sqbufp = seqFastqCreate(0, SEQTYP_FASTQ)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  sifp = seqIOopen(&errcode, infilnam, SEQIO_READ, 0);
  if (errcode) 
    ERRMSGNO(errmsgp, errcode);

  printf("# Find out maximum read length ...\n");
  errcode = basQualFindExtrema(sifp, sqbufp, &nreads, &maxlen, &minlen, &maxq, &minq);
  printf("# Number of reads: %llu\n", (long long unsigned int) nreads);
  printf("# Maximum read length: %u\n", maxlen);
  printf("# Minimum read length: %u\n", minlen);
  printf("# Maximum quality: %hi\n", maxq);
  printf("# Minimum quality: %hi\n", minq);

  if ((errcode = seqIOReset(sifp)))
    ERRMSGNO(errmsgp, errcode);

  if (!(bqfp = basQualFreqCreate(minq, maxq - minq + 1, maxlen)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
 
  printf("\nGetting base qualities ...\n");
  if((errcode = basQualFreqFromFastq(bqfp, sqbufp, sifp)))
    ERRMSGNO(errmsgp, errcode);

  if (seqIOstatus(sifp) && 
      seqIOstatus(sifp) != ERRCODE_EOF) 
    ERRMSGNO(errmsgp, seqIOstatus(sifp));

  printf("\nWriting base qualities to file ...\n");
  if((errcode = basQualFreqWrite(bqfp, oufilnam)))
    ERRMSGNO(errmsgp, errcode);
  
  printf("\nReading base qualities from file ...\n");
  bqfBp = basQualFreqRead(&errcode, oufilnam);
  if ((errcode))
    ERRMSGNO(errmsgp, errcode);
  
  if ((errcode = cmpBasQualFreq(bqfp, bqfBp)))
    ERRMSGNO(errmsgp, errcode);
  
  printf("Test ok.\n");
  basQualFreqDelete(bqfBp);
  bqfBp = NULL;

  printf("\nCalculate sums ...\n");
  if ((errcode = basQualFreqSum(bqfp)))
    ERRMSGNO(errmsgp, errcode);

  printf("\nSample from empirical base quality profile ...\n");
  if ((errcode = seqIOReset(sifp)))
    ERRMSGNO(errmsgp, errcode);

  sofp = seqIOopen(&errcode, oufilnam, SEQIO_WRITE_FASTQ, 0);
  if (errcode) 
    ERRMSGNO(errmsgp, errcode);

  srand(1);
  while(!seqIOstatus(sifp) && !seqIOstatus(sofp)) {
    //printf("Reading sequence %i ...\n", sctr);
    if ((errcode = seqFastqRead(sqbufp, sifp)))
      break;
    qualp = seqFastqGetQualityFactors(sqbufp, &rlen, NULL);

    if ((errcode = basQualFreqSimulate(qualp, rlen, bqfp)))
      ERRMSGNO(errmsgp, errcode);
    
    if ((errcode = seqFastqWrite(sofp, sqbufp, 0)))
      ERRMSGNO(errmsgp, errcode);
  }
  if (seqIOstatus(sofp))
    ERRMSGNO(errmsgp, seqIOstatus(sofp));
  seqIOclose(sofp);

  if (seqIOstatus(sifp) && 
      seqIOstatus(sifp) != ERRCODE_EOF) 
    ERRMSGNO(errmsgp, seqIOstatus(sifp));
  seqIOclose(sifp);

  sifp = seqIOopen(&errcode, oufilnam, SEQIO_READ, 0);
  if (errcode) 
    ERRMSGNO(errmsgp, errcode);

  printf("# Find out maximum read length ...\n");
  errcode = basQualFindExtrema(sifp, sqbufp, &nreads, &maxlen, &minlen, &maxq, &minq);
  printf("# Number of reads: %llu\n", (long long unsigned int) nreads);
  printf("# Maximum read length: %u\n", maxlen);
  printf("# Minimum read length: %u\n", minlen);
  printf("# Maximum quality: %hi\n", maxq);
  printf("# Minimum quality: %hi\n", minq);
  
  if ((errcode = seqIOReset(sifp)))
    ERRMSGNO(errmsgp, errcode);

  if (!(bqfBp = basQualFreqCreate(minq, maxq - minq + 1, maxlen)))
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
 
  printf("\nGetting base qualities ...\n");
  if((errcode = basQualFreqFromFastq(bqfBp, sqbufp, sifp)))
    ERRMSGNO(errmsgp, errcode);

  if (seqIOstatus(sifp) && 
      seqIOstatus(sifp) != ERRCODE_EOF) 
    ERRMSGNO(errmsgp, seqIOstatus(sifp));
  seqIOclose(sifp);

  printf("\nCalculate sums ...\n");
  if ((errcode = basQualFreqSum(bqfBp)))
    ERRMSGNO(errmsgp, errcode);

  printf("\nCompare frequency profiles ...\n");
  if ((errcode = cmpBasQualFreqWithTolerance(bqfp, bqfBp, MAX_TOLERANCE)))
    ERRMSGNO(errmsgp, errcode);
  printf("Test ok.\n");

  basQualFreqDelete(bqfBp);
  basQualFreqDelete(bqfp);
  seqFastqDelete(sqbufp);
  
  ERRMSG_END(errmsgp);
  exit(errcode);
}
