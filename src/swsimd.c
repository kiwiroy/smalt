

/**< Smith-Waterman alignment using SIMD instructions */

/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 Genome Research Ltd.                                  * 
 *                                                                           *        
 *  Author: Hannes Ponstingl (hp3@sanger.ac.uk)                              *
 *                                                                           *
 *  This file is part of SMALT.                                              *
 *                                                                           *
 *  SMALT is free software: you can redistribute it and/or modify it under   *
 *  the terms of the GNU General Public License as published by the Free     *
 *  Software Foundation, either version 3 of the License, or (at your        *
 *  option) any later version.                                               *
 *                                                                           *
 *  This program is distributed in the hope that it will be useful, but      *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of               *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU         *
 *  General Public License for more details.                                 *
 *                                                                           *
 *  You should have received a copy of the GNU General Public License along  *
 *  with this program.  If not, see <http://www.gnu.org/licenses/>.          *
 *                                                                           *
 *****************************************************************************
 *****************************************************************************/

#include <string.h>
#include <limits.h>
#include <emmintrin.h>

#include "swsimd.h"
#include "alibuffer_struct.h"

enum {
  BYTEMASK = 0x00ff,
  BIAS_SHORT = 32768,
  NBITS_PER_BYTE = 8,
  NBYTES_REGISTER = 16, /**< NBYTES_REGISTER * 8 = 128 bits */
};

typedef unsigned char UCHAR;

struct SwsBuff {
  int qlen_max;          /**< Max query length for which buffer is currently initialised */
  UCHAR *datap;          /**< point of memory allocation */
  size_t allocsiz;       /**< number of ALIDPMSCOR_t elements currently 
			  * allocated per buffer */
  int blocksiz;          /**< Block size (granularity) for memory allocation
			  * as the number of ALIDPMSCOR_t elements per buffer. */
  __m128i *H1v;
  __m128i *H2v;
  __m128i *Ev;
};

/******************************************************************************
 ********************************** Macros ************************************
 ******************************************************************************/

#define ZERO_REGISTER_VARIABLE(v) \
  (v) = _mm_setzero_si128()

/* memset to stop compiler/leak checker comlaining about uninitialised values */
/* #define ZERO_REGISTER_VARIABLE(v) \ */
/*   memset(&(v), 0, NBYTES_REGISTER); \ */
/*   (v) = _mm_xor_si128 ((v), (v)) */

#define COPY_SHORT_TO_SSE_VARIABLE(b, v) \
  ZERO_REGISTER_VARIABLE((v)); \
  (v) = _mm_insert_epi16 ((v), b, 0); \
  (v) = _mm_shufflelo_epi16 ((v), 0); \
  (v) = _mm_shuffle_epi32 ((v), 0);

#define ZERO_NEGBIAS_SSE_VARIABLE(v) \
  ZERO_REGISTER_VARIABLE((v)); \
  (v) = _mm_cmpeq_epi16 ((v), (v)); \
  (v) = _mm_slli_epi16 ((v), 15); 

#define COPY_BYTE_TO_REGISTER_VALUE(b, v) \
  ZERO_REGISTER_VARIABLE((v)); \
  intval = ((b) << NBITS_PER_BYTE) | (b & BYTEMASK); \
  (v) = _mm_insert_epi16 ((v), intval, 0); \
  (v) = _mm_shufflelo_epi16 ((v), 0); \
  (v) = _mm_shuffle_epi32 ((v), 0);


/******************************************************************************
 ********************** Private SIMD Alignment Methods ************************
 ******************************************************************************/

#ifdef alignment_matrix_debug
static void printfStripedShortVector(short *vp, unsigned int qlen, int segsiz)
{
  unsigned int j , k, nseg;

  if (segsiz <= 0)
    return;
  nseg = (qlen + segsiz - 1)/segsiz;
  
  for (j=0; j<nseg; j++)
    for (k = (unsigned int) j; k<qlen; k += SCORSIMD_NSHORTS)
      printf("%3hi|", vp[k]+BIAS_SHORT);
      //printf("%3u|", k);
}
static void printfStripedByteVector(unsigned char *vp, unsigned int qlen, int segsiz)
{
  unsigned int j , k, nseg;

  if (segsiz <= 0)
    return;
  nseg = (qlen + segsiz - 1)/segsiz;
  
  for (j=0; j<nseg; j++)
    for (k = (unsigned int) j; k<qlen; k += SCORSIMD_NBYTES)
      printf("%3hu|", (unsigned short) vp[k]);
      //printf("%3u|", k);
}
#endif
static int alignSmiWatShortStriped(unsigned short *maxscor,
			    AliBuffer *abp, 
			    const ScoreProfile *spp,
#ifdef alignment_matrix_debug
			    const char *psqp,
			    const SeqCodec *codecp,
#endif
			    const char *usqp,
			    int uslen)
     /* Striped Smith-Waterman using SSE2 instructions.
      * adapted from M. Farrar (2007) Bioinformatics 23, 156 - 161.
      * http://farrar.michael.googlepages.com/Smith-waterman
      */
{
  int errcode, i, j, segsiz, score, cmpval;
  short tmp;
  unsigned short gap_init, gap_ext;
  const __m128i *vScorep, *vProfp;
  __m128i vH, vE, vF, vMax, vMin, vGapI, vGapE, vTmp, *vp;
  __m128i *vEp = abp->Ev;
  __m128i *vHSp = abp->H1v;
  __m128i *vHLp = abp->H2v;
#ifdef alignment_matrix_debug
  char unprof_symbol;
  const char *decoderp = seqCodecGetDecoder(codecp, NULL);
  unsigned short maximum_score = 0;
  unsigned int q, qlen;
  __m128i vMaxBuf, vTmpBuf;
#endif


  *maxscor = 0;
  if (!((vEp) && (vHSp) && (vHLp)))
    return ERRCODE_ASSERT;

  vProfp = (const __m128i *) scoreGetStripedProfile(NULL, 
#ifdef alignment_matrix_debug
				  &qlen,
#else
				  NULL, 
#endif
				  &gap_init, &gap_ext, NULL, &segsiz,
				  SCORPROF_STRIPED_16, spp);
  if (!vProfp)
    return ERRCODE_SWATSTRIP;

  /* Set SSE constants */
  /* should be able to do this with _mm_set1_epi16(short b) */
  COPY_SHORT_TO_SSE_VARIABLE(gap_init, vGapI);
  COPY_SHORT_TO_SSE_VARIABLE(gap_ext, vGapE);
  
  /*  load vMaxScore with the zeros.  since we are using signed */
  /*  math, we will bias the maxscore to -32768 so we have the */
  /*  full range of the short.
   */
  ZERO_NEGBIAS_SSE_VARIABLE(vMax)

  /* initialize elements of vMin to 4 */
  /* should be able to do this with vMin = _mm_shuffle_epi32 (vGapI, 0); */
  vMin = _mm_shuffle_epi32 (vMax, 0);
  vMin = _mm_srli_si128 (vMin, 14); /* shift right by 14 bytes, fill with 0s */

  /* initialize storage vector to 0 (biased to -32768) */
  for (i=0; i<segsiz; i++) {
    _mm_store_si128 (vEp + i, vMax);
    _mm_store_si128 (vHSp + i, vMax);
  }

  for (i=0; i<uslen; i++) {
    vScorep = vProfp + (usqp[i]&SEQCOD_ALPHA_MASK) * segsiz;
    
    /* zero out F */
    ZERO_NEGBIAS_SSE_VARIABLE(vF);
    
    /* load the next h value */
    vH = _mm_load_si128 (vHSp + segsiz - 1);
    vH = _mm_slli_si128 (vH, 2); /* shift left 2 bytes (for short) */
    vH = _mm_or_si128 (vH, vMin); /* initialise new short on right to -32768 (conceptual 0) */
	
    /* swap the two H vectors */
    vp = vHLp;
    vHLp = vHSp;
    vHSp = vp;

    for (j = 0; j < segsiz; j++) {
      /* load values of vF and vH from previous row (one unit up) */
      vE = _mm_load_si128 (vEp + j);

      /* add score to vH */
      vH = _mm_adds_epi16 (vH, vScorep[j]);

      /* Update highest score encountered this far */
      vMax = _mm_max_epi16 (vMax, vH);

#ifdef alignment_matrix_debug
      /* find largest score in the maxscorv vector */
      vMaxBuf = _mm_shufflehi_epi16 (vMax, 0xE4); /* copy instruction */
      vTmpBuf = _mm_srli_si128 (vMaxBuf, 8);
      vMaxBuf = _mm_max_epi16 (vMaxBuf, vTmpBuf);
      vTmpBuf = _mm_srli_si128 (vMaxBuf, 4);
      vMaxBuf = _mm_max_epi16 (vMaxBuf, vTmpBuf);
      vTmpBuf = _mm_srli_si128 (vMaxBuf, 2);
      vMaxBuf = _mm_max_epi16 (vMaxBuf, vTmpBuf);
      score = _mm_extract_epi16 (vMaxBuf, 0);
      score += BIAS_SHORT;
      //printf("sse_score = %i\n", score);
      if (score > maximum_score) {
	maximum_score = score;
	printf("sse_max_scor(%i,%i) = %i\n", i, j, maximum_score);
      }
#endif

      /* get max from vH, vE and vF */
      vH = _mm_max_epi16 (vH, vE);
      vH = _mm_max_epi16 (vH, vF);

      /* save vH values */
      _mm_store_si128 (vHSp + j, vH);

      /* update vE value */
      vH = _mm_subs_epi16 (vH, vGapI);
      vE = _mm_subs_epi16 (vE, vGapE);
      vE = _mm_max_epi16 (vE, vH);

      /* update vF value */
      vF = _mm_subs_epi16 (vF, vGapE);
      vF = _mm_max_epi16 (vF, vH);

      /* save vE values */
      _mm_store_si128 (vEp + j, vE);

      /* load the next h value */
      vH = _mm_load_si128 (vHLp + j);
    }

    /* reset pointers to the start of the saved data */
    j = 0;
    vH = _mm_load_si128 (vHSp + j);

    /*  the computed vF value is for the given column.  since */
    /*  we are at the end, we need to shift the vF value over */
    /*  to the next column. */
    vF = _mm_slli_si128 (vF, 2);
    vF = _mm_or_si128 (vF, vMin);
    vTmp = _mm_subs_epi16 (vH, vGapI);
    vTmp = _mm_cmpgt_epi16(vF, vTmp);
    cmpval  = _mm_movemask_epi8 (vTmp); /* Creates a 16-bit mask from the most significant bits of 
					 * the 16 signed or unsigned 8-bit integers in vTemp and 
					 * zero extends the upper bits. */

    /* lazy F-loop */
    while (cmpval != 0x0000) { /* not all 8-bit integers == 0 */
      vE = _mm_load_si128 (vEp + j);
      vH = _mm_max_epi16 (vH, vF);

      /* save vH values */
      _mm_store_si128 (vHSp + j, vH);
      
      /*  update vE incase the new vH value would change it */
      vH = _mm_subs_epi16 (vH, vGapI);
      vE = _mm_max_epi16 (vE, vH);
      _mm_store_si128 (vEp + j, vE);

      /* update vF value */
      vF = _mm_subs_epi16 (vF, vGapE);

      j++;
      if (j >= segsiz) {
	j = 0;
	vF = _mm_slli_si128 (vF, 2);
	vF = _mm_or_si128 (vF, vMin);
      }

      vH = _mm_load_si128 (vHSp + j);

      vTmp = _mm_subs_epi16 (vH, vGapI);
      vTmp = _mm_cmpgt_epi16 (vF, vTmp);
      cmpval  = _mm_movemask_epi8 (vTmp);
    }
#ifdef alignment_matrix_debug
    printf("  [%i,%i]: ", i, 0);
    for (q=0; q<qlen; q++)
      printf("%4i", q);
    unprof_symbol = decoderp[(UCHAR) usqp[i]];
    printf("\n%c [%i,%i]: |", unprof_symbol, i, 0);
    for (q=0; q<qlen; q++)
      printf(" %c |", decoderp[(UCHAR) psqp[q]]);
    printf("\nHp[%i,%i]: |", i, 0);
    printfStripedShortVector((short *) vHSp, qlen, segsiz);
    printf("\n");
#endif
  }

  /* find largest score in the maxscorv vector */
  vTmp = _mm_srli_si128 (vMax, 8);
  vMax = _mm_max_epi16 (vMax, vTmp);
  vTmp = _mm_srli_si128 (vMax, 4);
  vMax = _mm_max_epi16 (vMax, vTmp);
  vTmp = _mm_srli_si128 (vMax, 2);
  vMax = _mm_max_epi16 (vMax, vTmp);

  /* return */
  tmp = _mm_extract_epi16 (vMax, 0);
  score = ((int) tmp) + BIAS_SHORT;

  if (score >= USHRT_MAX)
    errcode = ERRCODE_SWATEXCEED;
  else if (score < 0)
    errcode = ERRCODE_ASSERT;
  else {
    errcode = ERRCODE_SUCCESS;
    *maxscor = (unsigned short) score;
  }
  
  return errcode;
}

static int alignSmiWatByteStriped(UCHAR *maxscor,
			   AliBuffer *abp,
			   const ScoreProfile *spp,
#ifdef alignment_matrix_debug
			    const char *psqp,
			    const SeqCodec *codecp,
#endif
			   const char *usqp,
			   int uslen)
     /* Striped Smith-Waterman using SSE2 instructions.
      * adapted from M. Farrar (2007) Bioinformatics 23, 156 - 161.
      * http://farrar.michael.googlepages.com/Smith-waterman
      *
      * \param maxscor Returns maximum score
      * \param abp Buffers used for dynamic programming
      * \param ssp Sequence profile (query).
      * \param psqp Profiled (query) sequence in SEQCOD_MANGLED encoding.
      * \param codec Sequence de/encoder.
      * \param usqp Unprofiled (subject) sequence in SEQCOD_MANGLED encoding.
      * \param uslen Length of the unprofiled (subject) sequence.
      */
{
  int errcode, i, j, score = 0, intval, segsiz, cmpval = 0;
  unsigned short gap_init, gap_ext, bias = 0;
  const __m128i *vScorep, *vProfp;
  __m128i vH, vE, vF, vMax, vBias, vZero, vGapI, vGapE, vTmp, *vp;
  __m128i *vEp = abp->Ev;
  __m128i *vHSp = abp->H1v;
  __m128i *vHLp = abp->H2v;
#ifdef alignment_matrix_debug
  char unprof_symbol;
  const char *decoderp = seqCodecGetDecoder(codecp, NULL);
  int maximum_score = 0;
  unsigned int q, qlen;
  __m128i vMaxBuf, vTmpBuf;
#endif

  *maxscor = 0;
  if (!((vEp) && (vHSp) && (vHLp)))
    return ERRCODE_ASSERT;

  vProfp = (const __m128i *) scoreGetStripedProfile(NULL, 
#ifdef alignment_matrix_debug
				  &qlen,
#else
				  NULL, 
#endif
				  &gap_init, &gap_ext, &bias, &segsiz,
				  SCORPROF_STRIPED_8, spp);
  if (!vProfp)
    return ERRCODE_SWATSTRIP;

  /* Load constants to SSE2 register value */
  /* should be able to do this with _mm_set1_epi8(char b) */
  COPY_BYTE_TO_REGISTER_VALUE(bias, vBias);
  COPY_BYTE_TO_REGISTER_VALUE(gap_init, vGapI);
  COPY_BYTE_TO_REGISTER_VALUE(gap_ext, vGapE);
  
  /* variables initialised to 0 */
  ZERO_REGISTER_VARIABLE(vMax);
  ZERO_REGISTER_VARIABLE(vZero);

  /* Zero out the storage vector */
  for (i=0; i<segsiz; i++) {
    _mm_store_si128 (vEp + i, vZero);
    _mm_store_si128 (vHSp + i, vZero);
  }

  for (i=0; i<uslen; i++) {
    vScorep = vProfp + (usqp[i]&SEQCOD_ALPHA_MASK) * segsiz;
    
    /* zero out F */
    ZERO_REGISTER_VARIABLE(vF);
    
    /* load the next h value */
    vH = _mm_load_si128 (vHSp + segsiz - 1);
    vH = _mm_slli_si128 (vH, 1);
	
    /* swap the two H vectors */
    vp = vHLp;
    vHLp = vHSp;
    vHSp = vp;

    for (j = 0; j < segsiz; j++) {
      /* load values of vF and vH from previous row (one unit up) */
      vE = _mm_load_si128 (vEp + j);

      /* add score to vH */
      vH = _mm_adds_epu8 (vH, vScorep[j]);
      vH = _mm_subs_epu8 (vH, vBias);

      /* Update highest score encountered this far */
      vMax = _mm_max_epu8 (vMax, vH);

#ifdef alignment_matrix_debug
      /* find largest score in the maxscorv vector */
      vMaxBuf = _mm_shufflehi_epi16 (vMax, 0xE4); /* copy instruction */
      vTmpBuf = _mm_srli_si128 (vMaxBuf, 8);
      vMaxBuf = _mm_max_epu8 (vMaxBuf, vTmpBuf);
      vTmpBuf = _mm_srli_si128 (vMaxBuf, 4);
      vMaxBuf = _mm_max_epu8 (vMaxBuf, vTmpBuf);
      vTmpBuf = _mm_srli_si128 (vMaxBuf, 2);
      vMaxBuf = _mm_max_epu8(vMaxBuf, vTmpBuf);
      vTmpBuf = _mm_srli_si128 (vMaxBuf, 1);
      vMaxBuf = _mm_max_epu8 (vMaxBuf, vTmpBuf);
      
      score = _mm_extract_epi16 (vMaxBuf, 0);
      score &= BYTEMASK;
      //printf("sse_score = %i\n", score);
      if (score > maximum_score) {
	maximum_score = score;
	printf("sse_max_scor(%i,%i) = %i\n", i, j, maximum_score);
      }
#endif

      /* get max from vH, vE and vF */
      vH = _mm_max_epu8 (vH, vE);
      vH = _mm_max_epu8 (vH, vF);

      /* save vH values */
      _mm_store_si128 (vHSp + j, vH);

      /* update vE value */
      vH = _mm_subs_epu8 (vH, vGapI);
      vE = _mm_subs_epu8 (vE, vGapE);
      vE = _mm_max_epu8  (vE, vH);

      /* update vF value */
      vF = _mm_subs_epu8 (vF, vGapE);
      vF = _mm_max_epu8 (vF, vH);

      /* save vE values */
      _mm_store_si128 (vEp + j, vE);

      /* load the next h value */
      vH = _mm_load_si128 (vHLp + j);
    } /* for (j = 0; j < segsiz; j++) */

    /* reset pointers to the start of the saved data */
    j = 0;
    vH = _mm_load_si128 (vHSp + j);

    /*  the computed vF value is for the given column.  since */
    /*  we are at the end, we need to shift the vF value over */
    /*  to the next column. */
    vF = _mm_slli_si128 (vF, 1);
    vTmp = _mm_subs_epu8 (vH, vGapI);
    vTmp = _mm_subs_epu8 (vF, vTmp);
    vTmp = _mm_cmpeq_epi8 (vTmp, vZero); /* Where 8-bit integers are equal, all bits are set (0xff)
					  * otherwise unset (0x00). */
    cmpval = _mm_movemask_epi8 (vTmp);  /* Creates a 16-bit mask from the most significant bits of 
					 * the 16 signed or unsigned 8-bit integers in vTemp and 
					 * zero extends the upper bits. */

    /* lazy F-loop */
    while (cmpval != 0xffff) { /* not all 8-bit integers == 0 */
      vE = _mm_load_si128 (vEp + j);
      vH = _mm_max_epu8 (vH, vF);

      /* save vH values */
      _mm_store_si128 (vHSp + j, vH);
      
      /*  update vE incase the new vH value would change it */
      vH = _mm_subs_epu8 (vH, vGapI);
      vE = _mm_max_epu8 (vE, vH);
      _mm_store_si128 (vEp + j, vE);

      /* update vF value */
      vF = _mm_subs_epu8 (vF, vGapE);

      j++;
      if (j >= segsiz) {
	j = 0;
	vF = _mm_slli_si128 (vF, 1);
      }

      vH = _mm_load_si128 (vHSp + j);

      vTmp = _mm_subs_epu8 (vH, vGapI);
      vTmp = _mm_subs_epu8 (vF, vTmp);
      vTmp = _mm_cmpeq_epi8 (vTmp, vZero);
      cmpval  = _mm_movemask_epi8 (vTmp);
    }
#ifdef alignment_matrix_debug
    printf("  [%i,%i]: ", i, 0);
    for (q=0; q<qlen; q++)
      printf("%4i", q);
    unprof_symbol = decoderp[(UCHAR) usqp[i]];
    printf("\n%c [%i,%i]: |", unprof_symbol, i, 0);
    for (q=0; q<qlen; q++)
      printf(" %c |", decoderp[(UCHAR) psqp[q]]);
    printf("\nHp[%i,%i]: |", i, 0);
    printfStripedByteVector((unsigned char *) vHSp, qlen, segsiz);
    printf("\n");
#endif
  }

  /* find largest score in the maxscorv vector */
  vTmp = _mm_srli_si128 (vMax, 8);
  vMax = _mm_max_epu8 (vMax, vTmp);
  vTmp = _mm_srli_si128 (vMax, 4);
  vMax = _mm_max_epu8 (vMax, vTmp);
  vTmp = _mm_srli_si128 (vMax, 2);
  vMax = _mm_max_epu8 (vMax, vTmp);
  vTmp = _mm_srli_si128 (vMax, 1);
  vMax = _mm_max_epu8 (vMax, vTmp);

  /* store in temporary variable */
  score = _mm_extract_epi16 (vMax, 0);
  score &= BYTEMASK;

  /*  check if we might have overflowed */
  if (score + bias >= UCHAR_MAX) {
    errcode = ERRCODE_SWATEXCEED;
  } else {
    errcode = ERRCODE_SUCCESS;
    *maxscor = (UCHAR) score; /* return largest score */
  }
  return errcode;
}


/******************************************************************************
 ********************** Public SIMD Alignment Methods *************************
 ******************************************************************************/

int swAlignStripedSSE2(ALIDPMSCOR_t *maxscor, 
		       AliBuffer *abp,
		       const ScoreProfile *profp, 
#ifdef alignment_matrix_debug 
		       const SeqCodec *codecp,
		       const char *profiled_seqp,
#endif
		       const char *unprofiled_seqp,
		       int unprofiled_seqlen)
{
  int errcode;
  unsigned char bytscor;
  unsigned short shortscor;

  *maxscor = 0;
  errcode = alignSmiWatByteStriped(&bytscor, 
				   abp, 
				   profp,
#ifdef alignment_matrix_debug
				   profiled_seqp,
				   codecp,
#endif 
				   unprofiled_seqp,
				   unprofiled_seqlen);
  if (errcode == ERRCODE_SWATEXCEED) {
    errcode = alignSmiWatShortStriped(&shortscor, 
				      abp, 
				      profp,
#ifdef alignment_matrix_debug
				      profiled_seqp,
				      codecp,
#endif 
				      unprofiled_seqp,
				      unprofiled_seqlen);
    if (!errcode) 
      *maxscor = (int) shortscor;
  } else if (!errcode) {
    *maxscor = (int) bytscor;
  }

  return errcode;
}

#ifdef alignment_timing
int swAlingStripedDirect(AliRsltSet *rssp, AliBuffer *bufp,
			 const ScoreProfile *profp,
			 const SeqCodec *codecp,
			 const char *psqp,
			 const char *usqp, int us_len,
			 ALIDPMSCOR_t *maxscor)
{
  int errcode;
  unsigned short  fastscor;
  unsigned char fastbytscor;

  *maxscor = 0;

  if (us_len < 1) 
    return ERRCODE_ASSERT;

  errcode = alignSmiWatByteStriped(&fastbytscor, bufp, profp,
#ifdef alignment_matrix_debug
				   psqp,
				   codecp,
#endif 
				   usqp, us_len);
  if (errcode == ERRCODE_SWATEXCEED) {
    errcode = alignSmiWatShortStriped(&fastscor, bufp, profp,
#ifdef alignment_matrix_debug
				      psqp,
				      codecp,
#endif
				      usqp, us_len);
    *maxscor = (ALIDPMSCOR_t) fastbytscor;
  } else {
    *maxscor = (ALIDPMSCOR_t) fastscor;
  }
  return errcode;
}
#endif //#ifdef alignment_timing
