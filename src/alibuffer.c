/**< Buffers for pairwise alignments */

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

#include <stdlib.h>
#include <limits.h>
#include <config.h>
#ifdef HAVE_EMMINTRIN_H
#include <emmintrin.h>
#endif

#include "elib.h"
#include "score.h"
#include "alibuffer.h"
#include "alibuffer_struct.h"

enum {
  SLACK_BYTES = 32,
  ALIBUFFER_DEFAULT_BLKSZ = 256,
};

/****************************************************************************
 ******************************* Macros *************************************
 ****************************************************************************/

#define ALIGN_16BYTE(p) (((size_t) (p) + SCORSIMD_ROUNDMASK) & (~SCORSIMD_ROUNDMASK));
/**< Align memory to 16 byte boundary */

/******************************************************************************
 ********************** Public Methods of Type AliBuffer **********************
 ******************************************************************************/

AliBuffer *aliBufferCreate(int blocksiz)
{
  AliBuffer *p;
  EMALLOCP0(p);
  if (!p) return 0;

  if (blocksiz < 1) blocksiz = ALIBUFFER_DEFAULT_BLKSZ;

  ECALLOCP(blocksiz, p->datap);
  if (!p->datap) {
    aliBufferDelete(p);
    return 0;
  }

  p->allocsiz = (size_t ) blocksiz;
  p->blocksiz = blocksiz;
  p->qlen_max = 0; /* not yet initialised (pointers not set) */
  return p;
}

void aliBufferDelete(AliBuffer *p)
{
  if (p)
    free(p->datap);
  free(p);
}

int aliBufferInit(AliBuffer *p, unsigned int qlen)
{
  unsigned int newlen;
  size_t siz, newsiz;
#ifdef HAVE_EMMINTRIN_H
  int qlenSHRT;
#endif 
  if (qlen <= p->qlen_max)
    return ERRCODE_SUCCESS;

  if (qlen > INT_MAX)
    return ERRCODE_SEQLEN;

  newlen = (qlen + p->blocksiz - 1)/p->blocksiz;
  newlen *= p->blocksiz;
  if (newlen > INT_MAX)
    newlen = INT_MAX;
  siz = (newlen+1)*2*sizeof(int);

#ifdef HAVE_EMMINTRIN_H
  qlenSHRT = (newlen + SCORSIMD_NSHORTS - 1)/SCORSIMD_NSHORTS;
  newsiz = SLACK_BYTES + 3 * qlenSHRT * sizeof(__m128i);
  if (siz < newsiz)
    siz = newsiz;
#endif  

  if (siz > p->allocsiz) {
    void *hp;
    newsiz = (siz + p->blocksiz - 1)/p->blocksiz;
    newsiz *= p->blocksiz;

    hp = EREALLOCP(p->datap, newsiz);
    if (!hp)
      return ERRCODE_NOMEM;
    p->datap = hp;
    p->allocsiz = newsiz;
  }
  p->baseHp = (ALIDPMSCOR_t *) ALIGN_16BYTE(p->datap);
  p->baseEp = p->baseHp + newlen + 1;
  p->baseEp = (ALIDPMSCOR_t *) ALIGN_16BYTE(p->baseEp);
  

#ifdef HAVE_EMMINTRIN_H
  /* align to 16 byte boundary for 128 bit SIMD register */
  p->H1v = (__m128i *) ALIGN_16BYTE(p->datap);
  p->H2v = p->H1v + qlenSHRT;
  p->Ev  = p->H2v + qlenSHRT;
#endif

  p->qlen_max = (int) newlen;

  return ERRCODE_SUCCESS;
}
