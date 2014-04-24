/** Test generic interface to sequence input formats
 */

/*****************************************************************************
 *****************************************************************************
 *                                                                           *
 *  Copyright (C) 2010 - 2014 Genome Research Ltd.                           *
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

#include <stdio.h>
#include <stdlib.h>

#include "infmt.h"

int main(int argc, char *argv[]) 
{
  int errcode = ERRCODE_SUCCESS;
  char *infilnam;
#ifdef HAVE_BAMBAMC
  char *tmpdir = NULL;
#endif
  char *infilnam_mates = NULL;
  unsigned char isPair = 0;
  INFMT_t fmt = INFMT_UNKNOWN;
  size_t sctr = 0, pctr = 0;
  SeqFastq *readAp, *readBp;
  SeqIO *sfop;
  InFmtReader *ifrp;
  ErrMsg *errmsg = 0;

  ERRMSG_CREATE(errmsg);

  if (argc < 3 || argc > 5) {
    printf("usage: %s <input file> <%s>"	\
	   "[<tmpdir> [<input file B>]]\n",
#ifdef HAVE_BAMBAMC
	   "f|b|u (FASTA|BAM|UNKNOWN format)"
#else
	   "f|u (FASTA|UNKNOWN format)"
#endif
	   , argv[0]);
    
    exit(EXIT_FAILURE);
  }

  infilnam = argv[1];
  if (argv[2][0] == 'f') 
    fmt = INFMT_FASTQ;
#ifdef HAVE_BAMBAMC
  else if (argv[2][0] == 'b')
    fmt = INFMT_BAM;

  if (argc > 3) 
    tmpdir = argv[3];
#endif
  if (argc > 4)
    infilnam_mates = argv[4];

  if (NULL == (readAp = seqFastqCreate(0, SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsg, ERRCODE_NOMEM);

  if (NULL == (readBp = seqFastqCreate(0, SEQTYP_UNKNOWN)))
    ERRMSGNO(errmsg, ERRCODE_NOMEM);
  
  ifrp = infmtCreateReader(&errcode, infilnam, infilnam_mates,
#ifdef HAVE_BAMBAMC
			   tmpdir, 
#endif
			   fmt);
  if (NULL == ifrp)
    ERRMSGNO(errmsg, errcode);

  sfop = seqIOopen(&errcode, "-", SEQIO_WRITE_FASTQ, 0);
  if (NULL == sfop)
    ERRMSGNO(errmsg, errcode);

  while (!(errcode = infmtRead(ifrp, readAp, readBp, &isPair))) {

    if ((errcode = seqFastqWrite(sfop, readAp, 0)))
      ERRMSGNO(errmsg, errcode);

    if (isPair) {
      pctr++;
      if ((errcode = seqFastqWrite(sfop, readBp, 0)))
	ERRMSGNO(errmsg, errcode);
    } else {
      sctr++;
    }
  }
  if (errcode != ERRCODE_EOF)
    ERRMSGNO(errmsg, errcode);

  printf("%llu read pairs.\n", (unsigned long long) pctr);
  printf("%llu single reads.\n", (unsigned long long) sctr);
  if (seqIOstatus(sfop)) 
    ERRMSGNO(errmsg, seqIOstatus(sfop));

  seqIOclose(sfop);
  infmtDeleteReader(ifrp);
  seqFastqDelete(readBp);
  seqFastqDelete(readAp);

  ERRMSG_END(errmsg);

  exit(EXIT_SUCCESS);
}
