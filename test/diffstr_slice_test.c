/** Testing handling of alignment strings
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

#include "elib.h"
#include "diffstr.h"

enum {
  DIFFSTR_BLKSZ = 16,
};

static const char ALISTR[] = "32M16I8D8S";

int sliceDiffStr(FILE *fp,  DiffStr *dfsbfp, 
		 int start_unprof_target, int end_unprof_target,
		 const DIFFSTR_T *diffstrp) 
{
  int errcode;
  int start_unprof, end_unprof;
  int start_prof, end_prof;

  errcode = diffStrSegment(dfsbfp, diffstrp, 
			   start_unprof_target, end_unprof_target,
			   &start_unprof, &end_unprof,
			   &start_prof, &end_prof);
  if (errcode)
    return errcode;

  fprintf(fp, "start_unprof_target = %i\n", start_unprof_target);
  fprintf(fp, "end_unprof_target = %i\n", end_unprof_target);
  fprintf(fp, "start_unprof = %i\n", start_unprof);
  fprintf(fp, "end_unprof = %i\n", end_unprof);
  fprintf(fp, "start_prof = %i\n", start_prof);
  fprintf(fp, "end_prof = %i\n", end_prof);

  if ((errcode = diffStrPrintf(fp, dfsbfp->dstrp, dfsbfp->len, DIFFSTRFORM_RAW, 0, 0, 0)))
    return errcode;
  
  errcode = diffStrPrintf(fp, dfsbfp->dstrp, dfsbfp->len, DIFFSTRFORM_CIGNORM, 0, 0, 0);
  fprintf(fp, "\n");

  return errcode;
}

int main() 
{
  int errcode = ERRCODE_SUCCESS;
  DiffStr *dfsp, *dfsbfp;
  ErrMsg *errmsgp = 0;

  ERRMSG_CREATE(errmsgp);

  dfsp = diffStrCreate(DIFFSTR_BLKSZ);
  if (dfsp == NULL)
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);

  dfsbfp = diffStrCreate(DIFFSTR_BLKSZ);
  if (dfsbfp == NULL)
    ERRMSGNO(errmsgp, ERRCODE_NOMEM);
  
  if ((errcode = diffStrParsePlain(dfsp, ALISTR)))
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = diffStrPrintf(stdout, dfsp->dstrp, dfsp->len, 
			       DIFFSTRFORM_RAW, 0, 0, 0)))
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = diffStrPrintf(stdout, dfsp->dstrp, dfsp->len,
			       DIFFSTRFORM_CIGNORM, 0, 0, 0)))
    ERRMSGNO(errmsgp, errcode);
  
  if ((errcode = sliceDiffStr(stdout, dfsbfp, 3, 20, dfsp->dstrp)))
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = sliceDiffStr(stdout, dfsbfp, 10, 49, dfsp->dstrp)))
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = sliceDiffStr(stdout, dfsbfp, 49, 57, dfsp->dstrp)))
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = sliceDiffStr(stdout, dfsbfp, 48, 58, dfsp->dstrp)))
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = sliceDiffStr(stdout, dfsbfp, 48, 57, dfsp->dstrp)))
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = sliceDiffStr(stdout, dfsbfp, 49, 58, dfsp->dstrp)))
    ERRMSGNO(errmsgp, errcode);
  
  if ((errcode = sliceDiffStr(stdout, dfsbfp, 59, 66, dfsp->dstrp)))
    ERRMSGNO(errmsgp, errcode);
  
  if ((errcode = sliceDiffStr(stdout, dfsbfp, 58, 65, dfsp->dstrp)))
    ERRMSGNO(errmsgp, errcode);

  if ((errcode = sliceDiffStr(stdout, dfsbfp, 57, 65, dfsp->dstrp)))
    ERRMSGNO(errmsgp, errcode);
  
  if ((errcode = sliceDiffStr(stdout, dfsbfp, 57, 68, dfsp->dstrp)))
    ERRMSGNO(errmsgp, errcode);

  /* the following should fail */
  errcode = sliceDiffStr(stdout, dfsbfp, 57, 57, dfsp->dstrp);
  if (!errcode)
    ERRMSGNO(errmsgp, ERRCODE_DIFFSTR);

  if ((errcode = sliceDiffStr(stdout, dfsbfp, 58, 58, dfsp->dstrp)))
    ERRMSGNO(errmsgp, errcode);
  
  if ((errcode = sliceDiffStr(stdout, dfsbfp, 0, 0, dfsp->dstrp)))
    ERRMSGNO(errmsgp, errcode);
  
  diffStrDelete(dfsbfp);
  diffStrDelete(dfsp);

  ERRMSG_END(errmsgp);

  exit(errcode);
}
