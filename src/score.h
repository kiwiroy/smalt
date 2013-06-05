/** Score matrix and sequence profile for dynamic programming */

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

#ifdef __cplusplus
extern "C"
{
#endif

#ifndef score_h
#define score_h

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include "sequence.h"

  /****************************************************************************
   ****************************** Constants ***********************************
   ****************************************************************************/
  enum SCORE_PENALTIES {
    SCORPNLTY_MAXVAL = 127, /**< maximum of absolute value for an affine penalty score
			     * (e.g. for opening a gap or for a mismatch) */ 
  };

  enum SCORE_PENALTY_TYPES {
    SCORPNLTYP_MATCH = 0,
    SCORPNLTYP_MISMATCH = 1,
    SCORPNLTYP_GAPOPEN = 2,
    SCORPNLTYP_GAPEXT = 3,
    SCORPNLTYP_NUM = 4
  };

  enum SCORE_SIMD_CONST {
    SCORSIMD_NBYTES = 16,             /**< number of bytes in 128 bit SIMD register */
    SCORSIMD_NSHORTS = 8,             /**< number of short integers in 128 bit SIMD register */
    SCORSIMD_ROUNDMASK = 0x0f,       /**< Mask used for pointer alignment to 16-byte boundaries */
  };

  enum SCORE_PROFILE_MODES {
    SCORPROF_SCALAR = 0x01,     /**< for Smith-Waterman w/o SSE2 */
    SCORPROF_STRIPED_8 = 0x02,  /**< striped for 8-bit scores using SIMD */
    SCORPROF_STRIPED_16 = 0x04, /**< striped for 16-bit scores using SIMD */
  };
  
  /****************************************************************************
   ***************************** Opaque Types *********************************
   ****************************************************************************/
  typedef struct ScorePenalties_ ScorePenalties;
  /**< Alignment penalties (mismatch, gap opening, extension for dynamic programming */

  typedef struct _ScoreMatrix ScoreMatrix;
  /**< Matrix of alignment scores for dynamic programming */

  typedef struct _ScoreProfile ScoreProfile;
   /**< Sequence profile for dynamic programming */
  
  
   /******************************************************************************
   *********************** Methods of Type ScorePenalties ***********************
   ******************************************************************************/

  ScorePenalties *scorePenaltiesCreate(void);
  /**< Constructor 
   */
  void scorePenaltiesDelete(ScorePenalties *p);
  /**< Destructor 
   */
  int scoreSetPenalty(ScorePenalties *p, short typ, short penalty);
  /**< Set a penalty score for one of SCORE_PENALTY_TYPES.
   * \param p Set of penalties.
   * \param typ Penalty type (one of SCORE_PENALTY_TYPES).
   * \param penalty Value -SCORPLTY_MAXVAL <= penalty <= +SCORPLTY_MAXVAL.
   */

  short scoreGetPenalties(const ScorePenalties *p, 
			  short *mismatch, short *gapinit, short *gapext);
  /**< Accessor returning the alignment score for a match.
   * \param mismatch Returns the score for a mismatch (can be NULL)
   * \param gapinit Returns the score for opening a gap (can be NULL)
   * \param gapext Returns the score for extending a gap (can be NULL)
   */
 
  /******************************************************************************
   ************************* Methods of Type ScoreMatrix ************************
   ******************************************************************************/

  ScoreMatrix *scoreCreateMatrix(const SeqCodec *scp, const ScorePenalties *penp);
  /**< Constructor 
   */
  void scoreDeleteMatrix(ScoreMatrix *amp);
  /**< Destructor 
   */

  short scoreMatrixGetAlphabetSize(const ScoreMatrix *smp);
  /**< Accessor of the number of letters in the alphabet
   */

  double scoreMatrixCalcLambda(const ScoreMatrix *smp);
  /**< Calculate Poisson parameter \f$\lambda\f$.
   * as root of
   * \f$\sum_{ab} f_a f_b e^{\lambda s(a,b)} = 1\f$
   * this is for a,b in {1,..,4} and f_i = 1/4 for i={1,..4}
   */

  signed char scoreMatrixGetMinSubstScore(const ScoreMatrix *smp);
  /**< Returns the minimum substitution score.
   */

  signed char scoreMatrixGetAvgSubstScores(signed char *matchscor, const ScoreMatrix *smp);
  /**< Return the average substituiton score.
   * \param matchscor Returns average match score (can be NULL).
   * \param smp Score matrix.
   */

  short scoreGetDefaults(short *mismatch, short *gapinit, short *gapext);
  /**< Accessor returning the default Smith-Waterman score for a match.
   * \param mismatch Returns the default score for a mismatch (can be NULL)
   * \param gapinit Returns the default score for opening a gap (can be NULL)
   * \param gapext Returns the default score for extending a gap (can be NULL)
   */
  
  int scorePrintMatrix(const ScoreMatrix *amp, const SeqCodec *scp);
  /**< Print alignment matrix on standard output 
   */
  
  /******************************************************************************
   *********************** Methods of Type ScoreProfile *************************
   ******************************************************************************/
  
  ScoreProfile *scoreCreateProfile(int blocksize, const SeqCodec *codep, unsigned char mod);
  /**< Constructor.
   * \param blocksize Blocksize for memory allocation 
   * \param codep Sequence en-/decoder.
   * \param mod Striping mode, a combination of SCORE_PROFILE_MODES bit flags
   */

  void scoreDeleteProfile(ScoreProfile *app);
  /**< Destructor */

  int scoreMakeProfileFromSequence(ScoreProfile *app, const SeqFastq *sqp, 
				   const ScoreMatrix *amp, const SeqCodec *codep);
  /**< Make a profile for a sequence. 
   * \param app Pointer to profile.
   * \param sqp Pointer to sequence.
   * \param amp Alignment matrix.
   * \param codep Nucleotide code used.
   */

  signed char * const *scoreGetProfile(short *alphabetsiz, SEQLEN_t *seqlen, 
				       signed char *gap_init, signed char *gap_ext,  
				       const ScoreProfile *spp);

  /**< Accessor returning the sequence profile matrix score[seqlen][aphabetsiz].
   * \param alphabetsiz Returns size of the sequence alphabet (fast running index).
   * \param seqlen Returns sequence length (slow running index).
   * \param gap_init Gap opening score as positive number
   *        (sign inverted penalty).
   * \param gap_ext Gap extension score as a positive number
   *        (sign inverted penalty).
   * \param spp Sequence profile.
   */

  short scoreProfileGetAvgPenalties(short *msmatch_avg,
				    short *gap_init, short *gap_ext,
				    const ScoreProfile *spp);
  /**< Return average match score,
   * \param match_avg Average mismatch penalty.
   * \param gap_init Penalty score for opening a gap.
   * \param gap_ext Penalty score for extending a gap.
   * \param spp Score profile.
   */

#ifdef HAVE_EMMINTRIN_H
  const void *scoreGetStripedProfile(short *alphabetsiz, SEQLEN_t *seqlen, 
				     unsigned short *gap_init, unsigned short *gap_ext,
				     unsigned short *bias, int *segsiz, char mod,
				     const ScoreProfile *spp);
  /**< Accessor returning the striped sequence profile.
   * \param alphabetsiz Returns size of the sequence alphabet (fast running index).
   * \param seqlen Returns sequence length (slow running index).
   * \param gap_init Gap opening score as positive number
   *        (sign inverted penalty).
   * \param gap_ext Gap extension score as a positive number
   *        (sign inverted penalty).
   * \param bias Minimum score of the profile.
   * \param segsiz Size of a striped 'segment'.
   * \param mod Striping mode: byte (SCORPROF_STRIPED_8) or short words (SCORPROF_STRIPED_16)
   */
 
#endif

#endif
#ifdef __cplusplus
}
#endif
