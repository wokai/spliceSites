/**************************************************************************************************
 **************************************************************************************************
 *
 * Project	:	spliceSites
 * Created	:	01.07.2013
 * Author	:	W. Kaisers
 *
 * Content	:	Managing splice sites objects
 *
 * Version	:	0.4.2
 *
 * Changelog	:
 * 03.Jul.13	:	Added tryp_seq, trunc_pos and silic_tryp_pos functions.
 * 06.Jun.12	:
 * 13.Jun.12	:
 * 17.Okt.12	:
 * 07.May.13    :
 * 08.May.13    :
 *
 **************************************************************************************************
 **************************************************************************************************/


#ifndef SPLICESITES_H_
#define SPLICESITES_H_

#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <stdbool.h>		// bool (true,false)
#include <math.h>  			// log2 (-lm)
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h> // DllInfo
#include <Rdefines.h>
#include <R_ext/PrtUtil.h>
#include "alpha_idx.h"


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Declarations for trunc_pos and silic_tryp
//
///////////////////////////////////////////////////////////////////////////////////////////////////

SEXP trunc_pos(SEXP qid, SEXP qpos, SEXP qseq, SEXP pTrunc);
SEXP silic_tryp_pos(SEXP qid, SEXP qpos, SEXP qseq);

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function declarations:
// 		scoreconsensus5, maxent_score5
//		scoreconsensus3, maxent_score3
//
///////////////////////////////////////////////////////////////////////////////////////////////////


static inline double scoreconsensus5(unsigned char c3, unsigned char c4);
static inline double scoreconsensus3(unsigned char c18,unsigned char c19);

SEXP maxent_score5(SEXP pSeq, SEXP pPos, SEXP pMe2x5);
SEXP maxent_score3(SEXP pSeq, SEXP pPos, SEXP pMeList);

SEXP maxent_seq_score5(SEXP pSeq, SEXP pFrame,SEXP pMe2x5);
SEXP maxent_seq_score3(SEXP pSeq, SEXP pFrame,SEXP pMeList);

SEXP maxent_score2strand(SEXP pPscore, SEXP pMscore);
SEXP maxent_combine_strand(SEXP lStrand, SEXP rStrand);

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// 		hbond_score
//
///////////////////////////////////////////////////////////////////////////////////////////////////

SEXP hbond_score(SEXP pSeq,SEXP pHb);
SEXP create_dna_n_mers(SEXP pLen);

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Declarations for alt_group
//
///////////////////////////////////////////////////////////////////////////////////////////////////

SEXP alt_group(SEXP pId, SEXP pSeq, SEXP pGroup, SEXP pAlt);

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Declarations for alph_index
//
///////////////////////////////////////////////////////////////////////////////////////////////////


SEXP get_alph_index(SEXP pIdxLen, SEXP pAlpha);

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Declarations for R_registerRoutines
//
///////////////////////////////////////////////////////////////////////////////////////////////////

void R_init_spliceSites(DllInfo *info);

#endif /* SPLICESITES_H_ */

