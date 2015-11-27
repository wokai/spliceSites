/**************************************************************************************************
 **************************************************************************************************
 *
 * Project	:	spliceSites
 * Created	:	01.07.2013
 * Author	:	W. Kaisers
 *
 * Content	:	Managing splice sites objects
 *
 * Version	:	0.5.0
 *
 * Changelog	:
 * 03.Jul.13	:	Added tryp_seq, trunc_pos and silic_tryp_pos functions.
 * 09.Jul.13	:	Added maxEntScore functions: scoreconsensus5/3 and maxent_score5/3
 * 13.Jun.12	:
 * 17.Okt.12	:
 * 07.May.13    :
 * 08.May.13    :
 *
 **************************************************************************************************
 **************************************************************************************************/


#ifndef SPLICESITES_C_
#define SPLICESITES_C_
#include "spliceSites.h"

static const int buf_size=2048; // buffer size for printing ints into chars

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Definitions for trunc_pos and silic_tryp
//
///////////////////////////////////////////////////////////////////////////////////////////////////

// Is 1 for ASCII of RKrk
static const unsigned char rk[256] = {
		// 64 - block							// start index
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 0
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 16
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 32
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 48


		0,0,0,0,	0,0,0,0,	0,0,0,1,	0,0,0,0,	// 64
		0,0,1,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 80
		0,0,0,0,	0,0,0,0,	0,0,0,1,	0,0,0,0,	// 96
		0,0,1,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 112

		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 128
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 144
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 160
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 176

		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 192
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 208
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 224
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0		// 240
};

// Used for tryp_seq
// Is one for ASCII of Pp and for 0
static const unsigned char pp0[256] = {
		// 64 - block									// start index
		1,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 0
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 16
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 32
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 48


		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 64
		1,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 80
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 96
		1,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 112

		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 128
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 144
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 160
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 176

		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 192
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 208
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0,	// 224
		0,0,0,0,	0,0,0,0,	0,0,0,0,	0,0,0,0		// 240
};


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Static variable declararions for maxEntScores (scoreconsensus and maxent_score functions)
//
///////////////////////////////////////////////////////////////////////////////////////////////////

// Translation tables: DNA-Nuc -> integer index
# define a_val 0
# define c_val 1
# define g_val 2
# define t_val 3
# define zvl   4

// ASCII:	A 65, C 67, G 71, T 84, a 97, c 99, g 1zvl3, t 116
// Translation:	A->0, C->1, G->2, T->3
static const unsigned char ACGT[256] = {
		// 64 - block																		// start index
		zvl,zvl,zvl,zvl,		zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	// 0
		zvl,zvl,zvl,zvl,		zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	// 16
		zvl,zvl,zvl,zvl,		zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	// 32
		zvl,zvl,zvl,zvl,		zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	// 48

		zvl,a_val,zvl,c_val,	zvl,zvl,zvl,g_val,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	// 64
		zvl,zvl,zvl,zvl,		t_val,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	// 80
		zvl,a_val,zvl,c_val,	zvl,zvl,zvl,g_val,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	// 96
		zvl,zvl,zvl,zvl,		t_val,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	// 112

		zvl,zvl,zvl,zvl,		zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	// 128
		zvl,zvl,zvl,zvl,		zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	// 144
		zvl,zvl,zvl,zvl,		zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	// 160
		zvl,zvl,zvl,zvl,		zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	// 176

		zvl,zvl,zvl,zvl,		zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	// 192
		zvl,zvl,zvl,zvl,		zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	// 208
		zvl,zvl,zvl,zvl,		zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	// 224
		zvl,zvl,zvl,zvl,		zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl,	zvl,zvl,zvl,zvl		// 240
};

// Character encoding for create_dna_n_mers
static const unsigned char rev_ACGT[4] = { 65, 67, 71, 84 };

// Global score arrays:
//				A	C	G	T
static const double bgd[4]    = {	0.27  ,	0.23  ,	0.23  ,	0.27};
static const double cons51[4] = {	0.004 ,	0.0032,	0.9896,	0.0032};
static const double cons52[4] = {	0.0034,	0.0039,	0.0042,	0.9884};
static const double cons31[4] = {	0.9903,	0.0032,	0.0034,	0.0030};
static const double cons32[4] = {	0.0027,	0.0037,	0.9905,	0.0030};

///////////////////////////////////////////////////////////////////////////////////////////////////
// maxent_score5 specific definitions:
// const unsigned four5[7] = {4096,1024,256,64,16,4,1};
static const unsigned l5seq=9;		// Sequence length for score5 (extracted from pSeq)
static const unsigned n5ExonNucs=3;	// Number of exonic nucleotides which are read upstream of pos


///////////////////////////////////////////////////////////////////////////////////////////////////
// maxent_score3 specific definitions
static const unsigned l3seq=23;	// Sequence length for score3 (extracted from pSeq)
static const unsigned nsc3=9;		// Number of primary Scores in score3 function
static const unsigned n3ExonNucs=20;	// Number of exonic nucleotides which are read upstream of pos


// score3 square schema: 7 columns, 9 rows
// Used for calculation, done inside of the score3.pl "maxentscore" function (i.e. the substrings).
// The array includes a correction for usage of "getrest" function in score3.pl:
// Nucleotides 18 and 19 (0-based) are excluded (they are used in the scoreconsensus3 function.
// This allows usage of command "substr($seq,14,7))" in score3.pl, line 123.
// Because sc3c addresses each nucleotide explicitly, it's possible to
// include the "getrest"-induced shift here.

static const unsigned sc3c[63] = {
		 0,  1,  2,  3,  4,  5,  6,
		 7,  8,  9, 10, 11, 12, 13,
		14, 15, 16, 17, 20, 21, 22,
		 4,  5,  6,  7,  8,  9, 10,
		11, 12, 13, 14, 15, 16, 17,
		 4,  5,  6,  0,  0,  0,  0,
		 7,  8,  9, 10,  0,  0,  0,
		11, 12, 13,  0,  0,  0,  0,
		14, 15, 16, 17,  0,  0,  0
};

static const unsigned sc3l[9]  = { 7, 7, 7, 7, 7, 3, 4, 3, 4 }; // Number of entries in each sc3c line
static const unsigned sc3s[9]  = { 0, 7,14,21,28,35,42,49,56 }; // Index  of first entry in each sc3c line
static const unsigned four3[]  = { 1, 4,16,64,256,1024,4096,16384 };
///////////////////////////////////////////////////////////////////////////////////////////////////

static const unsigned strand_plus=1;
static const unsigned strand_minus=2;
static const unsigned strand_star=3;


///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function definitions: tryp_seq and silic_tryp
//
///////////////////////////////////////////////////////////////////////////////////////////////////

/*
 *  In silico trypsinisation regarding position in sequence string
 *
 *  Identification of Splice-site: [RK](?!P) ('\0' after [RK] is
 *  also accepted: Implements the "Keil rule"
 *
 *  The function cuts the sequence behind [RK] and before (?!P).
 *  Several cases are considered, so that the cut always contains
 *  the position marked nucleotide
 *
 *	The position coordinate is updated (when neccessary)
 *	seqlen is only used for returning value
 *
 *	The rationale of the function is as follows: The position
 *	represents a splice site (0-based). The (AA-) sequence
 *	is trypsinated and the fragment which contains the splice-site
 *	should be returned. When no truncation site is found, a copy
 *	of the input string is returned.
 *
 *
 */

char * tryp_seq(const char* seq, unsigned *pos, unsigned *seqlen)
{
	int i=0;
	int start=-1; // start <0: value is still unset
	char *res=0;

	while(seq[i])
	{
		// regex: [RK](?!P)
		// (terminal '\0' allowed)
		// This means: trypsinization site found.
		if(rk[(unsigned)seq[i]] & (!pp0[(unsigned)seq[i+1]]))
		{
			if(start<0)
			{
				if(i>=*pos)
				{
					// xxxxPxTx
					// Truncate right
					*seqlen=i+1;
					res= R_alloc(*seqlen+1,sizeof(char));
					strncpy(res,seq,*seqlen);
					res[*seqlen]='\0';
					return res;
				}
				// Set start (>=0)
				start=i+1;
			}
			else
			{
				// Now: start>=0 (i.e. has already been set)
				if(i>=*pos)
				{
					// xxTxPxTx
					// TxTxPxxx
					// Truncate left and right
					*seqlen=i-start+1;
					*pos-=start;
					res= R_alloc(*seqlen+1,sizeof(char));
					strncpy(res,seq+start,*seqlen);
					res[*seqlen]='\0';
					return res;
				}
				else
					start=i+1; // Shift start pos to right
			}
		}
		++i;
	}

	if(start<0)
	{
		// xxxxPxxx
		// Just copy (no trunc)
		*seqlen=i;
		res= R_alloc(*seqlen+1,sizeof(char));
		strcpy(res,seq);
		res[*seqlen]='\0';
		return res;
	}

	// xxTxPxxx
	// Truncate left
	*seqlen=i-start;
	*pos-=start;
	res= R_alloc(*seqlen+1,sizeof(char));
	strcpy(res,seq+start);
	res[*seqlen]='\0';
	return res;
}


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
//
// trunc_pos
//
// Traverses (amino-acid) sequence string (qseq) from left to right
// and truncates at first occurence of character encoded by pTrunc
// depending on position relative to qpos:
//
// When pTrunc is found left of qpos, an empty string is returned. When pTrunc is found right
// of qpos, a copy of the strung up to the identified position is returned.
//
// The function returns a data.frame where nrows=length(qid)
//
// The rationale of the function is as follows: qpos marks the position of the splice site inside
// of an amino acid sequence. pTrunc=42 which encodes '*', the stop-codon. When the stop-codon is
// on the left side of the splice site, the splice site will not be be translated into protein.
// So it can be omitted.
// Otherwise the amino acid sequence is truncated at the stop-codon site.
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

SEXP trunc_pos(SEXP qid, SEXP qpos, SEXP qseq, SEXP pTrunc)
{
	if(TYPEOF(qid)!=INTSXP)
		error("[trunc_pos] qid is not INT!\n");
	if(TYPEOF(qpos)!=INTSXP)
		error("[trunc_pos] qpos is not INT!\n");
	if(TYPEOF(qseq)!=STRSXP)
		error("[trunc_pos] qseq is not STRING!\n");
	if(TYPEOF(pTrunc)!=INTSXP)
		error("[trunc_pos] pTrunc is not INT!\n");

	unsigned nRows=LENGTH(qid);
	if(nRows!=LENGTH(qpos))
		error("[trunc_pos] length(qid)!=length(qpos)!\n");
	if(nRows!=LENGTH(qseq))
		error("[trunc_pos] length(qid)!=length(qseq)!\n");

	int trunc=INTEGER(pTrunc)[0];
	unsigned i,j, done=0;
	char * ret=0;
	unsigned nProtected=0;
	//Rprintf("[trunc_pos] nRows=%u\ttrunc=%i\n",nRows,trunc);


	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Create data.frame
	unsigned nCols=4;

	// Column 0: id
	SEXP rid;
	PROTECT(rid=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 1: pos
	SEXP rpos;
	PROTECT(rpos=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 2: seq
	SEXP rseq;
	PROTECT(rseq=allocVector(STRSXP,nRows));
	++nProtected;

	// Column 3: strlen
	SEXP rstrlen;
	PROTECT(rstrlen=allocVector(INTSXP,nRows));
	++nProtected;


	for(i=0;i<nRows;++i)
	{
		//Rprintf("[trunc_pos] i: %u\tpos: %i\n",i,INTEGER(qpos)[i]);
		INTEGER(rid)[i]=INTEGER(qid)[i];
		const char *seq=CHAR(STRING_ELT(qseq,i));
		j=0;
		while(seq[j])
		{
			if(((unsigned)seq[j])==trunc)
			{
				// trunc left of pos: return empty string
				if(j<INTEGER(qpos)[i])
				{
					//Rprintf("[trunc_pos] trunc left of pos.\n");
					SET_STRING_ELT(rseq,i,mkChar(""));
					INTEGER(rstrlen)[i]=0;
					INTEGER(rpos)[i]=0;
					done=1;
					break;
				}

				// trunc on or right of pos: return truncated string
				ret=R_alloc(j+1,sizeof(char));
				strncpy(ret,seq,j);
				ret[j]='\0';
				SET_STRING_ELT(rseq,i,mkChar(ret));
				INTEGER(rstrlen)[i]=j;
				INTEGER(rpos)[i]=INTEGER(qpos)[i];
				done=1;
				break;
			}
			++j;

		}
		if(done==0)
		{
			//Rprintf("[trunc_pos] trunc not found.\n");
			// trunc not found: return copy of seq
			SET_STRING_ELT(rseq,i,mkChar(seq));
			INTEGER(rstrlen)[i]=j;
			INTEGER(rpos)[i]=INTEGER(qpos)[i];
		}
		else
			done=0;

	}


	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Create data.frame
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));

	SET_VECTOR_ELT(dflist,0,rid);
	SET_VECTOR_ELT(dflist,1,rpos);
	SET_VECTOR_ELT(dflist,2,rseq);
	SET_VECTOR_ELT(dflist,3,rstrlen);

	// Column names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names,0,mkChar("id"));
	SET_STRING_ELT(col_names,1,mkChar("pos"));
	SET_STRING_ELT(col_names,2,mkChar("seq"));
	SET_STRING_ELT(col_names,3,mkChar("lseq"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	// Row names
	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nRows));
    ++nProtected;
    char c[20];
    for(i=1;i<=nRows;++i)
    {
    	sprintf(c,"%i",i);
    	SET_STRING_ELT(row_names,i-1,mkChar(c));
    }
    setAttrib(dflist,R_RowNamesSymbol,row_names);

    // return
	UNPROTECT(nProtected);
	return dflist;
}


SEXP silic_tryp_pos(SEXP qid, SEXP qpos, SEXP qseq)
{
	if(TYPEOF(qid)!=INTSXP)
		error("[silic_tryp_pos] qid is not INT!\n");
	if(TYPEOF(qpos)!=INTSXP)
		error("[silic_tryp_pos] qpos is not INT!\n");
	if(TYPEOF(qseq)!=STRSXP)
		error("[silic_tryp_pos] qseq is not STRING!\n");

	unsigned nRows=LENGTH(qid);
	if(nRows!=LENGTH(qpos))
		error("[silic_tryp_pos] length(qid)!=length(qpos)!\n");
	if(nRows!=LENGTH(qseq))
		error("[silic_tryp_pos] length(qid)!=length(qseq)!\n");


	unsigned i;
	char * ret=0;
	unsigned nProtected=0;

	unsigned *seqlen=(unsigned*)R_alloc(1,sizeof(unsigned));
	unsigned *pos=(unsigned*)R_alloc(1,sizeof(unsigned));

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Create data.frame
	unsigned nCols=4;

	// Column 0: id
	SEXP rid;
	PROTECT(rid=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 1: pos
	SEXP rpos;
	PROTECT(rpos=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 2: seq
	SEXP rseq;
	PROTECT(rseq=allocVector(STRSXP,nRows));
	++nProtected;

	// Column 3: strlen
	SEXP rstrlen;
	PROTECT(rstrlen=allocVector(INTSXP,nRows));
	++nProtected;

	for(i=0;i<nRows;++i)
	{
		INTEGER(rid)[i]=INTEGER(qid)[i];
		*pos=INTEGER(qpos)[i];

		const char *seq=CHAR(STRING_ELT(qseq,i));
		ret=tryp_seq(seq,pos,seqlen);
		INTEGER(rpos)[i]=*pos;
		SET_STRING_ELT(rseq,i,mkChar(ret));
		INTEGER(rstrlen)[i]=*seqlen;
	}

	// + + + + + + + + + + + + + + + + + + + + + + + + + +
	// Create data.frame
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));

	SET_VECTOR_ELT(dflist,0,rid);
	SET_VECTOR_ELT(dflist,1,rpos);
	SET_VECTOR_ELT(dflist,2,rseq);
	SET_VECTOR_ELT(dflist,3,rstrlen);

	// Column names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names,0,mkChar("id"));
	SET_STRING_ELT(col_names,1,mkChar("pos"));
	SET_STRING_ELT(col_names,2,mkChar("seq"));
	SET_STRING_ELT(col_names,3,mkChar("lseq"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	// Row names
	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nRows));
    ++nProtected;
    char c[20];
    for(i=1;i<=nRows;++i)
    {
    	sprintf(c,"%i",i);
    	SET_STRING_ELT(row_names,i-1,mkChar(c));
    }
    setAttrib(dflist,R_RowNamesSymbol,row_names);

    // return
	UNPROTECT(nProtected);
	return dflist;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// Function definitions:
// 		scoreconsensus5, maxent_score5
//		scoreconsensus3, maxent_score3
//
// Original algorithm downloaded from:
//     http://genes.mit.edu/burgelab/maxent/download/fordownload.tar.gz
//
///////////////////////////////////////////////////////////////////////////////////////////////////


static inline double scoreconsensus5(unsigned char c3, unsigned char c4)
{
	return (cons51[c3]*cons52[c4])/(bgd[c3]*bgd[c4]);
}

static inline double scoreconsensus3(unsigned char c18,unsigned char c19)
{
	return ((cons31[c18]*cons32[c19])/(bgd[c18]*bgd[c19]));
}

///////////////////////////////////////////////////////////////////////////////////////////////////

SEXP maxent_seq_score5(SEXP pSeq, SEXP pFrame,SEXP pMe2x5)
{
	if(TYPEOF(pSeq)!=STRSXP)
		error("[maxent_seq_score5] pSeq must be STRING!");
	if(TYPEOF(pFrame)!=INTSXP)
		error("[maxent_seq_score5] pFrame must be INTEGER!");
	if(TYPEOF(pMe2x5)!=REALSXP)
		error("[maxent_seq_score5] pMe2x5 must be REAL!");

	if(STRING_ELT(pSeq,0)==NA_STRING)
	{
		Rprintf("[maxent_seq_score5] pSeq is NA!");
		return R_NilValue;
	}

    ////////////////////////////////////////////////////////
    //
    //    pos=1-based last exon nuc = 3, l5seq=9
    //    |
    //  ATG | GTC | ATCGAA
    //  123   456   789012
    //    |     |
    //    |      end=6
    //    start=3
    //
    ////////////////////////////////////////////////////////

	unsigned i,j,k,pos,npos;
	const char     *inseq=CHAR(STRING_ELT(pSeq,0));
	const unsigned  start=INTEGER(pFrame)[0];          // 1-based position of leftmost  score5 position
	const unsigned  end  =INTEGER(pFrame)[1];          // 1-based position of rightmost score5
	const unsigned  nchar=strlen(inseq);

	if(end<=start)
		error("[maxent_seq_score5] end must be > start!");
	if(start<n5ExonNucs)
		error("[maxent_seq_score5] start must be greater than %u!",n5ExonNucs);

	pos=end-n5ExonNucs+l5seq; // 1-based index of last examined nuc when pos=end is scored
	if(nchar<pos)
		error("[maxent_seq_score5] length of inseq must be greater than %u (or reduce end)!",pos);

	// number of positions for which score5 can be calculated
	npos=end-start+1;

	SEXP ans;
	PROTECT(ans=allocVector(REALSXP,npos));

	unsigned seq[l5seq];
	unsigned idx;

	for(pos=start,i=0; pos<=end; ++pos,++i)
	{
		// calculate score5
		// pos = 1-based position of last exonic nucleotide
		// Start position (pos-n5ExonNucs):
		// 			Includes last n5ExonNucs exonic nucleotides (when using 0-based array index)
		// Reads values for l5seq nuc's.
		for(j=0,k=pos-n5ExonNucs;j<l5seq;++j,++k)
		{
			seq[j]=ACGT[(unsigned)inseq[k]];
			if(seq[j]==zvl)
			{
				seq[0]=zvl;
				break;
			}
		}

		if(seq[0]!=zvl)
		{
			// numeric factors are statically used in this function (-> no variable)
			idx=seq[0]*4096+seq[1]*1024+seq[2]*256+seq[5]*64+seq[6]*16+seq[7]*4+seq[8];
			REAL(ans)[i]=log2(scoreconsensus5(seq[3],seq[4])*REAL(pMe2x5)[idx]);
		}
		else
			REAL(ans)[i]=NA_REAL;

	}
	UNPROTECT(1);
	return ans;
}


SEXP maxent_score5(SEXP pSeq, SEXP pPos, SEXP pMe2x5)
{
	if(TYPEOF(pSeq)!=STRSXP)
		error("[maxent_score5] pSeq must be STRING!");
	if(TYPEOF(pMe2x5)!=REALSXP)
		error("[maxent_score5] pMe2x5 must be REAL!");
	if(TYPEOF(pPos)!=INTSXP)
		error("[maxent_score5] pPos must be INT!");

	unsigned i,j,k,n=LENGTH(pSeq);
	if(n!=LENGTH(pPos))
		error("[maxent_score5] LENGTH(pSeq) must equal LENGTH(pPos)!");

	unsigned pos=INTEGER(pPos)[0];

	//Rprintf("[maxent_score5] nchar: %u\n",LENGTH(STRING_ELT(pSeq,0)));
	if(pos<n5ExonNucs)
		error("[score5.maxEntScore] pos must be >=%u! At least %u exon nucleotides needed!",n5ExonNucs,n5ExonNucs);
	unsigned seq[l5seq];
	unsigned idx;


	// For each sequence:
	SEXP ans;
	PROTECT(ans=allocVector(REALSXP,n));
	for(i=0;i<n;++i)
	{
		if(LENGTH(STRING_ELT(pSeq,i))<l5seq)
			error("[score5.maxEntScore] Sequence must at least contain %u nucs!",l5seq);


		if(STRING_ELT(pSeq,i)==NA_STRING)
			REAL(ans)[i]=NA_REAL;
		else
		{
			const char *inseq=CHAR(STRING_ELT(pSeq,i));

			// pos = 1-based position of last exonic nucleotide
			// Start position (pos-n5ExonNucs):
			// 			Includes last n5ExonNucs exonic nucleotides (when using 0-based array index)
			// Reads values for l5seq nuc's.
			for(j=0,k=pos-n5ExonNucs;j<l5seq;++j,++k)
			{
					seq[j]=ACGT[(unsigned)inseq[k]];
				if(seq[j]==zvl)
				{
					seq[0]=zvl;
					break;
				}
			}
			if(seq[0]!=zvl)
			{
				// numeric factors are statically used in this function (-> no variable)
				idx=seq[0]*4096+seq[1]*1024+seq[2]*256+seq[5]*64+seq[6]*16+seq[7]*4+seq[8];
				REAL(ans)[i]=log2(scoreconsensus5(seq[3],seq[4])*REAL(pMe2x5)[idx]);
			}
			else
				REAL(ans)[i]=NA_REAL;
		}
	}
	UNPROTECT(1);
	return ans;
}


SEXP maxent_seq_score3(SEXP pSeq, SEXP pFrame,SEXP pMeList)
{
	if(TYPEOF(pSeq)!=STRSXP)
		error("[maxent_seq_score3] pSeq must be STRING!");
	if(TYPEOF(pFrame)!=INTSXP)
		error("[maxent_seq_score3] pFrame must be INTEGER!");
	if(TYPEOF(pMeList)!=VECSXP)
		error("[maxent_seq_score3] pMeList must be VECSXP!");

	if(STRING_ELT(pSeq,0)==NA_STRING)
	{
		Rprintf("[maxent_seq_score3] pSeq is NA!");
		return R_NilValue;
	}

    ////////////////////////////////////////////////////////
    //
    //    pos=1-based last exon nuc = 3, l5seq=9
    //    |
    //  ATG | GTC | ATCGAA
    //  123   456   789012
    //    |     |
    //    |      end=6
    //    start=3
    //
    ////////////////////////////////////////////////////////

	unsigned i,j,k,l,pos,npos;
	const char     *inseq=CHAR(STRING_ELT(pSeq,0));
	const unsigned  start=INTEGER(pFrame)[0];          // 1-based position of leftmost  score3 position
	const unsigned  end  =INTEGER(pFrame)[1];          // 1-based position of rightmost score3
	const unsigned  nchar=strlen(inseq);

	if(end<=start)
		error("[maxent_seq_score3] end must be > start!");
	if(start<n3ExonNucs)
		error("[maxent_seq_score3] start must be greater than %u!",n3ExonNucs);

	pos=end-n3ExonNucs+l3seq; // 1-based index of last examined nuc when pos=end is scored
	if(nchar<pos)
		error("[maxent_seq_score3] length of inseq (%u) must be greater than %u (or reduce end)!",nchar,pos);

	// number of positions for which score5 can be calculated
	npos=end-start+1;

	SEXP ans;
	PROTECT(ans=allocVector(REALSXP,npos));

	SEXP pMe;				// R-list with pre-defined score values
	unsigned seq[l3seq];	// int-translated character values
	unsigned scidx;			// score-index, determines score to be read from pMe
	double   score[nsc3]; 	// stores pMe score result, finalscore is derived from entries
	double   finalscore;	// name equals return variable of function "maxentscore" in score3.pl

	for(pos=start,i=0; pos<=end; ++pos,++i)
	{
		// calculate score3
		// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
		// Translate sequence into numeric input values
		// Start position (pos-n3ExonNucs):
		// 			Includes last n3ExonNucs exonic nucleotides (when using 0-based array index)
		// Reads values for l5seq nucleotides (static string length).

		if(STRING_ELT(pSeq,i)==NA_STRING)
			REAL(ans)[i]=NA_REAL;
		else
		{
			const char* inseq=CHAR(STRING_ELT(pSeq,0));
			for(j=0,k=pos-n3ExonNucs;j<l3seq;++j,++k)
			{
				seq[j]=ACGT[(unsigned)inseq[k]];
				if(seq[j]==zvl)
				{
					seq[0]=zvl;
					break;
				}
			}
			if(seq[0]!=zvl)
			{
				// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
				// Calculate sequence score:
				for(j=0;j<nsc3;++j)
				{
					// calculate score-index from sequence
					// ("hashseq" function in score3.pl)
					pMe=VECTOR_ELT(pMeList,j);
					// score_columns values are stored in sc3c.
					scidx=0;
					for(k=0,l=(sc3l[j]-1);k<sc3l[j];++k,--l)
					scidx+=four3[l] * seq[sc3c[sc3s[j]+k]];

					// Read score value from input table
					score[j]=REAL(pMe)[scidx];
				}
				finalscore=(score[0]*score[1]*score[2]*score[3]*score[4])/
						(score[5]*score[6]*score[7]*score[8]);
				REAL(ans)[i]=log2(scoreconsensus3((unsigned char)seq[18],(unsigned char)seq[19])*finalscore);
			}
			else
				REAL(ans)[i]=NA_REAL;
		}
	}
	UNPROTECT(1);
	return ans;
}


SEXP maxent_score3(SEXP pSeq, SEXP pPos, SEXP pMeList)
{
	if(TYPEOF(pSeq)!=STRSXP)
		error("[maxent_score3] pSeq must be STRING!");
	if(TYPEOF(pMeList)!=VECSXP)
		error("[maxent_score3] pMeList must be VECSXP!");
	if(TYPEOF(pPos)!=INTSXP)
		error("[maxent_score3] pPos must be INT!");
	if(LENGTH(pSeq)!=LENGTH(pPos))
		error("[maxent_score3] LENGTH(pSeq) must equal LENGTH(pPos)!");

	int i,j,k,l,pos,n=LENGTH(pSeq);
	SEXP pMe;				// R-list with pre-defined score values
	unsigned seq[l3seq];	// int-translated character values
	unsigned scidx;			// score-index, determines score to be read from pMe
	double   score[nsc3]; 	// stores pMe score result, finalscore is derived from entries
	double   finalscore;	// name equals return variable of function "maxentscore" in score3.pl

	// Stored returned values for each sequence
	SEXP ans;
	PROTECT(ans=allocVector(REALSXP,n));
	for(i=0;i<n;++i)
	{
		// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
		// Check pos:
		// pos = 1-based position of last exonic nucleotide
		pos=INTEGER(pPos)[i];
		if(pos<n3ExonNucs)
			error("[score3.maxEntScore] pos must be >=%u! At least %u exon nucleotides needed!",n3ExonNucs,n3ExonNucs);

		// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
		// Translate sequence into numeric input values
		// Start position (pos-n3ExonNucs):
		// 			Includes last n3ExonNucs exonic nucleotides (when using 0-based array index)
		// Reads values for l5seq nucleotides (static string length).
		const char* inseq=CHAR(STRING_ELT(pSeq,i));
		for(j=0,k=pos-n3ExonNucs;j<l3seq;++j,++k)
		{
			seq[j]=ACGT[(unsigned)inseq[k]];
			if(seq[j]==zvl)
			{
				seq[0]=zvl;
				break;
			}
		}
		if(seq[0]!=zvl)
		{
			// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
			// Calculate sequence score:
			for(j=0;j<nsc3;++j)
			{
				// calculate score-index from sequence
				// ("hashseq" function in score3.pl)
				pMe=VECTOR_ELT(pMeList,j);
				// score_columns values are stored in sc3c.
				scidx=0;
				for(k=0,l=(sc3l[j]-1);k<sc3l[j];++k,--l)
					scidx+=four3[l] * seq[sc3c[sc3s[j]+k]];
				//if(scidx>=LENGTH(pMe))
				//	error("[maxent_score3] scidx error\n");

				// Read score value from input table
				score[j]=REAL(pMe)[scidx];
			}
			finalscore=(score[0]*score[1]*score[2]*score[3]*score[4])/
					(score[5]*score[6]*score[7]*score[8]);
			REAL(ans)[i]=log2(scoreconsensus3((unsigned char)seq[18],(unsigned char)seq[19])*finalscore);
		}
		else
			REAL(ans)[i]=NA_REAL;
	}
	UNPROTECT(1);
	return ans;
}


SEXP maxent_score2strand(SEXP pPscore, SEXP pMscore)
{
	if(TYPEOF(pPscore)!=REALSXP)
		error("[maxent_score2strand] pPscore must be REAL!");
	if(TYPEOF(pMscore)!=REALSXP)
		error("[maxent_score2strand] pMscore must be REAL!");
	if(LENGTH(pPscore)!=LENGTH(pMscore))
		error("[maxent_score2strand] pPscore and pMscore must have same length!");

	SEXP res;
	PROTECT(res=allocVector(INTSXP,LENGTH(pPscore)));
	unsigned i;
	for(i=0;i<LENGTH(pPscore);++i)
	{
		if(ISNA(REAL(pPscore)[i]) || ISNA(REAL(pMscore)[i]))
			INTEGER(res)[i]=strand_star;
		else
		{
			if(REAL(pPscore)[i]>= REAL(pMscore)[i])
				INTEGER(res)[i]=strand_plus;
			else
				INTEGER(res)[i]=strand_minus;
		}
	}

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Convert return type to factor

	SEXP levs;
	PROTECT(levs=allocVector(STRSXP,3));

	SET_STRING_ELT(levs,0,mkChar("+"));
	SET_STRING_ELT(levs,1,mkChar("-"));
	SET_STRING_ELT(levs,2,mkChar("*"));
	setAttrib(res,R_LevelsSymbol,levs);

	SEXP csymb;
	PROTECT(csymb=mkString("factor"));
	setAttrib(res,R_ClassSymbol,csymb);
	UNPROTECT(3);
	return res;
}

SEXP maxent_combine_strand(SEXP lStrand, SEXP rStrand)
{
	if(TYPEOF(lStrand)!=INTSXP)
		error("[maxent_combine_strand] lStrand must be INT!");
	if(TYPEOF(lStrand)!=INTSXP)
		error("[maxent_combine_strand] rStrand must be INT!");
	if(LENGTH(lStrand)!=LENGTH(rStrand))
		error("[maxent_combine_strand] LENGTH(lStrand) must equal LENGTH(rStrand)!");

	SEXP res;
	PROTECT(res=allocVector(INTSXP,LENGTH(lStrand)));
	unsigned i;
	for(i=0;i<LENGTH(lStrand);++i)
	{
		if(INTEGER(lStrand)[i]==strand_plus && INTEGER(rStrand)[i]==strand_plus)
			INTEGER(res)[i]=strand_plus;
		else
		{
			if(INTEGER(lStrand)[i]==strand_minus && INTEGER(rStrand)[i]==strand_minus)
				INTEGER(res)[i]=strand_minus;
			else
				INTEGER(res)[i]=strand_star;
		}
	}

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Convert res to factor
	SEXP levs;
	PROTECT(levs=allocVector(STRSXP,3));
	SET_STRING_ELT(levs,0,mkChar("+"));
	SET_STRING_ELT(levs,1,mkChar("-"));
	SET_STRING_ELT(levs,2,mkChar("*"));
	setAttrib(res,R_LevelsSymbol,levs);

	SEXP csymb;
	PROTECT(csymb=mkString("factor"));
	setAttrib(res,R_ClassSymbol,csymb);

	UNPROTECT(3);
	return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// 		hbond_score
//
///////////////////////////////////////////////////////////////////////////////////////////////////

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

static R_INLINE bool do_count_nMers(char const * const inseq, int * const array, unsigned const n, unsigned seqlen)
{
	// n = length of nMer
	if(n>16)
		error("[count_n_mers] n must be <= 16!");

	//Rprintf("[do_count_nMers] n=%u\n",n);

	// array length= 4^sqlen
	// unsigned nchar= 1 << (2*n);	// = 4^len

	// Points to start of n-mer
	char const * iter=inseq;
	unsigned array_idx=0;

	// iter traverses seq from *inseq for seqlen positions
	unsigned i,j,cval;
	// seqlen-n is the last start position for counting
	for(i=0;i<(seqlen-n+1);++i)
	{
		// Calculate array_index for each nMer at each seq-position
		// and increase value at index-position by 1
		//Rprintf("[count_n_mers] seq= '%s'\n",iter);
		for(j=0;j<n;++j)
		{
			// Throw error when traversing '\0' !
			if(iter[j]=='\0')
				error("[count_n_mers] Found string terminating '\0'!");
			cval=ACGT[(unsigned)iter[j]];
			if(cval==zvl)
				return false;

			//Rprintf("[count_n_mers] j: %u\tchar: '%c'\tcval: %u\n",j,iter[j],cval);
			// Index reverse order (because dna nMers do it this way)
			array_idx+= (cval << (2*((n-1)-j)));
		}
		//Rprintf("[count_n_mers] array_idx: %u\n",array_idx);
		++array[array_idx];
		++iter;
		array_idx=0;
	}

	return true;
}



SEXP count_nMers(SEXP pSeq, SEXP pN)
{
	if(TYPEOF(pSeq)!=STRSXP)
		error("[count_nMers] pSeq must be String!");
	if(TYPEOF(pN)!=INTSXP)
		error("[count_nMers] pN must be INT!");

	unsigned nSeq=LENGTH(pSeq);
	int n=INTEGER(pN)[0];

	///////////////////////////////////////////////////////////////////////////
	// Create data.frame
	unsigned nCols=nSeq;
	unsigned nProtected=0;
	//Rprintf("[count_nMers] nSeq=%u\n",nSeq);

	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;

	//SET_VECTOR_ELT(dflist,0,vov);

	// array length= 4^sqlen
	// unsigned nchar= 1 << (2*n);	// = 4^len
	unsigned array_length= 1 << (2*n); // nRows

	//Rprintf("[count_nMers] array_length=%u\n",array_length);

	// For each sequence

	unsigned i;
	for(i=0;i<nSeq;++i)
	{
		char const * const seq=CHAR(STRING_ELT(pSeq,i));
		// A) Count nchar
		unsigned nchar=strlen(seq);
		// B) do_count_nMers in given seq
		SEXP colvec;
		PROTECT(colvec=allocVector(INTSXP,array_length));
		++nProtected;

		memset(INTEGER(colvec),0,array_length*sizeof(int));

		// C) Add column to result data.frame
		//char const * const inseq, unsigned * const array, unsigned const n, unsigned seqlen)
		if(!do_count_nMers(seq,INTEGER(colvec),n,nchar))
			error("[count_nMers] character mismatch!");

		SET_VECTOR_ELT(dflist,i,colvec);
	}

	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	int buf_size=1024;
	char *buf=R_alloc(buf_size,sizeof(char));
    for(i=0;i<nCols;++i)
    {
    	sprintf(buf,"%i",i);
    	SET_STRING_ELT(col_names,i,mkChar(buf));
    }
	setAttrib(dflist,R_NamesSymbol,col_names);


	// Add sequence chars as row.names
    setAttrib(dflist,R_RowNamesSymbol,create_dna_n_mers(pN));

    // Make list to data.frame
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	UNPROTECT(nProtected);

	return dflist;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //


static const unsigned n_seqIdx=9;
static const unsigned seq_idx[]	={0,1,2,5,6,7,8,9,10};
static const unsigned hbs_shift[]={16,14,12,10,8,6,4,2,0};

static R_INLINE double do_calc_hbs(const char * const inseq, const double * const hbs)
{
	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Calculate array index in hbs from sequence
	// 3 exon nucs and 6 intron nucs (GT excluded)
	unsigned j,cval,hbs_idx=0;
	for(j=0;j<n_seqIdx;++j)
	{
		// Look at positions 0,1,2,5,6,...
		// and translate into int value with ACTG array
		cval=ACGT[(unsigned)inseq[seq_idx[j]]];
		if(cval==zvl)
			return 0;
		// Insert bit-mask into hbs_idx
		hbs_idx+=cval << hbs_shift[j];
	}
	// H-Bond score
	return hbs[hbs_idx];
}

SEXP hbond_score(SEXP pSeq,SEXP pHb)
{
	const unsigned nh=262144;
	const unsigned ns=LENGTH(pSeq);
	unsigned i;
	SEXP res;

	if(TYPEOF(pSeq)!=STRSXP)
		error("[hbond_score] pSeq must be STRING!");
	if(TYPEOF(pHb)!=REALSXP)
		error("[hbond_score] pHb must be REAL!");

	if(LENGTH(pHb)!=nh)
		error("[hbond_score] pHb must have length 262144!");
	PROTECT(res=allocVector(REALSXP,ns));

	for(i=0;i<ns;++i)
	{
		if(STRING_ELT(pSeq,i)==NA_STRING)
			REAL(res)[i]=0;
		else
		{
			const char *inseq=CHAR(STRING_ELT(pSeq,i));
			// Score 0 when first exon nucs unequal 'GT'
			if((ACGT[(unsigned)inseq[3]]!=g_val) | (ACGT[(unsigned)inseq[4]]!=t_val))
				REAL(res)[i]=0;
			else
				REAL(res)[i]=do_calc_hbs(inseq,REAL(pHb));
		}
	}
	UNPROTECT(1);
	return res;
}


// Creates string vector with all ATCG combinations for given
// string length (up to 32)
SEXP create_dna_n_mers(SEXP pLen)
{
	if(TYPEOF(pLen)!=INTSXP)
		error("[create_dna_n_mers] pLen must be INT!");
	if(LENGTH(pLen)!=1)
		error("[create_dna_n_mers] Length(pLen) must be 1!");

	int len=INTEGER(pLen)[0];

	if(len>32)
		error("[create_dna_n_mers] Maximum value for pLen is 32!");
	if(len<=0)
		error("[create_dna_n_mers] pLen must be > 0!");

	unsigned nstr= 1 << (2*len);	// = 4^len

	unsigned long long i;
	unsigned j,k;

	char *c = (char*) R_alloc(len+1,sizeof(char));
	c[len]='\0';

	SEXP res;
	PROTECT(res=allocVector(STRSXP,nstr));

	for(i=0;i<nstr;++i)
	{
		for(j=0,k=(2*(len-1));j<len;++j)
		{
			c[j]=(char) rev_ACGT[(i >> k) & 3];
			k-=2;
		}
		SET_STRING_ELT(res,i,mkChar(c));
	}
	UNPROTECT(1);
	return res;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
//
// alt_group
// Calculate alternative gap-site (splice-site) values
// alt_groups share identical sequence- and group-values and have pairwise different alt-values.
//
// The function expects that input values are sorted by pSeq and pGroup and pAlt
// A) When one of both values changes, a new alt-group is defined
// B) All precedent and successive group-identical rows are treated as *one Group*.
// C) The function defines for each group more than one element
// 				altid  :   unique id value for each grop
//				gapdiff:   Positive difference of gaplen between presend and next lower gaplen
//				altrank:   Positive rank for gaplen value
// 				nalt   :   Number of different alt-positions in group
//
// D) For singular rows, the aforementioned values are set to 0
//
///////////////////////////////////////////////////////////////////////////////////////////////////

SEXP alt_group(SEXP pId, SEXP pSeq, SEXP pGroup, SEXP pAlt)
{
	if(TYPEOF(pId)!=INTSXP)
		error("[alt_group] pInt must be INT!");
	if(TYPEOF(pSeq)!=INTSXP)
		error("[alt_group] pSeq must be INT!");
	if(TYPEOF(pGroup)!=INTSXP)
		error("[alt_group] pGroup must be INT!");
	if(TYPEOF(pAlt)!=INTSXP)
		error("[alt_group] pAlt must be INT!");

	unsigned nRows=LENGTH(pId);
	if(LENGTH(pSeq)!=nRows)
		error("[alt_group] length(pSeq) must equal length(pId)!");
	if(LENGTH(pGroup)!=nRows)
		error("[alt_group] length(pGroup) must equal length(pId)!");
	if(LENGTH(pAlt)!=nRows)
		error("[alt_group] length(pAlt) must equal length(pId)!");

	if(nRows<2)
		return R_NilValue;

	int *id=INTEGER(pId);
	int *seq=INTEGER(pSeq);
	int *group=INTEGER(pGroup);
	int *alt=INTEGER(pAlt);

	// create data.frame
	unsigned nProtected=0;
	unsigned nCols=5;

	// Column 0: id
	SEXP id_vec;
	PROTECT(id_vec=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 1: alt_id
	SEXP altid_vec;
	PROTECT(altid_vec=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 2: alt_rank
	SEXP altrank_vec;
	PROTECT(altrank_vec=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 3: gap_diff
	SEXP gapdiff_vec;
	PROTECT(gapdiff_vec=allocVector(INTSXP,nRows));
	++nProtected;

	// Column 4: n_alt
	SEXP nalt_vec;
	PROTECT(nalt_vec=allocVector(INTSXP,nRows));
	++nProtected;

	int i,j,pos,nalt=1; // must be signed
	unsigned altid=1;
	unsigned altrank=1;
	//unsigned nSites=1;

	// + + + + + + + + + + + + + + + + + //
	// Set values for first row
	INTEGER(id_vec)[0]=id[0];
	INTEGER(altid_vec)[0]=altid;
	INTEGER(altrank_vec)[0]=altrank;
	INTEGER(gapdiff_vec)[0]=0;

	// + + + + + + + + + + + + + + + + + //
	// Cycle rows
	for(i=1;i<nRows;++i)
	{
		// Copy id value
		INTEGER(id_vec)[i]=id[i];

		// Last and current row belong to same group
		if((seq[i]==seq[i-1]) & (group[i]==group[i-1]))
		{
			++altrank;
			++nalt;
			INTEGER(altid_vec)[i]=altid;
			INTEGER(altrank_vec)[i]=altrank;
			INTEGER(gapdiff_vec)[i]=alt[i]-alt[i-1];

		}
		else
		{
			// + + + + + + + + + + + + + + + + + //
			// Group change:
			if(nalt==1)
			{
				// A) Precedent row is singular
				// No new altid value is created
				// Row values are set to 0
				pos=i-1;
				INTEGER(altid_vec)[pos]=0;
				INTEGER(altrank_vec)[pos]=0;
				INTEGER(gapdiff_vec)[pos]=0;
				INTEGER(nalt_vec)[pos]=0;
			}
			else
			{
				// B) Precedent row is part of group block
				// Create new altid value
				// Step back for inserting n_alt values
				pos=i-1;
				for(j=pos;j>(pos-nalt);--j)
					INTEGER(nalt_vec)[j]=nalt;
				++altid;
			}

			// C) Prepare for first row of next group
			// Set current row values
			//++nSites;
			nalt=1;
			altrank=1;
			INTEGER(altid_vec)[i]=altid;
			INTEGER(altrank_vec)[i]=altrank;
			INTEGER(gapdiff_vec)[i]=0;
		}
	}

	// + + + + + + + + + + + + + + + + + //
	// Process last group
	if(nalt==1)
	{
		// A) Last row is singular
		pos=i-1;
		INTEGER(altid_vec)[pos]=0;
		INTEGER(altrank_vec)[pos]=0;
		INTEGER(gapdiff_vec)[pos]=0;
		INTEGER(nalt_vec)[pos]=0;
	}
	else
	{
		// B) Last row is part of group block
		// Step back for inserting n_alt values
		--i;
		for(j=i;j>(i-nalt);--j)
			INTEGER(nalt_vec)[j]=nalt;
	}


	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	SEXP dflist;
	PROTECT(dflist=allocVector(VECSXP,nCols));
	++nProtected;

	SET_VECTOR_ELT(dflist,0,id_vec);
	SET_VECTOR_ELT(dflist,1,altid_vec);
	SET_VECTOR_ELT(dflist,2,altrank_vec);
	SET_VECTOR_ELT(dflist,3,gapdiff_vec);
	SET_VECTOR_ELT(dflist,4,nalt_vec);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Column Names
	SEXP col_names;
	PROTECT(col_names=allocVector(STRSXP,nCols));
	++nProtected;

	SET_STRING_ELT(col_names,0,mkChar("id"));
	SET_STRING_ELT(col_names,1,mkChar("alt_id"));
	SET_STRING_ELT(col_names,2,mkChar("diff_ranks"));
	SET_STRING_ELT(col_names,3,mkChar("gap_diff"));
	SET_STRING_ELT(col_names,4,mkChar("nr_alt"));
	setAttrib(dflist,R_NamesSymbol,col_names);

	// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
	// Row Names
	SEXP row_names;
    PROTECT(row_names=allocVector(STRSXP,nRows));
    ++nProtected;
	char *buf=R_alloc(buf_size,sizeof(char));
    for(i=0;i<nRows;++i)
    {
    	sprintf(buf,"%i",i);
    	SET_STRING_ELT(row_names,i,mkChar(buf));
    }
    //free(buf);

    setAttrib(dflist,R_RowNamesSymbol,row_names);
	setAttrib(dflist,R_ClassSymbol,mkString("data.frame"));
	UNPROTECT(nProtected);
	//Rprintf("[alt_group] Found %u sites.\n",nSites);
	return dflist;

}

SEXP get_alph_index(SEXP pIdxLen, SEXP pAlpha)
{
	if(TYPEOF(pIdxLen)!=INTSXP)
		error("[get_alph_index] pIdxLen must be INT!");
	if(TYPEOF(pAlpha)!=STRSXP)
		error("[get_alph_index] pAlpha must be STR!");
	if(LENGTH(pAlpha)>1)
		error("[get_alph_index] pAlpha must have length 1!");

	const char *alpha=CHAR(STRING_ELT(pAlpha, 0));
	int * idx_len=INTEGER(pIdxLen);
	unsigned n_char;
	char *idx=create_idx(*idx_len, alpha,strlen(alpha), &n_char);

	SEXP ans;
	PROTECT(ans=allocVector(STRSXP,*idx_len));
	unsigned i, str_diff=n_char+1;

	char *iter=idx;
	for(i=0;i<(*idx_len);++i)
	{
		SET_STRING_ELT(ans,i,mkCharLen(iter,n_char));
		iter+=str_diff;
	}
	UNPROTECT(1);
	return ans;
}


void R_init_spliceSites(DllInfo *info)
{
	R_CallMethodDef cmd[] ={
			{ "trunc_pos", 				(DL_FUNC) &trunc_pos,				4},
			{ "silic_tryp_pos",			(DL_FUNC) &silic_tryp_pos,			3},
			{ "maxent_score5",			(DL_FUNC) &maxent_score5,			3},
			{ "maxent_score3",			(DL_FUNC) &maxent_score3,			3},
			{ "maxent_seq_score5",		(DL_FUNC) &maxent_seq_score5,		3},
			{ "maxent_seq_score3",		(DL_FUNC) &maxent_seq_score3,		3},
			{ "maxent_score2strand",	(DL_FUNC) &maxent_score2strand,		2},
			{ "maxent_combine_strand",	(DL_FUNC) &maxent_combine_strand,	2},
			{ "alt_group",				(DL_FUNC) &alt_group,				4},
			{ "get_alph_index",			(DL_FUNC) &get_alph_index,			2},
			{ "hbond_score",            (DL_FUNC) &hbond_score,             2},
			{ "create_dna_n_mers",		(DL_FUNC) &create_dna_n_mers,       1},
			{ "count_nMers",			(DL_FUNC) &count_nMers,				2},
			{NULL, NULL, 0}
	};
	//			{ "",	(DL_FUNC) &,	}
	R_registerRoutines(info, NULL, cmd, NULL, NULL);
}

#endif /* SPLICESITES_C_ */

