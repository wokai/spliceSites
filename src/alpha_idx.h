/*
 * alpha_idx.h
 *
 *  Created on: 15.08.2013
 *      Author: kaisers
 */

#ifndef ALPHA_IDX_H_
#define ALPHA_IDX_H_

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// 	Definition: str_chain's:
//	Multiple equally sized (and '\0' separated) strings
//	combined together in one char-array.
// 	e.g. "aa\0ab\0ba\0ab\0"
//
//	Usage:
//	char *res=chop_str("abc",3);
//	Result: "a\0b\0c\0"
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

char * chop_str(const char * str, unsigned str_len)
{
	char* res=R_alloc(2*str_len,sizeof(char));
	unsigned i,j;
	for(i=0,j=0;i<str_len;++i,++j)
	{
		res[j]=str[i];
		++j;
		res[j]='\0';
	}
	return res;
}


// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
// cmb_str:
// 	Combine each char in lhs with rhs entry as suffix.
// 	Input format	: lhs="ab\0", rhs="aa\0ab\0ba\0ab\0" (str_chain)
// 	rn_char		: Word length in rhs
//	n_str		: Number of words in rhs
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

void cmb_str(const char *lhs,unsigned llen, char **rhs,unsigned *rn_char, unsigned *rn_str)
{
	char *res=R_alloc(((*rn_char)+2)*llen*(*rn_str),sizeof(char));
	unsigned i,j;
	char *res_iter=res;
	const char *lhs_iter=lhs, *rhs_iter=*rhs;

	for(i=0;i<llen;++i)
	{
		for(j=0;j<(*rn_str);++j)
		{
			// Prefix: lhs
			strncpy(res_iter,lhs_iter,1);
			++res_iter;
			// Suffix: rhs
			strncpy(res_iter,rhs_iter,*rn_char);
			res_iter+=((*rn_char)+1);
			rhs_iter+=((*rn_char)+1);
		}
		++lhs_iter;
		rhs_iter=*rhs;
	}

	// Update values in rhs-pointers
	++(*rn_char);
	(*rn_str)*=llen;
	*rhs=res;
}

// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //
//
// create_idx:	Creates str_chain of given length 
//		from combinations of given alphabet with unique entries.
//
// + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + //

char * create_idx(int idx_len, const char *alpha,unsigned alpha_len, unsigned *mx_nchar)
{
	*mx_nchar=ceil(log(idx_len)/log(alpha_len));
	char *res;
	unsigned n_char, n_str,i;

	res=chop_str(alpha,alpha_len);
	n_char=1;
	n_str=alpha_len;

	if((*mx_nchar)>1)
	{
		char **cmb=(char**)R_alloc(1,sizeof(char*));
		*cmb=res;
		for(i=1;i<(*mx_nchar);++i)
			cmb_str(alpha,alpha_len,cmb,&n_char,&n_str);

		res=*cmb;
		*cmb=0;
	}
	// Trunc string array to appropriate size:
	unsigned ret_size=idx_len*(n_char+1);
	char *ret=R_alloc(ret_size,sizeof(char));
	memcpy(ret,res,ret_size);
	return ret;
}


#endif /* ALPHA_IDX_H_ */
