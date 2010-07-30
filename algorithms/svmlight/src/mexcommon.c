/************************************************************************/
/*                                                                      */
/*   mexcommon.c                                                        */
/*                                                                      */
/*   Common definitions and functions in support of the MATLAB/MEX      */
/*                                                                      */
/*   Author: Tom Briggs                                                 */
/*   Date: January 1, 2005                                              */
/*                                                                      */
/*                                                                      */
/* Based on SVMLite:                                                    */
/*   Copyright (c) 2004  Thorsten Joachims - All rights reserved        */
/*                                                                      */
/*   This software is available for non-commercial use only. It must    */
/*   not be modified and distributed without prior permission of the    */
/*   author. The author is not responsible for implications from the    */
/*   use of this software.                                              */
/*                                                                      */
/*   January 1, 2005 - Modifications by Tom Briggs for MATLAB           */
/************************************************************************/
#ifdef MATLAB_MEX    /* if not in MATLAB - skip the whole file */

#include <stdio.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"
#include "math.h"
#include "svm_common.h"
#include "svm_learn.h"
#include "global.h"

#include "mexsvmlearn.h"
#include "mexcommon.h"


/**
 * storeValue( ) - store the given double value in the mxStruct
 * structure at field fieldName */
void storeValue(mxArray *mxStruct, char *fieldName, double value)
{
  char msg[1025];
  mxArray *fieldValue; 

  int fieldNum = mxGetFieldNumber(mxStruct, fieldName); 
  if (fieldNum < 0) {
    sprintf(msg,"Error: unknown field: (%s)\n", fieldName);
    mexErrMsgTxt(msg);
  }

  
  fieldValue = mxCreateDoubleMatrix(1,1,mxREAL); 
  *mxGetPr(fieldValue) = value; 
  mxSetFieldByNumber(mxStruct, 0, fieldNum, fieldValue); 

}

/**
 * storeString( ) - store the given string value in the 
 * mxStruct structure in field named fieldName 
 */
void storeString(mxArray *mxStruct, char *fieldName, char *value)
{
  char msg[1025];
  mxArray *fieldValue;
  int fieldNum = mxGetFieldNumber(mxStruct, fieldName); 
  if (fieldNum < 0) {
    sprintf(msg,"Error: unknown field: (%s)\n", fieldName);
    mexErrMsgTxt(msg);
  }


  fieldValue = mxCreateString(value); 
  mxSetFieldByNumber(mxStruct, 0, fieldNum, fieldValue); 

}

/**
 * storeArray( ) - store the given array double array with n elements
 * into the structure mxStruct at field named name.
 */
void storeArray(mxArray *mxStruct, char *name, double *srcArray, int n)
{
  mxArray *fieldValue;
  double *tmpArray;
  int i;

  fieldValue = mxCreateDoubleMatrix(n,1,mxREAL); 
  tmpArray = mxGetPr(fieldValue);

  for (i = 0; i < n; i++) 
    tmpArray[i] = srcArray[i];

  mxSetField(mxStruct, 0, name, fieldValue); 
}

/**
 * storeArrayLong( ) - store the given long array with n elements
 * into the structure mxStruct at field named name.
 */
void storeArrayLong(mxArray *mxStruct, char *name, long *srcArray, int n)
{
  mxArray *fieldValue;
  double *tmpArray;
  int i;

  fieldValue = mxCreateDoubleMatrix(n,1,mxREAL); 
  tmpArray = mxGetPr(fieldValue);

  for (i = 0; i < n; i++)
    tmpArray[i] = (double) srcArray[i];
  mxSetField(mxStruct, 0, name, fieldValue); 
}

/**
   TODO add support for sparse arrays  - MATLAB has sparse matrices
    that are similar to SVM's internal format for docs.  This method
    doesn't work with sparse vectors.
**/

/**
 * storeDocs - store the given DOC named docs with ndocs entries of nwords
 * width in the structure named mxStruct at field named name.
 */
 void storeDocs(mxArray *mxStruct, char *name, DOC **docs, int ndocs, int nwords)
 {
 	if (input_is_sparse)
 	 storeDocsSparse(mxStruct, name, docs, ndocs, nwords);
 	else
 	 storeDocsFull(mxStruct,name,docs,ndocs,nwords);
 }
 
 
 
 WORD *getWords(DOC **docs, int row)
 {
 	DOC *doc;
 	
 	
	if (docs == NULL) return NULL;
	
	doc = docs[row];
	
	if (doc == NULL) return NULL;
	if (doc->fvec == NULL) return NULL;
	if (doc->fvec->words == NULL) return NULL;
	return docs[row]->fvec->words;
 }
 
 
 int numOccupied(DOC **docs, int ndocs)
 {
 	int total = 0;
 	int i,j ;
 	WORD *words;
 	DOC *doc;

	if (docs == NULL) return 0;
	
 	for (i = 0; i < ndocs; i++) {
 			

		doc = docs[i];
		if (doc == NULL) continue;
		
		if (doc->fvec == NULL) continue;
		if (doc->fvec->words == NULL) continue;
		
		words = doc->fvec->words;			 	

	    for (j = 0; words[j].wnum != 0; j++) 
 			total++; 
 	}
 	
 	return total;
 }
 
  
 void storeDocsSparse(mxArray *mxStruct, char *name, DOC **docs, int ndocs, int nwords)
 {
  mxArray *fieldValue;
  double *sr;
  int i,j,k, count,*irs,*jcs, numOcc;
  SVECTOR *fvec;
  WORD *words;
  DOC *doc;

	if (docs == NULL) return;

	numOcc = numOccupied(docs,ndocs); 
	printf("Using %d occupied cells\n");
	fieldValue = mxCreateSparse(ndocs, nwords, numOcc,mxREAL);
  
  sr = mxGetPr(fieldValue);
  irs = mxGetIr(fieldValue);
  jcs = mxGetJc(fieldValue);

  count = 0;					/* track # of items added */
  for(i = 0; i < nwords; i++)   /* rank by column/word */
  {
  	for (j = 0; j < ndocs; j++)  /* handle each row */
  	{
  		WORD *w = getWords(docs,j);	/* get the word list for this row */
  		
  		if (w == NULL)  continue;   /* skip if empty */
  		
	  	for (k = 0; w[k].wnum != 0; k++) /* iterate over each word */
  		{   
  			/* this row contains the word indicated by the outermost loop */
  			if ((w[k].wnum ) == i) 
  			{
  				sr[count] = w[k].weight;
  				irs[count] = j;	
				count++;
  			}
  		} /* end handle word array */

  	} /* end handle row */
    jcs[i+1] = count;
	jcs[i] = count;
  }
  		
    mxSetField(mxStruct,0, name, fieldValue);
}


 void storeDocsFull(mxArray *mxStruct, char *name, DOC **docs, int ndocs, int nwords)
 {
  mxArray *fieldValue;
  double *tmpArray;
  int i,j;
  SVECTOR *fvec;
  WORD *words;
  DOC *doc;

  if (docs == NULL) return;

  fieldValue = mxCreateDoubleMatrix(ndocs, nwords, mxREAL);
  tmpArray = mxGetPr(fieldValue);

	/* zero out tmparray */
	for (i = 0; i < (ndocs*nwords); i++)
		tmpArray[i] = 0;
		
		


  for (i = 0; i < ndocs; i++) {
    doc = docs[i];
    if (doc == NULL) continue;

    fvec = doc->fvec;
    if (fvec == NULL) continue;

    words = fvec->words;
    if (words == NULL) continue;

    for (j = 0; words[j].wnum != 0; j++) {
      WORD w = words[j];
      int index = computeOffset(ndocs, nwords, i,(w.wnum-1));
      tmpArray[index] = w.weight;
    }
  }

  mxSetField(mxStruct,0, name, fieldValue);
  

}


/**
 * store_kern_parms() - store the given kernel parameters in kerenel_parm
 * into the given mxStruct structure.
 */
void store_kern_parms(mxArray *mxStruct, char *name,KERNEL_PARM *kernel_parm)
{
  const char *fields[] = { "kernel_type", "poly_degree", "rbf_gamma",
			   "coef_lin", "coef_const", "custom" };

  mxArray *fieldStruct = mxCreateStructMatrix(1,1,6,fields);
  storeValue(fieldStruct, "kernel_type", kernel_parm->kernel_type);
  storeValue(fieldStruct, "poly_degree", kernel_parm->poly_degree);
  storeValue(fieldStruct, "rbf_gamma", kernel_parm->rbf_gamma);
  storeValue(fieldStruct, "coef_lin", kernel_parm->coef_lin);
  storeValue(fieldStruct, "coef_const", kernel_parm->coef_const);
  storeString(fieldStruct, "custom", kernel_parm->custom);

  mxSetField(mxStruct,0,name, fieldStruct);
}
  



/**
 * store_model() - store the given model into the given mxStruct structure
 */
void store_model(struct model *model, mxArray *mxOut[] )
{
  int i = 0;
  int numfields;
  int dims[] = { 1,1 };
  mxArray *mxStruct;
  
  const char *fields[] = { "sv_num", "upper_bound", "b", "totwords", 
			   "totdoc", "loo_error", "loo_recall", "loo_precision", 
			   "xa_error","xa_recall", "xa_precision", "maxdiff", 
			   "r_delta_sq", "r_delta_avg", "model_length", "loss", 
			   "vcdim", "alpha", "lin_weights", "index",
			   "supvec" , "kernel_parm","example_length","a" };
  
  numfields = 24;
  mxStruct = mxCreateStructMatrix(1,1, numfields, fields);
  
  storeValue(mxStruct, "sv_num", model->sv_num);
  storeValue(mxStruct, "upper_bound", model->at_upper_bound);
  storeValue(mxStruct, "b", model->b);
  storeValue(mxStruct, "totwords", model->totwords);
  storeValue(mxStruct, "totdoc", model->totdoc);
  storeValue(mxStruct, "loo_error", model->loo_error);
  storeValue(mxStruct, "loo_recall", model->loo_recall);
  storeValue(mxStruct, "loo_precision", model->loo_precision);
  storeValue(mxStruct, "xa_error", model->xa_error);
  storeValue(mxStruct, "xa_recall", model->xa_recall);
  storeValue(mxStruct, "xa_precision", model->xa_precision);
  storeValue(mxStruct, "maxdiff", model->maxdiff);
  storeValue(mxStruct, "r_delta_sq", model->r_delta_sq);
  storeValue(mxStruct, "r_delta_avg", model->r_delta_avg);
  storeValue(mxStruct, "model_length", model->model_length);
  storeValue(mxStruct, "example_length", model->example_length);
  storeValue(mxStruct, "loss", model->loss);
  storeValue(mxStruct, "vcdim", model->vcdim);
  
  storeArray(mxStruct, "alpha", model->alpha, model->totdoc);
  storeArray(mxStruct, "a", model->a, model->totdoc);

  if (model->lin_weights != NULL)
    storeArray(mxStruct, "lin_weights", model->lin_weights, model->totdoc+2);
  else {
  	if (verbosity > 0) 
	  	printf("Skipping lin_weights (array is empty)\n");
    }

  if (model->index != NULL)
    storeArrayLong(mxStruct, "index", model->index, (int) model->totdoc);
  else {
  	if (verbosity > 0)
  	 printf("Skipping model->index (array is empty)\n");
  }

  if (model->supvec != NULL)
     storeDocs(mxStruct,"supvec", model->supvec, model->sv_num, model->totwords);
  else {
  	if (verbosity > 0) printf("Skipping supvec (array is empty)\n");
  }
  
  store_kern_parms(mxStruct, "kernel_parm", &(model->kernel_parm));

  mxOut[0] = mxStruct;


}


/**
 * mexToDOC() - convert the MATLAB/MEX array mxData and mxLabels into 
 * the SVMLite formatted docs array and label array.  Note that 
 * MATLAB uses column major ordering, and SVMLite (and most programs)
 * use row major ordering.  This method unravels that complication.
 */
void mexToDOC(mxArray *mxData, mxArray *mxLabels, DOC ***docs, double **label, 
		    long int *totwords, long int *totdoc)
{
  int i;
  int rows, cols;
  double *yvals;
  WORD *words;

  if (mxData == NULL) {
    printf("WARNING: mexToDoc : mxData is NULL");
    return;
  }

  /* retrieve the rows and columns */
  rows = mxGetM(mxData);
  cols = mxGetN(mxData);

  /* allocate memory for the DOC rows */
  (*docs) = (DOC **)my_malloc(sizeof(DOC *) * rows);

  /* allocate memory for the labels */
  if (mxLabels != NULL)
    (*label) = (double *)my_malloc(sizeof(double)* rows);

  /* allocate a single buffer in memory for the words (hold n columns) */
  words = (WORD *)my_malloc(sizeof(WORD)*cols);

  /* store the number of words and docs */
  if ((totwords != NULL) && (totdoc != NULL)) {
    (*totwords) = cols;
    (*totdoc) = rows;
  }

  /* load the yvals from the mxLabels array */
  if (mxLabels != NULL)
    yvals = mxGetPr(mxLabels);
  
  /* for each row, create a corresponding vector of *DOC and store it */
  for (i = 0; i < rows; i++) {
    SVECTOR *fvec = NULL;
    int j;

	/* parse and copy the given mxData into a WORD array words */
    parse_mxEntry(i, mxData, mxLabels, words);

	/* create the intermediate structure (svector in svm_common.c) */
    fvec = create_svector(words,"",1.0);


	for (j = 0; j < 2; j++) {
		(*docs)[i] = create_example(i, 0, 0, 1.0, fvec);
		if (mxLabels != NULL)
		(*label)[i] = yvals[i];
	}

  }

}


/**
 * parse_mxEntry() - parse the given array mxData into the words WORD array
 */
int parse_mxEntry(int row, mxArray *mxData, mxArray *mxLabels, WORD *words)
{
	
	if (mxIsSparse(mxData))
		return parse_mxSparseEntry(row, mxData,mxLabels, words);
	else
		return parse_mxFullEntry(row,mxData,mxLabels, words);
}


int parse_mxFullEntry(int row, mxArray *mxData, mxArray *mxLabels, WORD *words)
{
  int i;
  int cols,rows;
  int index;
  double *data;

  data = mxGetPr(mxData);

  cols = mxGetN(mxData);
  rows = mxGetM(mxData);

  /* here we assume that the mxData is full, not sparse,
   * and relying on that assumption, we assign the word number
   * to be the position within the full array */
  for (i = 0; i < cols; i++)
    {
      index = computeOffset(rows, cols, row, i);

      (words[i]).wnum=(i+1);  /* assign the word number  */
      (words[i]).weight=(CFLOAT)data[index];

    }

  /* assign the last word to be 0 (null) to note the end of the list */
  (words[i]).wnum = 0;
	
  return 0;
}

/**
 * Extract a single row of entries out of the sparse matrix.
 * This is a little ugly, since MATLAB uses column major ordering,
 * and we have to stride across all columns, then iterate over each row
 * to find if a particular row holds a particular column. */
 
 /** this code is based on examples from explore.c */
 
int parse_mxSparseEntry(int row, mxArray *mxData, mxArray *mxLabels, WORD *words)
{
  int col,n,idx,total;
  int start_row_index, end_row_index;
  int curr_row;
  double *pr;
  int *ir, *jc;
  
     input_is_sparse = 1;
  
	pr = mxGetPr(mxData);
	ir = mxGetIr(mxData);
	jc = mxGetJc(mxData);


	n = mxGetN(mxData);
	total = 0;
	idx = 0; /* word vector index  */
	
 	for (col = 0; col < n; col++)
    {
		start_row_index = jc[col];
		end_row_index = jc[col+1];
		
		if (start_row_index == end_row_index) 
			continue; /* no rows for this column */
			
		else {
			for (curr_row = start_row_index; curr_row < end_row_index; curr_row++)
			{
				if ((ir[curr_row]) == row) { /* if this is the desired row... */
					(words[idx]).wnum = col+1;
					(words[idx]).weight = (CFLOAT)pr[total];
					idx++;  /* increment pointer within words array */
				}
				total++;				
			}
		}
    }

  /* assign the last word to be 0 (null) to note the end of the list */
  (words[idx]).wnum = 0;
}



/**
 * computeOffset() - given an i (row) and j (col) in row-major order,
 *  compute the corresponding index into a column-major single dimensioned
 *  array.
 */
int computeOffset(int numRows, int numCols, int i, int j)
{
  int index = (numRows * j) + i;
  if (index >= numRows * numCols) 
      mexErrMsgTxt("Error: index exceeds array length"); 

  return index;
}


/**
 * myBzero() - clear out an area of memory by zero'ing out the 
 * range of values.  
 */
void myBzero(char *ptr, int numbytes)
{
  int i;
  if (ptr == NULL) return;
  if (numbytes == 0) return;

  for (i = 0; i < numbytes; i++)
    {
      ptr[i] = 0x00;
    }
}
 
/**
 * restoreArray() - given a structure in mxStruct, retrieve a double
 * array from the field named "name".
 * TODO handle sparse arrays
 */
double *restoreArray(const mxArray *mxStruct, char *name)
{
  char msg[255];
  mxArray *mxValue;

  /* retrieve the array element from the struct */
  if (! mxIsStruct(mxStruct)) {
    sprintf(msg,"mxStruct is not a struct, it is a %s\n", 
	    mxGetClassName(mxStruct));
    mexErrMsgTxt(msg);
  }

  mxValue = mxGetField(mxStruct, 0, name);

  if (mxValue == NULL) {
    /*printf("Warning: could not retrieve %s from the given struct.\n", name); */
    return NULL;
  }

  return mxGetPr(mxValue); 
}
 
/**
 * restoreArrayLong() - given a structure in mxStruct, retrieve
 * a long integer array from the field named name in mxStruct.
 */
long *restoreArrayLong(const mxArray *mxStruct, char *name)
{
  /* retrieve the array element from the struct */
  long *lArray;
  double *dArray;
  int m,i;

  mxArray *mxValue = mxGetField(mxStruct, 0, name);

  if (mxValue == 0) {
    printf("Warning: could not retrieve %s from the given struct.\n", name);
    return NULL;
  }

   m = mxGetM(mxValue);

   lArray = (long *) my_malloc(m * sizeof(long));
   dArray = mxGetPr(mxValue);

   if ((lArray == NULL) || (dArray == NULL))
     mexErrMsgTxt("Error allocating lArray or dArray");

   for (i = 0; i < m; i++) {
     lArray[i] = (long) dArray[i];
   }

   return lArray;
}
  

/**
 * restoreValue() - return a the given double value from the mxStructure
 * from field named name.
 */
double restoreValue(const mxArray *mxStruct, char *name)
{
  double *array = restoreArray(mxStruct,name);
  if (array == NULL) 
    mexErrMsgTxt("Error: restoreArray returned NULL array to restoreValue");

  return array[0];
}


/**
 * restore_model() - construct an SVM "MODEL" structure from the given
 * MEX structure.
 */
MODEL *restore_model(const mxArray *mxStruct )
{
  char strPtr[KPARM_CUSTOM_LEN],*buff;
  long words, doc; 
  mxArray *kstruct,*kstring,*svec;
  MODEL *model = NULL;
  KERNEL_PARM *kparm;
  int buflen;

  model = (MODEL *)my_malloc(sizeof(MODEL)); 

  model->sv_num = (long) restoreValue(mxStruct,"sv_num");
  model->at_upper_bound = (long) restoreValue(mxStruct,"upper_bound");
  model->b = restoreValue(mxStruct,"b");
  model->totwords = (long) restoreValue(mxStruct, "totwords");
  model->totdoc = (long) restoreValue(mxStruct, "totdoc");
  model->loo_error = restoreValue(mxStruct, "loo_error");
  model->loo_recall = restoreValue(mxStruct, "loo_recall");
  model->xa_error = restoreValue(mxStruct, "xa_error");
  model->xa_recall = restoreValue(mxStruct, "xa_recall");
  model->xa_precision = restoreValue(mxStruct, "xa_precision");
  model->maxdiff = restoreValue(mxStruct, "maxdiff");
  model->r_delta_sq = restoreValue(mxStruct, "r_delta_sq");
  model->r_delta_avg = restoreValue(mxStruct, "r_delta_avg");
  model->model_length = restoreValue(mxStruct, "model_length");
  model->loss = restoreValue(mxStruct, "loss");
  model->vcdim = restoreValue(mxStruct, "vcdim");
  model->alpha = restoreArray(mxStruct, "alpha");
  model->lin_weights = restoreArray(mxStruct, "lin_weights");
  model->index = restoreArrayLong(mxStruct, "index");

  kstruct = mxGetField(mxStruct, 0, "kernel_parm"); 
  kparm = &(model->kernel_parm);
  kparm->kernel_type = (long) restoreValue(kstruct,"kernel_type"); 
  kparm->poly_degree = (long) restoreValue(kstruct,"poly_degree"); 
  kparm->rbf_gamma = (double) restoreValue(kstruct,"rbf_gamma"); 
  kparm->coef_lin = (double) restoreValue(kstruct,"coef_lin"); 
  kparm->coef_const = (double) restoreValue(kstruct,"coef_const"); 

  kstring = mxGetField(kstruct,0,"custom");
  buflen = (mxGetM(kstring) * mxGetN(kstring) * sizeof(mxChar)) + 2;
  buff = my_malloc(buflen);
  mxGetString(kstring,buff, buflen);

  strcpy((*kparm).custom, buff);
  my_free(buff);

  /*  if (model->supvec != NULL) { */
  svec = mxGetField(mxStruct,0,"supvec");
  if (svec == NULL) {
    printf("Warning: supvec field was not found.");
  }
  else {
    mexToDOC(svec, NULL, &model->supvec, NULL, NULL, NULL); 
  }

  return model;
}

#else
#include <stdio.h>
void _QWERTY( )
{
  char *x = "ANSI C requires something here. This is never really used.";
  printf("%s\n", x);
}

#endif
