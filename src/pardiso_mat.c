/*
********************************************************************************
*                              INTEL CONFIDENTIAL
*   Copyright(C) 2004-2011 Intel Corporation. All Rights Reserved.
*   The source code contained  or  described herein and all documents related to
*   the source code ("Material") are owned by Intel Corporation or its suppliers
*   or licensors.  Title to the  Material remains with  Intel Corporation or its
*   suppliers and licensors. The Material contains trade secrets and proprietary
*   and  confidential  information of  Intel or its suppliers and licensors. The
*   Material  is  protected  by  worldwide  copyright  and trade secret laws and
*   treaty  provisions. No part of the Material may be used, copied, reproduced,
*   modified, published, uploaded, posted, transmitted, distributed or disclosed
*   in any way without Intel's prior express written permission.
*   No license  under any  patent, copyright, trade secret or other intellectual
*   property right is granted to or conferred upon you by disclosure or delivery
*   of the Materials,  either expressly, by implication, inducement, estoppel or
*   otherwise.  Any  license  under  such  intellectual property  rights must be
*   express and approved by Intel in writing.
*
********************************************************************************
*   Content : MKL PARDISO C example
*
********************************************************************************
*/
/* -------------------------------------------------------------------- */
/* Example program to show the use of the "PARDISO" routine */
/* on symmetric linear systems */
/* -------------------------------------------------------------------- */
/* This program can be downloaded from the following site: */
/* www.pardiso-project.org */
/* */
/* (C) Olaf Schenk, Department of Computer Science, */
/* University of Basel, Switzerland. */
/* Email: olaf.schenk@unibas.ch */
/* -------------------------------------------------------------------- */
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>

//#include "mkl_pardiso.h"
//#include "mkl_types.h"
#include "syspara.h"


MKL_INT init_pardiso( void ) {
	
	MKL_INT count,m,p,cflag;
	long long n;
	MKL_INT64 i,j,k;
	n = MAT_SIZE;
	printf("into init_pardiso function\n");
/* -------------------------------------------------------------------- */
/* .. Setup Pardiso control parameters. */
/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		var.iparm[i] = 0;
	}
	var.iparm[0] = 1; /* No solver default */
	var.iparm[1] = 2; /* Fill-in reordering from METIS */
	var.iparm[2] = 8; /* Numbers of processors, value of OMP_NUM_THREADS */
	var.iparm[3] = 0; /* No iterative-direct algorithm */
	var.iparm[4] = 0; /* No user fill-in reducing permutation */
	var.iparm[5] = 0; /* Write solution into x */
	var.iparm[6] = 0; /* Not in use */
	var.iparm[7] = 2; /* Max numbers of iterative refinement steps */
	var.iparm[8] = 0; /* Not in use */
	var.iparm[9] = 13; /* Perturb the pivot elements with 1E-13 */
	var.iparm[10] = 1; /* Use nonsymmetric permutation and scaling MPS */
	var.iparm[11] = 0; /* Not in use */
	var.iparm[12] = 1; /* Maximum weighted matching algorithm is switched-on (default for non-symmetric) */
	var.iparm[13] = 0; /* Output: Number of perturbed pivots */
	var.iparm[14] = 0; /* Not in use */
	var.iparm[15] = 0; /* Not in use */
	var.iparm[16] = 0; /* Not in use */
	var.iparm[17] = -1; /* Output: Number of nonzeros in the factor LU */
	var.iparm[18] = -1; /* Output: Mflops for LU factorization */
	var.iparm[19] = 0; /* Output: Numbers of CG Iterations */
	var.maxfct = 1; /* Maximum number of numerical factorizations. */
	var.mnum = 1; /* Which factorization to use. */
	var.msglvl = 0; /* Print statistical information in file */
	var.error = 0; /* Initialize error flag */

/* -------------------------------------------------------------------- */
/* .. Initialize the internal solver memory pointer. This is only */
/* necessary for the FIRST call of the PARDISO solver. */
/* -------------------------------------------------------------------- */
	for (i = 0; i < 64; i++) {
		var.pt[i] = 0;
	}

/* -------------------------------------------------------------------- */
/* .. Reordering and Symbolic Factorization. This step also allocates */
/* all memory that is necessary for the factorization. */
/* -------------------------------------------------------------------- */
	var.phase = 11;
	//PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	PARDISO_64 (var.pt, &var.maxfct, &var.mnum, &var.mtype, &var.phase, &n, var.AT, var.ia, var.ja, &var.idum, &var.nrhs, var.iparm, &var.msglvl, &var.ddum, &var.ddum, &var.error);
	
	if (var.error != 0) {
		printf("\nERROR during symbolic factorization: %ld", var.error);
		exit(1);
	}
	//printf("\nReordering completed ... ");
	//printf("\nNumber of nonzeros in factors = %ld", iparm[17]);
	//printf("\nNumber of factorization MFLOPS = %ld", iparm[18]);

/* -------------------------------------------------------------------- */
/* .. Numerical factorization. */
/* -------------------------------------------------------------------- */
	var.phase = 22;
	//PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, &ddum, &ddum, &error);
	PARDISO_64 (var.pt, &var.maxfct, &var.mnum, &var.mtype, &var.phase, &n, var.AT, var.ia, var.ja, &var.idum, &var.nrhs, var.iparm, &var.msglvl, &var.ddum, &var.ddum, &var.error);
	
	if (var.error != 0) {
		printf("\nERROR during numerical factorization: %ld", var.error);
		exit(2);
	}

	return 0;
}
