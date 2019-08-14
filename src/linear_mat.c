/* solving the matrix equation A*x=b using PARDISO */
 
#include "syspara.h"

void linear_coeff()
{
	int i,j,l; 
	int count,m,p,cflag;

	count=0;

	for(i=0;i<MEDIA_SITE;i++){ // first cell
		if(i==0){
			var.AT[0] = -1.0-var.p1-var.p3;		var.ja[0]=1;		var.ia[count]=1; count++;
			var.AT[1] = -var.p2+var.p3;			var.ja[1]=3*MEDIA_SITE;

			var.AT[2] = var.p1;					var.ja[2]=1;		var.ia[count]=3; count++;
			var.AT[3] = -1.0;					var.ja[3]=2;		
			var.AT[4] = var.p1;					var.ja[4]=3;
			var.AT[5] = var.p2;					var.ja[5]=4;
			var.AT[6] = var.p2;					var.ja[6]=3*MEDIA_SITE;

			var.AT[7] = -1.0-var.p1-var.p3;		var.ja[7]=3;		var.ia[count]=8; count++;
			var.AT[8] = -var.p2+var.p3;			var.ja[8]=4;

		} else if(i==MEDIA_SITE-1){ // last cell
			j = MEDIA_PATCH*(MEDIA_SITE-1);
			l = 9*i;
			var.AT[l+0] = -var.p2+var.p3;		var.ja[l+0]=j;		var.ia[count]=l+1; count++;
			var.AT[l+1] = -1.0-var.p1-var.p3;	var.ja[l+1]=j+1;

			var.AT[l+2] = var.p2;				var.ja[l+2]=1;		var.ia[count]=l+3; count++;
			var.AT[l+3] = var.p2;				var.ja[l+3]=j;		
			var.AT[l+4] = var.p1;				var.ja[l+4]=j+1;
			var.AT[l+5] = -1.0;					var.ja[l+5]=j+2;
			var.AT[l+6] = var.p1;				var.ja[l+6]=j+3;

			var.AT[l+7] = -var.p2;				var.ja[l+7]=1;		var.ia[count]=l+8; count++;
			var.AT[l+8] = -1.0-var.p1-var.p3;	var.ja[l+8]=j+3;	
																	var.ia[count]=l+10; // var.par_k+1;
		} else {
			j = MEDIA_PATCH*i;
			l = 9*i;
			var.AT[l+0] = -var.p2+var.p3;		var.ja[l+0]=j;		var.ia[count]=l+1; count++;
			var.AT[l+1] = -1.0-var.p1-var.p3;	var.ja[l+1]=j+1;

			var.AT[l+2] = var.p2;				var.ja[l+2]=j;		var.ia[count]=l+3; count++;
			var.AT[l+3] = var.p1;				var.ja[l+3]=j+1;
			var.AT[l+4] = -1.0;					var.ja[l+4]=j+2;
			var.AT[l+5] = var.p1;				var.ja[l+5]=j+3;
			var.AT[l+6] = var.p2;				var.ja[l+6]=j+4;

			var.AT[l+7] = -1.0-var.p1-var.p3;	var.ja[l+7]=j+3;	var.ia[count]=l+8; count++;
			var.AT[l+8] = -var.p2+var.p3;		var.ja[l+8]=j+4;
		}
	}

}


/* solving the matrix equation A*x=b using LAPACK */

void linear_vec()
{

	int i,j,l;

	for(i=0;i<MEDIA_SITE;i++){
		if(i==0){
			var.b[0] =  var.gx0*(var.u[0][0][0] - var.u[0][0][1]) - var.ggx*(var.u[0][MEDIA_SITE-1][2] - var.u[0][0][0]);
			var.b[1] = -var.gx0*(var.u[0][0][0] - 2.0*var.u[0][0][1] + var.u[0][0][2]);
			var.b[2] = -var.gx0*(var.u[0][0][1] - var.u[0][0][2]) + var.ggx*(var.u[0][0][2] - var.u[0][1][0]);
		} else if(i==MEDIA_SITE-1){
			l = 3*i;
			var.b[l+0] =  var.gx0*(var.u[0][i][0] - var.u[0][i][1]) -var.ggx*(var.u[0][i-1][2] -var.u[0][i][0]);
			var.b[l+1] = -var.gx0*(var.u[0][i][0] - 2.0*var.u[0][i][1] + var.u[0][i][2]);
			var.b[l+2] = -var.gx0*(var.u[0][i][1] - var.u[0][i][2]) + var.ggx*(var.u[0][i][2] - var.u[0][0][0]);
		} else {
			l = 3*i;
			var.b[l+0] =  var.gx0*(var.u[0][i][0] - var.u[0][i][1]) -var.ggx*(var.u[0][i-1][2] -var.u[0][i][0]);
			var.b[l+1] = -var.gx0*(var.u[0][i][0] - 2.0*var.u[0][i][1] + var.u[0][i][2]);
			var.b[l+2] = -var.gx0*(var.u[0][i][1] - var.u[0][i][2]) + var.ggx*(var.u[0][i][2] - var.u[0][i+1][0]);
		}
	}
	//for(i=0;i<MEDIA_SITE;i++) printf("b[%d]=%e\n",i,var.b[i]);

}

void linear_solve()
{

	long long n;
	MKL_INT i,j,k;

	n=MAT_SIZE;
	var.phase = 33;
	var.iparm[7] = 2; /* Max numbers of iterative refinement steps. */

	//PARDISO (pt, &maxfct, &mnum, &mtype, &phase, &n, a, ia, ja, &idum, &nrhs, iparm, &msglvl, b, x, &error);
	PARDISO_64 (var.pt, &var.maxfct, &var.mnum, &var.mtype, &var.phase, &n, var.AT, var.ia, var.ja, &var.idum, &var.nrhs, var.iparm, &var.msglvl, var.b, var.solv, &var.error);
		if (var.error != 0) {
			printf("\nERROR during solution: %ld", var.error);
			exit(3);
		}

/*	for (i = 0; i < n; i++) {
		printf("\n k=%d x [%d] = %16.15e",i,i,var.solv[i]);
	}
	printf ("\n");
	*/
//	exit(0);

}

void release_of_memory()
{
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
	long long n;
	n=MAT_SIZE;

	var.phase = -1; /* Release internal memory. */
	PARDISO_64 (var.pt, &var.maxfct, &var.mnum, &var.mtype, &var.phase, &n, &var.ddum, var.ia, var.ja, &var.idum, &var.nrhs, var.iparm, &var.msglvl, &var.ddum, &var.ddum, &var.error);
//	return 0;

}

