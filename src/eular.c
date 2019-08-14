
#include "syspara.h"

void eular(int n, double h, double t)
{
    
	int i,j,k;  

	// Calculation of input current into each unit.
	linear_vec();
	linear_solve();

	// Store of cleft potential
		for(i=0;i<MEDIA_SITE-1;i++){
			var.cleft_potential[i] = var.Rj*(var.solv[(MEDIA_PATCH-1)+MEDIA_PATCH*i] + var.solv[MEDIA_PATCH+MEDIA_PATCH*i]);
			var.cleft_axial[i][0] = 0.5*var.Rd*var.solv[(MEDIA_PATCH-1)+MEDIA_PATCH*i];
			var.cleft_axial[i][1] = 0.5*var.Rd*var.solv[MEDIA_PATCH+MEDIA_PATCH*i];
		}
	// Determination of membrane potential in JM from intracellular - extracellular potential
		for(i=0;i<MEDIA_SITE;i++){
			for(j=0;j<MEDIA_PATCH;j++){
				if(j==0 && i>0){
					var.inner_v[i][j] = var.u[0][i][j] + 0.5*var.Rd*var.solv[MEDIA_PATCH*i] + var.Rj*(var.solv[-1+MEDIA_PATCH*i]+var.solv[MEDIA_PATCH*i]);
				} else if(j==MEDIA_PATCH-1 && i<MEDIA_SITE-1){
					var.inner_v[i][j] = var.u[0][i][j] + 0.5*var.Rd*var.solv[j+MEDIA_PATCH*i] + var.Rj*(var.solv[j+MEDIA_PATCH*i]+var.solv[MEDIA_PATCH*(i+1)]);
				} else{
					var.inner_v[i][j] = var.u[0][i][j];
				}
			}
		}

	// Calculation of Ion channel currents
	for(i=0;i<MEDIA_SITE;i++){
		for(j=0;j<MEDIA_PATCH;j++){
			function(i,j,t);
			for(k=0;k<n;k++){
				if(k == 0){
					//var.inner_current[i][j] = 1000*h/(var.RGC*var.s2[j])*var.solv[j+MEDIA_PATCH*i];
					var.inner_current[i][j] = 1000*h/(var.acap[j])*var.solv[j+MEDIA_PATCH*i];
					var.u[k][i][j] = var.u[k][i][j] + h*var.k1[k][i][j] + var.inner_current[i][j];
					//printf("ion_current[%d][%d]=%e, in_current[%d][%d]=%e\n",i,j,xtemp[0][i][j],i,j,var.inner_current[i][j]);
				} else {
					var.u[k][i][j] = var.u[k][i][j] + h*var.k1[k][i][j];
				}
			}
		}
	}

// Determination of membrane potential in JM from intracellular - extracellular potential
/*		for(i=0;i<MEDIA_SITE-1;i++){
			for(j=0;j<MEDIA_PATCH;j++){
				if(j==0 && i>0){
					var.inner_v[i][j] = var.u[0][i][j] + 0.5*var.Rd[i]*var.solv[MEDIA_PATCH*i] + var.Rj[i]*(var.solv[-1+MEDIA_PATCH*i]+var.solv[MEDIA_PATCH*i]);
				} else if(j==MEDIA_PATCH-1 && i<MEDIA_SITE-1){
					var.inner_v[i][j] = var.u[0][i][j] + 0.5*var.Rd[i]*var.solv[j+MEDIA_PATCH*i] + var.Rj[i]*(var.solv[j+MEDIA_PATCH*i]+var.solv[MEDIA_PATCH*(i+1)]);
				} else{
					var.inner_v[i][j] = var.u[0][i][j];
				}
			}
		}
*/
}
