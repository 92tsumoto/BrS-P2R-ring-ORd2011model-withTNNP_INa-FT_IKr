#include <string.h>
#include "syspara.h"

void vm_data(FILE *fp7, double time)
{

	int i,j;

	fprintf(fp7,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp7,"%lf ",var.u[0][i][j]);
		}
	}
	fprintf(fp7,"\n");

} 

void vm1_data(FILE *fp18, double time)
{

	int i,j;

	fprintf(fp18,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp18,"%lf ",var.u[0][i][j]);
		}
	}
	fprintf(fp18,"\n");

} 

void vm2_data(FILE *fp19, double time)
{

	int i,j;

	fprintf(fp19,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp19,"%lf ",var.u[0][i][j]);
		}
	}
	fprintf(fp19,"\n");

} 

void ina_data(FILE *fp8, double time)
{

	int i,j;
	
	fprintf(fp8,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp8,"%lf ",ina.total[i][j]);
		}
	}
	fprintf(fp8,"\n");
	
}

void ical_data(FILE *fp9, double time)
{

	int i,j;
	
	fprintf(fp9,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp9,"%lf ",ical.ica[i][j]);
		}
	}
	fprintf(fp9,"\n");
	 	
}

void ik1_data(FILE *fp10, double time)
{

	int i,j;

	fprintf(fp10,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp10,"%lf ",ik1.ik[i][j]);
		}
	}
	fprintf(fp10,"\n");
 	
}

void iks_data(FILE *fp11, double time)
{

	int i,j;

	fprintf(fp11,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp11,"%lf ",iks.ik[i][j]);
		}
	}
	fprintf(fp11,"\n");

}

void ikr_data(FILE *fp12, double time)
{

	int i,j;

	fprintf(fp12,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp12,"%lf ",ikr.ik[i][j]);
		}
	}
	fprintf(fp12,"\n");

}

void ito_data(FILE *fp13, double time)
{

	int i,j;

	fprintf(fp13,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp13,"%lf ",ito.ik[i][j]);
		}
	}
	fprintf(fp13,"\n");

}

void inak_data(FILE *fp14, double time)
{

	int i,j;

	fprintf(fp14,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp14,"%lf ",inak.inak[i][j]);
		}
	}
	fprintf(fp14,"\n");

}

void incx_data(FILE *fp15, double time)
{

	int i,j;

	fprintf(fp15,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp15,"%lf ",var.inaca[i][j]);
		}
	}
	fprintf(fp15,"\n");

}

void cleft(FILE *fp16, double time)
{

	int i;

	fprintf(fp16,"%lf ",time);
	for(i=0;i<MEDIA_SITE-1;i++){
		fprintf(fp16,"%lf ",var.cleft_potential[i]);
	}
	fprintf(fp16,"\n");

}

void intra_v_data(FILE *fp17, double time)
{

	int i,j;

	fprintf(fp17,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp17,"%lf ",var.inner_v[i][j]);
		}
	}
	fprintf(fp17,"\n");

}

void nai_data(FILE *fp18, double time)
{

	int i,j;

	fprintf(fp18,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp18,"%e ",var.u[32][i][j]);
		}
	}
	fprintf(fp18,"\n");

}

void ki_data(FILE *fp19, double time)
{

	int i,j;

	fprintf(fp19,"%lf ",time);
	for (i=0;i<MEDIA_SITE;i++){
		for (j=0;j<MEDIA_PATCH;j++){
			fprintf(fp19,"%e ",var.u[34][i][j]);
		}
	}
	fprintf(fp19,"\n");

}

/*
void cleft_INa_data(FILE *fp14)
{

	int i;

	for(i=0;i<MEDIA_SITE-1;i++){
		fprintf(fp14,"%d %d %e \n",var.beat,i+1,var.max_ina_id[i]);
	}

}
*/

