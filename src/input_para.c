#include <string.h>
#include "syspara.h"

void input_para(FILE *fpin)
{
	int i,ii;

	fscanf(fpin,"%d",&var.couple_type);
	fscanf(fpin,"%d",&var.model_type);
	fscanf(fpin,"%d",&var.drug_type);
	fscanf(fpin,"%d",&var.drug_mode);
	fscanf(fpin,"%d",&var.na_switch);
	fscanf(fpin,"%lf",&ina.Gna_fast);
	fscanf(fpin,"%lf",&var.jmrate1);
	fscanf(fpin,"%lf",&var.lmrate1);
	fscanf(fpin,"%lf",&var.jmrate2);
	fscanf(fpin,"%lf",&var.lmrate2);
	fscanf(fpin,"%lf",&var.a);
	fscanf(fpin,"%lf",&var.length);
	fscanf(fpin,"%lf",&var.Rmyo);
	fscanf(fpin,"%lf",&var.gj_rate);
	fscanf(fpin,"%lf",&var.cleft_R);
	fscanf(fpin,"%lf",&var.cleft_width);
	fscanf(fpin,"%lf",&var.BCL);
	fscanf(fpin,"%lf",&var.Istim_base);
	fscanf(fpin,"%d",&var.celltype);
	fscanf(fpin,"%lf",&var.gkr);
	fscanf(fpin,"%lf",&var.block_rate);
	fscanf(fpin,"%lf",&var.drug_conc);
	fscanf(fpin,"%lf",&var.ic50);
	fscanf(fpin,"%lf",&var.hill);
	fscanf(fpin,"%lf",&var.ec50);
	fscanf(fpin,"%lf",&var.ehill);
	fscanf(fpin,"%d",&var.datas);
	for (ii = 0; ii < var.datas; ii++){
		for (i=0;i<NN;i++){
			fscanf(fpin,"%lf",&var.x0[ii][i]);
		}
		fscanf(fpin, "%lf", &var.tsign[ii]);
		fscanf(fpin, "%lf", &var.tend[ii]);
	}
	fscanf(fpin,"%lf",&var.RGC);
	fscanf(fpin,"%lf",&ical.rate);
	fscanf(fpin,"%d",&var.stim_sw);
	fscanf(fpin,"%d",&var.l);
	fscanf(fpin,"%d",&var.m);
	fscanf(fpin,"%d",&var.write);
	fscanf(fpin,"%d",&var.out_data);
	fscanf(fpin,"%d",&var.out_data_plus);
	fscanf(fpin,"%d",&var.out_data_plus2);
	printf("RGC=%lf\n",var.RGC);
	printf("Pca_rate=%lf\n",ical.rate);
	printf("stim_sw=%d\n",var.stim_sw);
	printf("l=%d\n",var.l);
	printf("m=%d\n",var.m);
	printf("write=%d\n",var.write);
	printf("out data flag=%d\n",var.out_data);
	printf("out plus data flag=%d\n",var.out_data_plus);
	printf("out plus data flag=%d\n",var.out_data_plus2);

}
