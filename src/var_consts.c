#include "syspara.h"

void static_paras(FILE *fp1)
{
	int i,j,k;
	double r;
	double S_all,V_all;
	double d1,d2;
	double k1,k2;

	// Cell Geometry */

	//var.RGC = 2.0;  // The ratio of the actual surface area to the geometrical surface area (cm)
	d1 = var.length/(MEDIA_PATCH-2);
	d2 = var.cleft_width;

	//var.vcell = 1000*M_PI*var.a*var.a*var.length; // Cell Volume:3.801e-5 (uL)
	if(MEDIA_PATCH==5){
		for(i=0;i<MEDIA_PATCH;i++){
			if(i==0 || i==MEDIA_PATCH-1) {
				var.vcell[i] = 1000*M_PI*var.a*var.a*((var.length/(MEDIA_PATCH-2))*(1.0/4.0)); // (length*(1/(unit-2)*(1/4)))
			} else if (i==1 || i== MEDIA_PATCH -2 ){
				var.vcell[i] = 1000*M_PI*var.a*var.a*((var.length/(MEDIA_PATCH-2))*(3.0/4.0)); // (lenght*(1/(unit-2)*(3/4)))
			} else {
				var.vcell[i] = 1000*M_PI*var.a*var.a*var.length/(MEDIA_PATCH-2); // (lenght*(1/(unit-2)))
			}
			}
	} else if(MEDIA_PATCH==3){
		for(i=0;i<MEDIA_PATCH;i++){
			if(i==0 || i==MEDIA_PATCH-1) {
				var.vcell[i] = 1000*M_PI*var.a*var.a*var.length*(1.0/12.0); // (length*(1/(unit-2)*(1/4)))
			} else {
				var.vcell[i] = 1000*M_PI*var.a*var.a*var.length*(5.0/6.0); // (lenght*(1/(unit-2)))
			}
		}
	}
	// geometric membrane area: 7.671e-5 (cm^2)
	//var.ageo = M_PI*var.a*var.a+2.0*M_PI*var.a*var.length;

	for(i=0;i<MEDIA_PATCH;i++){
		if(i==0 || i == MEDIA_PATCH-1){
			var.ageo[i] = M_PI*var.a*var.a;
		} else {
			var.ageo[i] = 2.0*M_PI*var.a*var.length/(MEDIA_PATCH-2);
		}
	}
		
	for(i=0;i<MEDIA_PATCH;i++){
		var.acap[i] = var.ageo[i]*var.RGC;          // Capacitive membrane area: 1.534e-4 cm^2 (cm^2)
		var.vmyo[i] = var.vcell[i]*0.68;      // Myoplasm volume (uL) = 68% for Cell volume
		var.vss[i] = var.vcell[i]*0.02;     // Mitochondria volume (uL) = 26% for cell volume
		var.vsr[i] = var.vcell[i]*0.06;       // SR volume (uL)
		var.vnsr[i] = var.vcell[i]*0.0552;    // NSR volume (uL)
		var.vjsr[i] = var.vcell[i]*0.0048;    // JSR volume (uL)
		var.vcleft[i] = var.vcell[i]*0.12/0.88;  // Cleft volume (uL)
	}

	// surface area of junctional and non-Junctional membran unit   
	S_all = 0; V_all = 0;
	for(i=0;i<MEDIA_PATCH;i++){
		if (i==0 || i== MEDIA_PATCH-1){
			var.s2[i] = M_PI*var.a*var.a;
		} else {
			var.s2[i] = 2.0*M_PI*var.a*d1;
		}
		S_all+=var.s2[i];
		printf("s2[%d]=%e\n",i,var.s2[i]);
		fprintf(fp1,"s2[%d]=%e\n",i,var.s2[i]);
	}
	//printf("total_Surface=%e:%e, Cm=%e\n",S_all,2*M_PI*var.a*var.length+2*M_PI*var.a*var.a,var.RGC*S_all*1E-6);
	//fprintf(fp1,"total_Surface=%e:%e, Cm=%e\n",S_all,2*M_PI*var.a*var.length+2*M_PI*var.a*var.a,var.RGC*S_all*1E-6);
	printf("total_Surface=%e:%e, Cm,jm1=%e Cm,lm=%e, Cm,jm2=%e\n",S_all,2*M_PI*var.a*var.length+2*M_PI*var.a*var.a,var.acap[0]*1E-6,var.acap[1]*1E-6,var.acap[2]*1E-6);
	fprintf(fp1,"total_Surface=%e:%e, Cm,jm1=%e Cm,lm=%e, Cm,jm2=%e\n",S_all,2*M_PI*var.a*var.length+2*M_PI*var.a*var.a,var.acap[0]*1E-6,var.acap[1]*1E-6,var.acap[2]*1E-6);

	for(i=0;i<MEDIA_PATCH;i++){
		V_all+=var.vcell[i];
		printf("vell[%d]=%e\n",i,var.vcell[i]);
		fprintf(fp1,"vell[%d]=%e\n",i,var.vcell[i]);
	}
	printf("total_cell_volume=%e:%e\n",V_all,1000*M_PI*var.a*var.a*var.length);
	fprintf(fp1,"total_cell_Volume=%e:%e\n",V_all,1000*M_PI*var.a*var.a*var.length);

	// Temputure factor
	var.Tf = pow(Q,((T-273.0)-25.0)/10.0);
	fprintf(fp1,"Temputure factor = %lf\n",var.Tf);

	// Ion Valences
	var.zna = 1;  // Na valence
	var.zk = 1;   // K valence
	var.zca = 2;  // Ca valence
	
	for(i=0;i<MEDIA_PATCH;i++){
		var.vr1[i] = var.acap[i]/(F*var.vmyo[i]);
		var.vr2[i] = var.vss[i]/var.vmyo[i];
		var.vr3[i] = var.acap[i]/(F*var.vss[i]);
		var.vr4[i] = var.acap[i]/(2.0*F*var.vmyo[i]);
		var.vr5[i] = var.vnsr[i]/var.vmyo[i];
		var.vr6[i] = var.acap[i]/(2.0*F*var.vss[i]);
		var.vr7[i] = var.vjsr[i]/var.vss[i];
		var.vr8[i] = var.vjsr[i]/var.vnsr[i];
	}

	// Extracellular ion concentrations
		//var.nao = 132.0;     // Na (mM) Original value
		//var.nao = 140.0;     // Na (mM) Correct value
		//var.ko = 4.5;      // K (mM)
		//var.cao = 1.8;     // Ca (mM)

	// general
	var.RTonF = R*T/F;
	var.RTon2F = R*T/(var.zca*F);

	// Extracellular Concentrations
	var.nao = 140.0;     // Initial Bulk Medium Na (mM)
	var.ko = 5.4;      // Initial Bulk Medium K (mM)
	var.cao = 1.8;     // Initial Bulk Medium Ca (mM)

	// Na/K Permiability Ratio
	var.prnak = 0.01833;     // Na/K Permiability Ratio

	// Fast sodium current
	//ina.Gna_fast = 75.0;	// (mS/uF).
	ina.Ah_fast = 0.99;
	ina.Ah_slow = 0.01; 	// (1.0-ina.Ah_fast)
	ina.Ah_CaMK_fast = ina.Ah_fast;
	ina.Ah_CaMK_slow = ina.Ah_slow;
		
	// Late sodium current
	ina.tauhl = 200.0;				// (ms).
	ina.tauhl_CaMK = 3.0*ina.tauhl;	// 600 (ms).

	if(var.celltype==1){	// Endo
		ina.Gna_late = 0.0075;		// (mS/uF).
	} else if (var.celltype==2){// Mid
		ina.Gna_late = 0.0075;		
	} else if (var.celltype==3){// Epi
		ina.Gna_late = 0.0075*0.6;
	}

	//k1=ina.Gna_fast*2.0*M_PI*var.a*var.length; // assuming 50% of entire Na channel conductance of LM
	//k2=2.0*k1; // entire Na channel conductance 
	//printf("entire Na channel conductance = %e\n",k2);

	printf("na_switch=%d\n",var.na_switch);
	switch(var.na_switch){
		case 1: // uniform %
			for(i=0;i<MEDIA_SITE;i++){
				for(j=0;j<MEDIA_PATCH;j++){
					ina.gnaf_local[i][j] = ina.Gna_fast;
					ina.gnal_local[i][j] = ina.Gna_late;
				}
			}
			break;
		case 2: // localization version1
			ina.gnaf_all = ina.Gna_fast*(2.0*M_PI*var.a*var.a + 2.0*M_PI*var.a*var.length);
			ina.gnal_all = ina.Gna_late*(2.0*M_PI*var.a*var.a + 2.0*M_PI*var.a*var.length);
			for(i=0;i<MEDIA_SITE;i++){
				for(j=0;j<MEDIA_PATCH;j++){
					if (j==0||j==MEDIA_PATCH-1){
						ina.gnaf_local[i][j] = (0.5*var.jmrate1*ina.gnaf_all)/(M_PI*var.a*var.a);
						ina.gnal_local[i][j] = (0.5*var.jmrate1*ina.gnal_all)/(M_PI*var.a*var.a);
					} else {
						ina.gnaf_local[i][j] = (var.lmrate1*ina.gnaf_all/(MEDIA_PATCH-2))/(2.0*M_PI*var.a*var.length/(MEDIA_PATCH-2));
						ina.gnal_local[i][j] = (var.lmrate1*ina.gnal_all/(MEDIA_PATCH-2))/(2.0*M_PI*var.a*var.length/(MEDIA_PATCH-2));
					}
				}
			}
			break;
		case 3: // heterogeneous myofiber and local differences of NaCh in a myofiber case  %
			printf("3 was not defined\n");
			break;
		case 4: // control case  %
			for(i=0;i<MEDIA_SITE;i++){
				for(j=0;j<MEDIA_PATCH;j++){
					if (j==0||j==MEDIA_PATCH-1){
						ina.gnaf_local[i][j] = var.jmrate1*ina.Gna_fast*4.0;
						ina.gnal_local[i][j] = var.jmrate1*ina.Gna_late*4.0;
					} else {
						ina.gnaf_local[i][j] = var.lmrate1*ina.Gna_fast;
						ina.gnal_local[i][j] = var.lmrate1*ina.Gna_late;
					}
				}
			}
			break;
		case 5: // test case  %
			for(i=0;i<200;i++){
				for(j=0;j<MEDIA_PATCH;j++){
					if (j==0||j==MEDIA_PATCH-1){
						ina.gnaf_local[i][j] = var.jmrate1*ina.Gna_fast*4.0;
						ina.gnal_local[i][j] = var.jmrate1*ina.Gna_late*4.0;
					} else {
						ina.gnaf_local[i][j] = var.lmrate1*ina.Gna_fast;
						ina.gnal_local[i][j] = var.lmrate1*ina.Gna_late;
					}
				}
			}
			for(i=200;i<MEDIA_SITE;i++){
				for(j=0;j<MEDIA_PATCH;j++){
					if (j==0||j==MEDIA_PATCH-1){
						ina.gnaf_local[i][j] = var.jmrate2*ina.Gna_fast*4.0;
						ina.gnal_local[i][j] = var.jmrate2*ina.Gna_late*4.0;
					} else {
						ina.gnaf_local[i][j] = var.lmrate2*ina.Gna_fast;
						ina.gnal_local[i][j] = var.lmrate2*ina.Gna_late;
					}
				}
			}
			break;
		case 6: // test case  %
			for(i=0;i<200;i++){
				for(j=0;j<MEDIA_PATCH;j++){
					if (j==0||j==MEDIA_PATCH-1){
						ina.gnaf_local[i][j] = var.jmrate1*ina.Gna_fast*4.0;
						ina.gnal_local[i][j] = var.jmrate1*ina.Gna_late*4.0;
					} else {
						ina.gnaf_local[i][j] = var.lmrate1*ina.Gna_fast;
						ina.gnal_local[i][j] = var.lmrate1*ina.Gna_late;
					}
				}
			}
			for(i=200;i<300;i++){
				for(j=0;j<MEDIA_PATCH;j++){
					if (j==0||j==MEDIA_PATCH-1){
						ina.gnaf_local[i][j] = var.jmrate2*ina.Gna_fast*4.0;
						ina.gnal_local[i][j] = var.jmrate2*ina.Gna_late*4.0;
					} else {
						ina.gnaf_local[i][j] = var.lmrate2*ina.Gna_fast;
						ina.gnal_local[i][j] = var.lmrate2*ina.Gna_late;
					}
				}
			}
			for(i=300;i<500;i++){
				for(j=0;j<MEDIA_PATCH;j++){
					if (j==0||j==MEDIA_PATCH-1){
						ina.gnaf_local[i][j] = var.jmrate1*ina.Gna_fast*4.0;
						ina.gnal_local[i][j] = var.jmrate1*ina.Gna_late*4.0;
					} else {
						ina.gnaf_local[i][j] = var.lmrate1*ina.Gna_fast;
						ina.gnal_local[i][j] = var.lmrate1*ina.Gna_late;
					}
				}
			}
			break;
	}

	for(i=0;i<MEDIA_SITE;i++){
		for(k=0;k<MEDIA_PATCH;k++){
			printf("gnaf[%d][%d]=%lf gnal[%d][%d]=%lf\n",i,k,ina.gnaf_local[i][k],i,k,ina.gnal_local[i][k]);
			//fprintf(fp1,"gnaf[%d][%d]=%lf gnal[%d][%d]=%lf\n",i,k,ina.gnaf_local[i][k],i,k,ina.gnal_local[i][k]);
		}
	}
	for(k=0;k<MEDIA_PATCH;k++){
		fprintf(fp1,"gnaf[%d][%d]=%lf gnal[%d][%d]=%lf\n",MEDIA_SITE/2,k,ina.gnaf_local[MEDIA_SITE/2][k],MEDIA_SITE/2,k,ina.gnal_local[MEDIA_SITE/2][k]);
	}
	
	// Ito
	if(var.celltype==1){ //Endo
		ito.Gto = 0.02;
	} else if(var.celltype==2){ //Mid
		ito.Gto = 4.0*0.02;
	} else if(var.celltype==3){ //Epi
		ito.Gto = 4.0*0.02;
	} // (mS/uF)

	// L-type calcium current
	ical.Af_fast = 0.6;
	ical.Af_slow = 1.0 - ical.Af_fast;
	ical.taujca = 75.0; // (ms).
	ical.Af_CaMK_fast = ical.Af_fast;
	ical.Af_CaMK_slow = ical.Af_slow;

	ical.kmn = 0.002;
	ical.kp2n = 1000.0;

	ical.gacai = 1.0;         					// Activity coefficient of Ca
	ical.gacao = 0.341;     					// Activity coefficient of Ca
	ical.ganai = 0.75;      					// Activity coefficient of Na
	ical.ganao = 0.75;      					// Activity coefficient of Na
	ical.gaki = 0.75;       					// Activity coefficient of K
	ical.gako = 0.75;       					// Activity coefficient of K
	
	if(var.celltype == 1){ // Endo
		ical.pca = ical.rate*0.0001;      					// Permiability of membrane to Ca (cm/s)
		ical.pcana = 0.00125*ical.pca;   			// Permiability of membrane to Na (cm/s)
		ical.pcak = 3.574E-4*ical.pca;  			// Permiability of membrane to K (cm/s)
		ical.pca_CaMK = 1.1*ical.pca;				// Permiability of membrane to Ca (cm/s)
		ical.pcana_CaMK = 0.00125*ical.pca_CaMK;	// Permiability of membrane to Ca (cm/s)
		ical.pcak_CaMK = 3.574E-4*ical.pca_CaMK;	// Permiability of membrane to Ca (cm/s)
	} else if (var.celltype == 2){ // Mid
		ical.pca = ical.rate*0.0001*2.5;      					// Permiability of membrane to Ca (cm/s)
		ical.pcana = 0.00125*ical.pca;   			// Permiability of membrane to Na (cm/s)
		ical.pcak = 3.574E-4*ical.pca;  			// Permiability of membrane to K (cm/s)
		ical.pca_CaMK = 1.1*ical.pca;				// Permiability of membrane to Ca (cm/s)
		ical.pcana_CaMK = 0.00125*ical.pca_CaMK;	// Permiability of membrane to Ca (cm/s)
		ical.pcak_CaMK = 3.574E-4*ical.pca_CaMK;	// Permiability of membrane to Ca (cm/s)
	} else if (var.celltype == 3){ // Epi 
		ical.pca = ical.rate*0.0001*1.2;      					// Permiability of membrane to Ca (cm/s)
		ical.pcana = 0.00125*ical.pca;   			// Permiability of membrane to Na (cm/s)
		ical.pcak = 3.574E-4*ical.pca;  			// Permiability of membrane to K (cm/s)
		ical.pca_CaMK = 1.1*ical.pca;				// Permiability of membrane to Ca (cm/s)
		ical.pcana_CaMK = 0.00125*ical.pca_CaMK;	// Permiability of membrane to Ca (cm/s)
		ical.pcak_CaMK = 3.574E-4*ical.pca_CaMK;	// Permiability of membrane to Ca (cm/s)
	}

	// Rapid delayed rectifier potassium current (Ikr)
	ikr.rategkr = sqrt(var.ko/5.4);
	if(var.celltype == 1){ // Endo
		//ikr.Gkr = 0.046;
		ikr.Gkr = var.gkr;
	} else if(var.celltype == 2){ // Mid
		//ikr.Gkr = 0.8*0.046;
		ikr.Gkr = 0.8*var.gkr;
	} else if(var.celltype == 3){ // Epi
		//ikr.Gkr = 1.3*0.046;
		ikr.Gkr = 1.3*var.gkr;
	} // (mS/uF)

	// Slow delayed rectifier potassium current (Ikr)
	if(var.celltype == 1){	// Endo
		iks.Gks = 0.0034;
	} else if(var.celltype ==2 ){ // Mid
		iks.Gks = 0.0034;
	} else if(var.celltype == 3){ // Epi
		iks.Gks = 1.4*0.0034;
	} // (mS/uF)
	
	// Inward rectifier K current: Ik1
	ik1.rategk1 = sqrt(var.ko);
	if(var.celltype == 1){ // Endo
		ik1.Gk1 = 0.1908;
	} else if(var.celltype == 2){ // Mid
		ik1.Gk1 = 1.3*0.1908;
	} else if(var.celltype == 3){ // Epi
		ik1.Gk1 = 1.2*0.1908;
	} // (mS/uF)

	// Sodium-Calcium Exchanger V-S
	var.kna1 = 15.0;		// (mM)
	var.kna2 =  5.0;		// (mM)
	var.kna3 = 88.12;		// (mM)
	var.kasym = 12.5;			// (mM)
	var.omega_na = 6.0E+4;		// (Hz)
	var.omega_ca = 6.0E+4;		// (Hz)
	var.omega_naca = 5.0E+3;	// (Hz)
	var.kca_on = 1.5E+6;		// (mM/ms)
	var.kca_off = 5.0E+3;		// (Hz)
	var.qna = 0.5224;   // 
	var.qca = 0.1670;   // 
	var.km_ca_act = 150.0E-6;		// (mM)

	if(var.celltype == 1){ // Endo
		var.Gnaca = 0.0008;
	} else if(var.celltype == 2){ // Mid
		var.Gnaca = 1.4*0.0008;
	} else if(var.celltype == 3){ // Epi
		var.Gnaca = 1.1*0.0008;
	} // (uC/uF)

	// Sodium-Potassium Pump
	inak.kp1 = 949.5;		// (Hz)
	inak.km1 = 182.4;		// (1/mM)
	inak.kp2 = 687.2;		// (Hz)
	inak.km2 = 39.4;			// (Hz)
	inak.kp3 = 1899.0;		// (Hz)
	inak.km3 = 79300.0;		// (Hz/mM/mM)
	inak.kp4 = 639.0;		// (Hz)
	inak.km4 = 40.0;			// (Hz)
	inak.ko_nai = 9.073;		// (mM)
	inak.ko_nao = 27.78;		// (mM)
	inak.delta = -0.1550;
	inak.kki = 0.5;			// (mM)
	inak.kko = 0.3582;		// (mM)
	inak.MgADP = 0.05;
	inak.MgATP = 9.8;
	inak.k_MgATP = 1.698E-7;	// (mM)
	inak.H = 1.0E-7;		// (mM) proton
	inak.SigP = 4.2;			// (mM)
	inak.HP = 1.698E-7;	// (mM)
	inak.nap = 224.0;		// (mM)
	inak.kp = 292.0;		// (mM)
	inak.b1 = inak.km1*inak.MgADP;
	inak.a2 = inak.kp2;
	inak.a4 = (inak.kp4*inak.MgATP/inak.k_MgATP)/(1.0 + inak.MgATP/inak.k_MgATP);
	if(var.celltype == 1){ // Endo
		inak.G = 30.0;
	} else if(var.celltype == 2){ // Mid
		inak.G = 0.7*30.0;
	} else if(var.celltype == 3){ // Epi
		inak.G = 0.9*30.0;
	}

	// Sarcolemmal Ca Pump
	ipca.G = 0.0005;		// Max. Ca current through sarcolemmal Ca pump (mS/uF)
	ipca.km = 0.0005;		// Half-saturation concentration of sarcolemmal Ca pump (mM)

	// K Background Current (Max. conductance of K background)
	if(var.celltype == 1){ // Endo
		ikb.G = 0.003;
	} else if(var.celltype == 2){ // Mid
		ikb.G = 0.003;
	} else if(var.celltype == 3){ // Epi
		ikb.G = 0.6*0.003;
	} // (mS/uF)

	// Ca Background Current 
	icab.pcab = 2.5E-8;		// (cm/s)
	icab.gacai = 1.0;
	icab.gacao = 0.341;

	// Na Background Current 
	inab.pnab = 3.75E-10;    // (cm/s)

	// Ca/CaM dependent protein kinese (CaMK)
	CaMK.a = 0.05;		// (1.0/ms)
	CaMK.b = 0.00068;	// (1.0/ms)
	CaMK.z = 0.05;
	CaM.Km = 0.0015;	// (mM)
	CaMK.Km = 0.15;
	
	// diffusion fluxes
	jdiff.tau_na = 2.0;	// (ms)
	jdiff.tau_k = 2.0;	// (ms)
	jdiff.tau_ca = 0.2;	// (ms)
	
	// SR calcium release flux, via RyR (Jrel)
	jrel.b_tau = 4.75;					// (ms)
	jrel.a = 0.5*jrel.b_tau; 			// (ms)
	jrel.b_tau_CaMK = 1.25*jrel.b_tau;	// (ms)
	jrel.a_CaMK = 0.5*jrel.b_tau_CaMK;	// (ms)

	if(var.celltype == 1){ // Endo
		jrel.p = 1.0;
	} else if(var.celltype == 2){ // Mid
		jrel.p = 1.7;
	} else if(var.celltype == 3){ // Epi
		jrel.p = 1.0;
	}

	// calcium uptake via SERCA pump (Jup)
	jup.dKm_PLB = 0.00017;		// (mM)
	jup.dCaMK = 1.75;
	if(var.celltype == 1){ // Endo
		jup.p = 1.0;
	} else if(var.celltype == 2){ // Mid
		jup.p = 1.0;
	} else if(var.celltype == 3){ // Epi
		jup.p = 1.3;
	}

	// Translocation of Ca Ions from NSR to JSR
	jtr.tau = 100.0;      // Time constant of Ca transfer from NSR to JSR (ms)

	// Myoplasmic Ca Ion Concentration Changes 
	// Max. [Ca] buffered in CMDN (mM)
	if(var.celltype == 1){ // Endo
		var.cmdnbar = 0.050;
	} else if(var.celltype == 2){
		var.cmdnbar = 0.050;
	} else if(var.celltype == 3){
		var.cmdnbar = 1.3*0.050;
	}

	var.trpnbar = 0.070;   // Max. [Ca] buffered in TRPN (mM)
	var.bsrbar = 0.047;		// Max.[Ca] buffered in BSR (mM)
	var.bslbar = 1.124;		// Max.[Ca] buffered in BSL (mM)
	var.csqnbar = 10.0;     // Max. [Ca] buffered in CSQN (mM)
	var.kmcmdn = 0.00238;  // Equalibrium constant of buffering for CMDN (mM)
	var.kmtrpn = 0.0005;   // Equalibrium constant of buffering for TRPN (mM)
	var.kmbsr = 0.00087;  // Equalibrium constant of buffering for BSR (mM)
	var.kmbsl = 0.0087;   // Equalibrium constant of buffering for BSL (mM)
	var.kmcsqn = 0.8;     // Equalibrium constant of buffering for CSQN (mM)

	// Gap junction parameter
	// Gap junction conductance in Endo region = 2.5 [uS] -> S/gj = [Ohm cm^2]
	// resistance for Gap junction [Ohm cm^2]
	//var.Gj = 150.0;
	var.Gj = var.Rmyo;

	// surface area of Junctional membran unit or cross-section area of myocyte 
	var.s1 = M_PI*var.a*var.a;

	var.Ri = var.Rmyo*d1/var.s1;
	printf("Rmyo=%e, Ri=%e\n",var.Rmyo,var.Ri);
	fprintf(fp1,"Rmyo=%e, Ri=%e\n",var.Rmyo,var.Ri);
	var.gx0 = 2.0/var.Ri;

	var.Rg = var.Gj*var.length/var.s1;
	//var.Rg = var.Rmyo*var.length/var.s1;
	var.ggx = var.gj_rate*(1.0/var.Rg);
	printf("Rg=%e Gg=%e\n",var.Rg,var.ggx);   
	fprintf(fp1,"Rg=%e Gg=%e\n",var.Rg,var.ggx);   

	//var.Rj = var.cleft_R;
	var.Rj = 150.0/(8.0*M_PI*var.cleft_width);
	//var.Rj = var.Gj/(8.0*M_PI*var.cleft_width);
	var.cleft_R=var.Rj;
	//var.Rd = (150.0*(150.0/(8*M_PI*var.cleft_R)))/var.s1;
	var.Rd = (var.Gj*(var.Gj/(8*M_PI*var.cleft_R)))/var.s1;
	printf("Rd=%e, Rj=%e, cw=%e\n",var.Rd,var.Rj,var.cleft_width);
	fprintf(fp1,"Rd=%e, Rj=%e, cw=%e\n",var.Rd,var.Rj,var.cleft_width);

	switch(var.couple_type){
		case 0: // gap junction alone
			var.p1 = 0.0;
			var.p2 = 0.0;
			var.p3 = 0.0;
			break;
		case 1: // EF+gap junction
			var.p1 = var.gx0*(var.Rj+0.5*var.Rd);
			var.p2 = var.gx0*var.Rj;
			var.p3 = var.ggx*var.Rd*0.5;
			break;
	}

	printf("radius = %e\n",var.a);
	printf("unit length = %e\n",var.length);
	printf("unit surface area 2*Pi*a*dx = %e\n",2.0*M_PI*var.a*var.length);
	printf("Cross section area Pi*a*a = %e\n",M_PI*var.a*var.a);
	printf("intracellular resistance = %e\n",var.Rmyo*var.length/(M_PI*var.a*var.a));
	printf("Ri = %e, (Gi = %e)\n",var.Ri,1.0/var.Ri);
	printf("Rg = %e, (Gg = %e)\n",var.Rg,var.ggx);
	printf("gx0*(Rj+Rd/2)=%e, gx0*Rj=%e, ggx*Rd/2=%e\n",var.p1,var.p2,var.p3);
	printf("Ena=%f, Ek=%f\n",log(var.u[33][0][1]/var.nao)*var.RTonF,log(var.u[35][0][1]/var.ko)*var.RTonF);
	switch(var.couple_type){
		case 0: // gap junction alone
			printf("Gap junction mechanism alone model\n");
			break;
		case 1: // EF+gap junction
			printf("EF+Gap junction mechanism model\n");
		break;
	}

	fprintf(fp1,"radius = %e\n",var.a);
	fprintf(fp1,"unit length = %e\n",var.length);
	fprintf(fp1,"unit surface area 2*Pi*a*dx = %e\n",2.0*M_PI*var.a*var.length);
	fprintf(fp1,"Cross section area Pi*a*a = %e\n",M_PI*var.a*var.a);
	fprintf(fp1,"intracellular resistance = %e\n",var.Rmyo*var.length/(M_PI*var.a*var.a));
	fprintf(fp1,"Ri = %e, (Gi = %e)\n",var.Ri,1.0/var.Ri);
	fprintf(fp1,"Rg = %e, (Gg = %e)\n",var.Rg,var.ggx);
	fprintf(fp1,"gx0*(Rj+Rd/2)=%e, gx0*Rj=%e, ggx*Rd/2=%e\n",var.p1,var.p2,var.p3);
	fprintf(fp1,"Ena=%f, Ek=%f\n",log(var.u[33][0][1]/var.nao)*var.RTonF,log(var.u[35][0][1]/var.ko)*var.RTonF);
	switch(var.couple_type){
		case 0: // gap junction alone
			fprintf(fp1,"Gap junction mechanism alone model\n");
			break;
		case 1: // EF+gap junction
			fprintf(fp1,"EF+Gap junction mechanism model\n");
		break;
	}


}

