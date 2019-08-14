#include "syspara.h"

void function(int site, int patch, double t)
{

	comp_reversal_potential(site,patch);
	comp_CaMK(site,patch);
	comp_ina(site,patch);
	comp_ito(site,patch);
	comp_ical(site,patch);
	comp_ikr(site,patch);
	comp_iks(site,patch);
	comp_ik1(site,patch);
	comp_inaca(site,patch);
	comp_inak(site,patch);
	comp_ipca(site,patch);
	comp_ikb(site,patch);
	comp_icab(site,patch);
	comp_inab(site,patch);
	comp_diff(site,patch);
	comp_jrel(site,patch);
	comp_jup(site,patch);
	comp_jtr(site,patch);
	comp_concentration(site,patch);
	
	var.Ina_i_total[site][patch] = ina.total[site][patch] + inab.na[site][patch] + 3.0*inak.inak[site][patch] + 3.0*var.inaca_i[site][patch];
	var.Ina_ss_total[site][patch] = ical.icana[site][patch] + 3.0*var.inaca_ss[site][patch];
	var.Ik_ss_total[site][patch] = ical.icak[site][patch];
	var.Ica_i_total[site][patch] = ipca.ca[site][patch] + icab.ca[site][patch] - 2.0*var.inaca_i[site][patch];
	var.Ica_ss_total[site][patch] = ical.ica[site][patch] - 2.0*var.inaca_ss[site][patch];
	
	if(var.stim_sw==0){
		//if (site == 0 && patch == 1){
		if (site == 299 && patch == 1){
			var.Ik_i_total[site][patch] = ito.ik[site][patch] + ikr.ik[site][patch] + iks.ik[site][patch] + ik1.ik[site][patch] + ikb.k[site][patch] - 2.0*inak.inak[site][patch] + var.Istim;
			var.Itotal[site][patch] = ina.total[site][patch] + ito.ik[site][patch] + ical.ica[site][patch] + ical.icana[site][patch] 
							+ ical.icak[site][patch] + ikr.ik[site][patch] + iks.ik[site][patch] + ik1.ik[site][patch] 
							+ var.inaca[site][patch] + inak.inak[site][patch] + inab.na[site][patch] + icab.ca[site][patch] + ikb.k[site][patch] + ipca.ca[site][patch] + var.Istim;
		} else {
			var.Ik_i_total[site][patch] = ito.ik[site][patch] + ikr.ik[site][patch] + iks.ik[site][patch] + ik1.ik[site][patch] + ikb.k[site][patch] - 2.0*inak.inak[site][patch];
			var.Itotal[site][patch] = ina.total[site][patch] + ito.ik[site][patch] + ical.ica[site][patch] + ical.icana[site][patch] 
							+ ical.icak[site][patch] + ikr.ik[site][patch] + iks.ik[site][patch] + ik1.ik[site][patch] 
							+ var.inaca[site][patch] + inak.inak[site][patch] + inab.na[site][patch] + icab.ca[site][patch] + ikb.k[site][patch] + ipca.ca[site][patch];
		}
	} else {
		if (site == MEDIA_SITE-1 && patch == 1){
			var.Ik_i_total[site][patch] = ito.ik[site][patch] + ikr.ik[site][patch] + iks.ik[site][patch] + ik1.ik[site][patch] + ikb.k[site][patch] - 2.0*inak.inak[site][patch] + var.Istim;
			var.Itotal[site][patch] = ina.total[site][patch] + ito.ik[site][patch] + ical.ica[site][patch] + ical.icana[site][patch] 
							+ ical.icak[site][patch] + ikr.ik[site][patch] + iks.ik[site][patch] + ik1.ik[site][patch] 
							+ var.inaca[site][patch] + inak.inak[site][patch] + inab.na[site][patch] + icab.ca[site][patch] + ikb.k[site][patch] + ipca.ca[site][patch] + var.Istim;
		} else {
			var.Ik_i_total[site][patch] = ito.ik[site][patch] + ikr.ik[site][patch] + iks.ik[site][patch] + ik1.ik[site][patch] + ikb.k[site][patch] - 2.0*inak.inak[site][patch];
			var.Itotal[site][patch] = ina.total[site][patch] + ito.ik[site][patch] + ical.ica[site][patch] + ical.icana[site][patch] 
							+ ical.icak[site][patch] + ikr.ik[site][patch] + iks.ik[site][patch] + ik1.ik[site][patch] 
							+ var.inaca[site][patch] + inak.inak[site][patch] + inab.na[site][patch] + icab.ca[site][patch] + ikb.k[site][patch] + ipca.ca[site][patch];
		}
	}

	// membrane potential
	var.k1[0][site][patch] = -var.Itotal[site][patch];
	//Fast sodium current
	var.k1[1][site][patch] = (ina.mss - var.u[1][site][patch])/ina.taum; // m
	var.k1[2][site][patch] = (ina.hss - var.u[2][site][patch])/ina.tauh_fast; // h_fast
	var.k1[3][site][patch] = (ina.hss - var.u[3][site][patch])/ina.tauh_slow; // h_slow
	var.k1[4][site][patch] = (ina.hss - var.u[4][site][patch])/ina.tauj; // j
	var.k1[5][site][patch] = (ina.hCaMKss - var.u[5][site][patch])/ina.tauh_CaMK_slow; // h_CaMK_slow
	var.k1[6][site][patch] = (ina.hss - var.u[6][site][patch])/ina.tauj_CaMK; // j_CaMK
	//late sodium current
	var.k1[7][site][patch] = (ina.mlss - var.u[7][site][patch])/ina.tauml; // ml
	var.k1[8][site][patch] = (ina.hlss - var.u[8][site][patch])/ina.tauhl; // hl
	var.k1[9][site][patch] = (ina.hlCaMKss - var.u[9][site][patch])/ina.tauhl_CaMK; // hl
	//Transient outward current
	var.k1[10][site][patch] = (ito.ass - var.u[10][site][patch])/ito.taua;
	var.k1[11][site][patch] = (ito.iss - var.u[11][site][patch])/ito.taui_fast;
	var.k1[12][site][patch] = (ito.iss - var.u[12][site][patch])/ito.taui_slow;
	var.k1[13][site][patch] = (ito.aCaMKss - var.u[13][site][patch])/ito.taua;
	var.k1[14][site][patch] = (ito.iss - var.u[14][site][patch])/ito.taui_CaMK_fast;
	var.k1[15][site][patch] = (ito.iss - var.u[15][site][patch])/ito.taui_CaMK_slow;
	// LTCC
	var.k1[16][site][patch] = (ical.dss - var.u[16][site][patch])/ical.taud;
	var.k1[17][site][patch] = (ical.fss - var.u[17][site][patch])/ical.tauf_fast;
	var.k1[18][site][patch] = (ical.fss - var.u[18][site][patch])/ical.tauf_slow;
	var.k1[19][site][patch] = (ical.fss - var.u[19][site][patch])/ical.taufca_fast;
	var.k1[20][site][patch] = (ical.fss - var.u[20][site][patch])/ical.taufca_slow;
	var.k1[21][site][patch] = (ical.fss - var.u[21][site][patch])/ical.taujca;
	var.k1[22][site][patch] = (ical.fss - var.u[22][site][patch])/ical.tauf_CaMK_fast;
	var.k1[23][site][patch] = (ical.fss - var.u[23][site][patch])/ical.taufca_CaMK_fast;
	var.k1[24][site][patch] = ical.alpha_n*ical.kp2n - var.u[24][site][patch]*ical.km2n;
	// Ikr
	var.k1[25][site][patch] = (ikr.xrss - var.u[25][site][patch])/ikr.tauxr_fast;
	var.k1[26][site][patch] = (ikr.xrss - var.u[26][site][patch])/ikr.tauxr_slow;
	// Iks
	var.k1[27][site][patch] = (iks.xs1ss - var.u[27][site][patch])/iks.tauxs1;
	var.k1[28][site][patch] = (iks.xs1ss - var.u[28][site][patch])/iks.tauxs2;
	// Ik1
	var.k1[29][site][patch] = (ik1.k1ss - var.u[29][site][patch])/ik1.tauk1;
	// CaMK
	var.k1[30][site][patch] = CaMK.a*CaMK.bound[site][patch]*(CaMK.bound[site][patch]+var.u[30][site][patch]) - CaMK.b*var.u[30][site][patch];
	// Jrel
	var.k1[31][site][patch] = (jrel.NPss[site][patch] - var.u[31][site][patch])/jrel.tau_NP[site][patch]; 
	var.k1[32][site][patch] = (jrel.CaMKss[site][patch] - var.u[32][site][patch])/jrel.tau_CaMK[site][patch]; 
	// [Na]i
	var.k1[33][site][patch] = -var.Ina_i_total[site][patch]*var.vr1[patch] + jdiff.na[site][patch]*var.vr2[patch];
	// [Na]ss
	var.k1[34][site][patch] = -var.Ina_ss_total[site][patch]*var.vr3[patch] - jdiff.na[site][patch];
	// [K]i
	var.k1[35][site][patch] = -var.Ik_i_total[site][patch]*var.vr1[patch] + jdiff.k[site][patch]*var.vr2[patch];
	// [K]ss
	var.k1[36][site][patch] = -var.Ik_ss_total[site][patch]*var.vr3[patch] - jdiff.k[site][patch];
	// [Ca]i
	var.k1[37][site][patch] = var.b_Ca_i[site][patch]*(-var.Ica_i_total[site][patch]*var.vr4[patch] - jup.ca[site][patch]*var.vr5[patch] + jdiff.ca[site][patch]*var.vr2[patch]);
	// [Ca]ss
	var.k1[38][site][patch] = var.b_Ca_ss[site][patch]*(-var.Ica_ss_total[site][patch]*var.vr6[patch] +jrel.ca[site][patch]*var.vr7[patch] - jdiff.ca[site][patch]);
	// [Ca]nsr
	var.k1[39][site][patch] = jup.ca[site][patch] - jtr.ca[site][patch]*var.vr8[patch];
	// [Ca]jsr
	var.k1[40][site][patch] = var.b_Ca_jsr[site][patch]*(jtr.ca[site][patch]-jrel.ca[site][patch]);
	// Ikr2 (facilitated fraction)
	var.k1[41][site][patch] = var.Tf*(ikr.xrss2 - var.u[41][site][patch])/ikr.tauxr_fast2;
	var.k1[42][site][patch] = var.Tf*(ikr.xrss2 - var.u[42][site][patch])/ikr.tauxr_slow2;


	//printf("NPss=%lf,NPp=%lf\n",jrel.NPss,jrel.ca[site][patch]MKss);
	//for(i=0;i<NN;i++){
	//	printf("var.u[%d]=%e\n",i,var.k1[i]);
	//}
}

void comp_ina(int site, int patch)
{
	//MKL_INT iV=0;
	int iV=0;
	double V1,V2,d1,d2;
	
	V1 = (var.u[0][site][patch]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;
	//printf("vm[%d][%d]=%lf, iV=%d,V1=%f,V2=%f,d1=%f,d2=%f\n",site,patch,var.u[0][site][patch],iV,V1,V2,d1,d2);

	if(var.model_type==1){ // original ORd
		ina.mss = ina.Tmss[iV]*d2 + ina.Tmss[iV+1]*d1;
		ina.taum = ina.Ttaum[iV]*d2 + ina.Ttaum[iV+1]*d1;
		ina.hss = ina.Thss[iV]*d2 + ina.Thss[iV+1]*d1;
		ina.tauh_fast = ina.Ttauh_fast[iV]*d2 + ina.Ttauh_fast[iV+1]*d1;
		ina.tauh_slow = ina.Ttauh_slow[iV]*d2 + ina.Ttauh_slow[iV+1]*d1;
		ina.h = ina.Ah_fast*var.u[2][site][patch] + ina.Ah_slow*var.u[3][site][patch];
		ina.jss = ina.Tjss[iV]*d2 + ina.Tjss[iV+1]*d1;
		ina.tauj = ina.Ttauj[iV]*d2 + ina.Ttauj[iV+1]*d1;
		ina.hCaMKss = ina.ThCaMKss[iV]*d2 + ina.ThCaMKss[iV+1]*d1;
	} else if(var.model_type==2){ // Dutta et al modification
		ina.mss = ina.Tmss[iV]*d2 + ina.Tmss[iV+1]*d1;
		ina.taum = ina.Ttaum[iV]*d2 + ina.Ttaum[iV+1]*d1;
		ina.hss = ina.Thss[iV]*d2 + ina.Thss[iV+1]*d1;
		ina.tauh_fast = ina.Ttauh_fast[iV]*d2 + ina.Ttauh_fast[iV+1]*d1;
		ina.tauh_slow = ina.Ttauh_slow[iV]*d2 + ina.Ttauh_slow[iV+1]*d1;
		ina.h = ina.Ah_fast*var.u[2][site][patch] + ina.Ah_slow*var.u[3][site][patch];
		ina.jss = ina.Tjss[iV]*d2 + ina.Tjss[iV+1]*d1;
		ina.tauj = ina.Ttauj[iV]*d2 + ina.Ttauj[iV+1]*d1;
		ina.hCaMKss = ina.ThCaMKss[iV]*d2 + ina.ThCaMKss[iV+1]*d1;
	} else if(var.model_type==3){ // TP INa substitution	
		ina.mss = ina.Tmss[iV]*d2 + ina.Tmss[iV+1]*d1;
		ina.taum = ina.Ttaum[iV]*d2 + ina.Ttaum[iV+1]*d1;
		ina.hss = ina.Thss[iV]*d2 + ina.Thss[iV+1]*d1;
		ina.tauh_fast = ina.Ttauh_fast[iV]*d2 + ina.Ttauh_fast[iV+1]*d1;
		ina.tauh_slow = ina.Ttauh_slow[iV]*d2 + ina.Ttauh_slow[iV+1]*d1;
		//ina.h = ina.Ah_fast*var.u[2][site][patch] + ina.Ah_slow*var.u[3][site][patch];
		ina.h = var.u[2][site][patch];
		ina.jss = ina.Tjss[iV]*d2 + ina.Tjss[iV+1]*d1;
		ina.tauj = ina.Ttauj[iV]*d2 + ina.Ttauj[iV+1]*d1;
		ina.hCaMKss = ina.ThCaMKss[iV]*d2 + ina.ThCaMKss[iV+1]*d1;
	}
	ina.tauh_CaMK_slow = 3.0*ina.tauh_slow;
	//ina.hCaMK = ina.Ah_CaMK_fast*var.u[2] + ina.Ah_CaMK_slow*var.u[5];	// h_CaMK_fast = h_fast (var.u[2])
	ina.hCaMK = ina.Ah_fast*var.u[2][site][patch] + ina.Ah_slow*var.u[5][site][patch];	// h_CaMK_fast = h_fast (var.u[2])

	//ina.jCaMKss = ina.jss; //jss=hss;
	ina.tauj_CaMK = 1.46*ina.tauj;
	
	ina.fast_pCaMK[site][patch] = 1.0/(1.0 + CaMK.Km/CaMK.active[site][patch]);
	ina.fast[site][patch] = ina.gnaf_local[site][patch]*(var.u[0][site][patch]-var.Ena[site][patch])
							*var.u[1][site][patch]*var.u[1][site][patch]*var.u[1][site][patch]*((1.0-ina.fast_pCaMK[site][patch])*ina.h*var.u[4][site][patch]
							+ina.fast_pCaMK[site][patch]*ina.hCaMK*var.u[6][site][patch]);

	ina.mlss = ina.Tmlss[iV]*d2 + ina.Tmlss[iV+1]*d1;
	ina.tauml = ina.Ttauml[iV]*d2 + ina.Ttauml[iV+1]*d1;
	ina.hlss = ina.Thlss[iV]*d2 + ina.Thlss[iV+1]*d1;
	ina.hlCaMKss = ina.ThlCaMKss[iV]*d2 + ina.ThlCaMKss[iV+1]*d1;

	ina.late_pCaMK[site][patch] = 1.0/(1.0 + CaMK.Km/CaMK.active[site][patch]);
	ina.late[site][patch] = ina.gnal_local[site][patch]*(var.u[0][site][patch]-var.Ena[site][patch])
							*var.u[7][site][patch]*((1.0-ina.late_pCaMK[site][patch])*var.u[8][site][patch] + ina.late_pCaMK[site][patch]*var.u[9][site][patch]);
	
	ina.total[site][patch] = ina.fast[site][patch] + ina.late[site][patch];
}

// Ito Transient Outward Current
void comp_ito (int site, int patch)
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (var.u[0][site][patch]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ito.ass = ito.Tass[iV]*d2 + ito.Tass[iV+1]*d1;
	ito.taua = ito.Ttaua[iV]*d2 + ito.Ttaua[iV+1]*d1;

	ito.iss = ito.Tiss[iV]*d2 + ito.Tiss[iV+1]*d1;
	if(var.celltype==3){
		ito.depi = ito.Tdepi[iV]*d2 + ito.Tdepi[iV+1]*d1;
	} else {
		ito.depi = 1.0;
	}
	ito.taui_fast = ito.depi*(ito.Ttaui_fast[iV]*d2 + ito.Ttaui_fast[iV+1]*d1);
	ito.taui_slow = ito.depi*(ito.Ttaui_slow[iV]*d2 + ito.Ttaui_slow[iV+1]*d1);
	ito.Ai_fast = ito.TAi_fast[iV]*d2 + ito.TAi_fast[iV+1]*d1;
	ito.Ai_slow = 1.0 - ito.Ai_fast;
	ito.i = ito.Ai_fast*var.u[11][site][patch] + ito.Ai_slow*var.u[12][site][patch];

	ito.aCaMKss = ito.TaCaMKss[iV]*d2 + ito.TaCaMKss[iV+1]*d1;
	//ito.taua_CaMK = ito.taua;

	//ito.iCaMKss = ito.iss;
	ito.deltaCaMK_dev = ito.TdeltaCaMK_dev[iV]*d2 + ito.TdeltaCaMK_dev[iV+1]*d1;
	ito.deltaCaMK_rec = ito.TdeltaCaMK_rec[iV]*d2 + ito.TdeltaCaMK_rec[iV+1]*d1;
	ito.taui_CaMK_fast = ito.taui_fast*ito.deltaCaMK_dev*ito.deltaCaMK_rec;
	ito.taui_CaMK_slow = ito.taui_slow*ito.deltaCaMK_dev*ito.deltaCaMK_rec;
	//ito.iCaMK = ito.Ai_CaMK_fast*var.u[14] + ito.Ai_CaMK_slow*var.u[15];
	ito.iCaMK = ito.Ai_CaMK_fast*var.u[14][site][patch] + ito.Ai_CaMK_slow*var.u[15][site][patch];

	ito.pCaMK[site][patch] = 1.0/(1.0 + CaMK.Km/CaMK.active[site][patch]);
	ito.ik[site][patch] = ito.Gto*(var.u[0][site][patch]-var.Ek[site][patch])*((1.0-ito.pCaMK[site][patch])*var.u[10][site][patch]*ito.i 
								+ ito.pCaMK[site][patch]*var.u[13][site][patch]*ito.iCaMK);

}


// L-type calcium current
void comp_ical(int site, int patch)
{

	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (var.u[0][site][patch]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

// VDA
	ical.dss = ical.Tdss[iV]*d2 + ical.Tdss[iV+1]*d1;
	ical.taud = ical.Ttaud[iV]*d2 + ical.Ttaud[iV+1]*d1;
// VDI 
	ical.fss = ical.Tfss[iV]*d2 + ical.Tfss[iV+1]*d1;
	ical.tauf_fast = ical.Ttauf_fast[iV]*d2 + ical.Ttauf_fast[iV+1]*d1;
	ical.tauf_slow = ical.Ttauf_slow[iV]*d2 + ical.Ttauf_slow[iV+1]*d1;
	ical.f = ical.Af_fast*var.u[17][site][patch] + ical.Af_slow*var.u[18][site][patch];

// CDI 
	//ical.fcass = ical.fss;
	ical.taufca_fast = ical.Ttaufca_fast[iV]*d2 + ical.Ttaufca_fast[iV+1]*d1;
	ical.taufca_slow = ical.Ttaufca_slow[iV]*d2 + ical.Ttaufca_slow[iV+1]*d1;
	ical.Afca_fast = ical.TAfca_fast[iV]*d2 + ical.TAfca_fast[iV+1]*d1;
	ical.Afca_slow = 1.0 - ical.Afca_fast;
	ical.fca = ical.Afca_fast*var.u[19][site][patch] + ical.Afca_slow*var.u[20][site][patch];

	//ical.jcass = ical.fcass;
	//ical.taujca = 75.0; // (ms).

// CaMK(VDI)
    //ical.f_CaMKss = ical.fss;
	ical.tauf_CaMK_fast = 2.5*ical.tauf_fast;
	//ical.Af_CaMK_fast = ical.Af_fast;
	//ical.Af_CaMK_slow = ical.Af_slow;
	//ical.f_CaMK_slow = var.u[18]; //(var.u[18] = f_slow )
	ical.f_CaMK = ical.Af_fast*var.u[22][site][patch] + ical.Af_slow*var.u[18][site][patch];

// CaMK(CDI)
	
	//ical.fca_CaMKss = ical.fss;
	ical.taufca_CaMK_fast = 2.5*ical.taufca_fast;
	//ical.Afca_CaMK_fast = ical.Afca_fast;
	//ical.Afca_CaMK_slow = ical.Afca_slow;
	//var.fca_CaMK_slow = var.u[20]; //(var.u[20] = fca_slow )
	ical.fca_CaMK = ical.Afca_fast*var.u[23][site][patch] + ical.Afca_slow*var.u[20][site][patch];

	ical.km2n = 1.0*var.u[21][site][patch]; //(jca: recovery from CDI for ICaL )
	ical.tmp = (1.0+ical.kmn/var.u[38][site][patch])*(1.0+ical.kmn/var.u[38][site][patch])*(1.0+ical.kmn/var.u[38][site][patch])*(1.0+ical.kmn/var.u[38][site][patch]);
	ical.alpha_n = 1.0/((ical.kp2n/ical.km2n)+ical.tmp);
	
	ical.exp_Ca = ical.Texp_Ca[iV]*d2 + ical.Texp_Ca[iV+1]*d1;
	ical.exp_Na = ical.Texp_Na[iV]*d2 + ical.Texp_Na[iV+1]*d1;
	ical.exp_K = ical.Texp_K[iV]*d2 + ical.Texp_K[iV+1]*d1;

	ical.phi_ca[site][patch] = (4.0*F*F*var.u[0][site][patch]/R/T)*(ical.gacai*var.u[38][site][patch]*ical.exp_Ca-ical.gacao*var.cao)/(ical.exp_Ca-1.0);
	ical.phi_na[site][patch] = (1.0*F*F*var.u[0][site][patch]/R/T)*(ical.ganai*var.u[34][site][patch]*ical.exp_Na-ical.ganao*var.nao)/(ical.exp_Na-1.0);
	ical.phi_k[site][patch] = (1.0*F*F*var.u[0][site][patch]/R/T)*(ical.gaki*var.u[36][site][patch]*ical.exp_K-ical.gako*var.ko)/(ical.exp_K-1.0);
	ical.ibarcal[site][patch] = ical.pca*ical.phi_ca[site][patch];
	ical.ibarcana[site][patch] = ical.pcana*ical.phi_na[site][patch];
	ical.ibarcak[site][patch] = ical.pcak*ical.phi_k[site][patch];

	ical.ibarcal_CaMK[site][patch] = ical.pca_CaMK*ical.phi_ca[site][patch];
	ical.ibarcana_CaMK[site][patch] = ical.pcana_CaMK*ical.phi_na[site][patch];
	ical.ibarcak_CaMK[site][patch] = ical.pcak_CaMK*ical.phi_k[site][patch];

	ical.pCaMK[site][patch] = 1.0/(1.0 + CaMK.Km/CaMK.active[site][patch]);

	ical.ica[site][patch] =ical.ibarcal[site][patch]*var.u[16][site][patch]*(1.0-ical.pCaMK[site][patch])
							*(ical.f*(1.0-var.u[24][site][patch])+ical.fca*var.u[24][site][patch]*var.u[21][site][patch])
							+ical.ibarcal_CaMK[site][patch]*var.u[16][site][patch]*ical.pCaMK[site][patch]*(ical.f_CaMK*(1.0-var.u[24][site][patch])
							+ical.fca_CaMK*var.u[24][site][patch]*var.u[21][site][patch]);

	ical.icana[site][patch] = ical.ibarcana[site][patch]*var.u[16][site][patch]*(1.0- ical.pCaMK[site][patch])
							*(ical.f*(1.0-var.u[24][site][patch])+ical.fca*var.u[24][site][patch]*var.u[21][site][patch])
							+ical.ibarcana_CaMK[site][patch]*var.u[16][site][patch]*ical.pCaMK[site][patch]*(ical.f_CaMK*(1.0-var.u[24][site][patch])
							+ical.fca_CaMK*var.u[24][site][patch]*var.u[21][site][patch]);

	ical.icak[site][patch] = ical.ibarcak[site][patch]*var.u[16][site][patch]*(1.0- ical.pCaMK[site][patch])
							*(ical.f*(1.0-var.u[24][site][patch])+ical.fca*var.u[24][site][patch]*var.u[21][site][patch])
							+ical.ibarcak_CaMK[site][patch]*var.u[16][site][patch]*ical.pCaMK[site][patch]*(ical.f_CaMK*(1.0-var.u[24][site][patch])
							+ical.fca_CaMK*var.u[24][site][patch]*var.u[21][site][patch]);

}

// Rapidly Activating Potassium Current 
void comp_ikr (int site, int patch)
{
	MKL_INT iV=0;	
	double V1,V2,d1,d2;
	
	V1 = (var.u[0][site][patch]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	// normal (unfacilitated fraction)
	ikr.xrss = ikr.Txrss[iV]*d2 + ikr.Txrss[iV+1]*d1;
	ikr.tauxr_fast = ikr.Ttauxr_fast[iV]*d2 + ikr.Ttauxr_fast[iV+1]*d1;
	ikr.tauxr_slow = ikr.Ttauxr_slow[iV]*d2 + ikr.Ttauxr_slow[iV+1]*d1;
	ikr.Axr_fast = ikr.TAxr_fast[iV]*d2 + ikr.TAxr_fast[iV+1]*d1;
	ikr.Axr_slow = 1.0 - ikr.Axr_fast;
	ikr.rkr = ikr.Trkr[iV]*d2 + ikr.Trkr[iV+1]*d1;

	ikr.xr = ikr.Axr_fast*var.u[25][site][patch] + ikr.Axr_slow*var.u[26][site][patch];

	// facilitated fraction
	ikr.xrss2 = ikr.Txrss2[iV]*d2 + ikr.Txrss2[iV+1]*d1;
	ikr.tauxr_fast2 = ikr.tauxr_fast;
	ikr.tauxr_slow2 = ikr.tauxr_slow;
	ikr.Axr_fast2 = ikr.Axr_fast;
	ikr.Axr_slow2 = 1.0 - ikr.Axr_fast2;
	ikr.rkr2 = ikr.rkr;

	ikr.xr2 = ikr.Axr_fast2*var.u[41][site][patch] + ikr.Axr_slow2*var.u[42][site][patch];

	// Drug effects
	if(var.drug_type == 0){
		ikr.without[site][patch] = var.normal_block*ikr.Gkr*ikr.rategkr*ikr.xr*ikr.rkr*(var.u[0][site][patch]-var.Ek[site][patch]);
		ikr.with[site][patch] = 0.0;
	} else if(var.drug_type == 1){
		ikr.without[site][patch] = var.normal_block*var.fraction_facil*ikr.Gkr*ikr.rategkr*ikr.xr*ikr.rkr*(var.u[0][site][patch]-var.Ek[site][patch]);
		ikr.with[site][patch] = var.normal_block*(1.0-var.fraction_facil)*ikr.Gkr*ikr.rategkr*ikr.xr2*ikr.rkr2*(var.u[0][site][patch]-var.Ek[site][patch]);
	}

	ikr.ik[site][patch] = ikr.with[site][patch]+ikr.without[site][patch];

}

// Slowly Activating Potassium Current 
void comp_iks (int site, int patch)
{
	
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (var.u[0][site][patch]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	iks.xs1ss = iks.Txs1ss[iV]*d2 + iks.Txs1ss[iV+1]*d1;
	//iks.xs2ss = iks.xs1ss;
	iks.tauxs1 = iks.Ttauxs1[iV]*d2 + iks.Ttauxs1[iV+1]*d1;
	iks.tauxs2 = iks.Ttauxs2[iV]*d2 + iks.Ttauxs2[iV+1]*d1;
	iks.KsCa = 1.0+0.6/(1.0+pow(3.8e-5/var.u[37][site][patch],1.4));
	iks.ik[site][patch] = iks.Gks*iks.KsCa*var.u[27][site][patch]*var.u[28][site][patch]*(var.u[0][site][patch]-var.Eks[site][patch]);

}

// Inward rectifier potassium current (Ik1)
void comp_ik1 (int site, int patch)
{
        
	MKL_INT iV=0;   
	double V1,V2,d1,d2;
	V1 = (var.u[0][site][patch]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ik1.k1ss = ik1.Tk1ss[iV]*d2 + ik1.Tk1ss[iV+1]*d1;
	ik1.tauk1 = ik1.Ttauk1[iV]*d2 + ik1.Ttauk1[iV+1]*d1;
	ik1.rk1 = ik1.Trk1[iV]*d2 + ik1.Trk1[iV+1]*d1;

	ik1.ik[site][patch] = ik1.Gk1*ik1.rategk1*var.u[29][site][patch]*ik1.rk1*(var.u[0][site][patch]-var.Ek[site][patch]);

}

// Sodium-Calcium Exchanger V-S

void comp_inaca (int site, int patch)
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (var.u[0][site][patch]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	var.hca=var.Thca[iV]*d2 + var.Thca[iV+1]*d1;
	var.hna=var.Thna[iV]*d2 + var.Thna[iV+1]*d1;

	// intracallular space
	ncxi.h1 = 1.0 + (var.u[33][site][patch]/var.kna3)*(1.0 + var.hna);
	ncxi.h2 = var.u[33][site][patch]*var.hna/(var.kna3*ncxi.h1);
	ncxi.h3 = 1.0/ncxi.h1;
	ncxi.h4 = 1.0 + (var.u[33][site][patch]/var.kna1)*(1.0+var.u[33][site][patch]/var.kna2);
	ncxi.h5 = var.u[33][site][patch]*var.u[33][site][patch]/(ncxi.h4*var.kna1*var.kna2);
	ncxi.h6 = 1.0/ncxi.h4;
	ncxi.h7 = 1.0+(var.nao/var.kna3)*(1.0+1.0/var.hna);
	ncxi.h8 = var.nao/(var.kna3*var.hna*ncxi.h7);
	ncxi.h9 = 1.0/ncxi.h7;
	ncxi.h10 = var.kasym+1.0+(var.nao/var.kna1)*(1.0+var.nao/var.kna2);
	ncxi.h11 = var.nao*var.nao/(ncxi.h10*var.kna1*var.kna2);
	ncxi.h12 = 1.0/ncxi.h10;

	ncxi.k1 = ncxi.h12*var.cao*var.kca_on;
	ncxi.k2 = var.kca_off;
	ncxi.k31 = ncxi.h9*var.omega_ca;
	ncxi.k32 = ncxi.h8*var.omega_naca;
	ncxi.k3 = ncxi.k31 + ncxi.k32;
	ncxi.k41 = ncxi.h3*var.omega_ca/var.hca;
	ncxi.k42 = ncxi.h2*var.omega_naca;
	ncxi.k4 = ncxi.k41 + ncxi.k42;
	ncxi.k5 = var.kca_off;
	ncxi.k6 = ncxi.h6*var.u[37][site][patch]*var.kca_on;
	ncxi.k7 = ncxi.h5*ncxi.h2*var.omega_na;
	ncxi.k8 = ncxi.h8*ncxi.h11*var.omega_na;
	
	ncxi.x1 = ncxi.k2*ncxi.k4*(ncxi.k7+ncxi.k6)+ncxi.k5*ncxi.k7*(ncxi.k2+ncxi.k3);
	ncxi.x2 = ncxi.k1*ncxi.k7*(ncxi.k4+ncxi.k5)+ncxi.k4*ncxi.k6*(ncxi.k1+ncxi.k8);
	ncxi.x3 = ncxi.k1*ncxi.k3*(ncxi.k7+ncxi.k6)+ncxi.k8*ncxi.k6*(ncxi.k2+ncxi.k3);
	ncxi.x4 = ncxi.k2*ncxi.k8*(ncxi.k4+ncxi.k5)+ncxi.k3*ncxi.k5*(ncxi.k1+ncxi.k8);

	ncxi.E1 = ncxi.x1/(ncxi.x1+ncxi.x2+ncxi.x3+ncxi.x4);	
	ncxi.E2 = ncxi.x2/(ncxi.x1+ncxi.x2+ncxi.x3+ncxi.x4);	
	ncxi.E3 = ncxi.x3/(ncxi.x1+ncxi.x2+ncxi.x3+ncxi.x4);	
	ncxi.E4 = ncxi.x4/(ncxi.x1+ncxi.x2+ncxi.x3+ncxi.x4);	

	ncxi.allo = 1.0/(1.0+(var.km_ca_act/var.u[37][site][patch])*(var.km_ca_act/var.u[37][site][patch]));

	ncxi.jnaca_na = 3.0*(ncxi.E4*ncxi.k7 - ncxi.E1*ncxi.k8) + ncxi.E3*ncxi.k42 - ncxi.E2*ncxi.k32;
	ncxi.jnaca_ca = ncxi.E2*ncxi.k2 - ncxi.E1*ncxi.k1;
	var.inaca_i[site][patch] = var.Gnaca*0.8*ncxi.allo*(var.zna*ncxi.jnaca_na + var.zca*ncxi.jnaca_ca);
	
	// subspace
	ncxss.h1 = 1.0 + (var.u[34][site][patch]/var.kna3)*(1.0 + var.hna);
	ncxss.h2 = var.u[34][site][patch]*var.hna/(var.kna3*ncxss.h1);
	ncxss.h3 = 1.0/ncxss.h1;
	ncxss.h4 = 1.0+var.u[34][site][patch]/var.kna1*(1.0+var.u[34][site][patch]/var.kna2);
	ncxss.h5 = var.u[34][site][patch]*var.u[34][site][patch]/(ncxss.h4*var.kna1*var.kna2);
	ncxss.h6 = 1.0/ncxss.h4;
	ncxss.h7 = 1.0+var.nao/var.kna3*(1.0+1.0/var.hna);
	ncxss.h8 = var.nao/(var.kna3*var.hna*ncxss.h7);
	ncxss.h9 = 1.0/ncxss.h7;
	ncxss.h10 = var.kasym+1.0+var.nao/var.kna1*(1.0+var.nao/var.kna2);
	ncxss.h11 = var.nao*var.nao/(ncxss.h10*var.kna1*var.kna2);
	ncxss.h12 = 1.0/ncxss.h10;

	ncxss.k1 = ncxss.h12*var.cao*var.kca_on;
	ncxss.k2 = var.kca_off;
	ncxss.k31 = ncxss.h9*var.omega_ca;
	ncxss.k32 = ncxss.h8*var.omega_naca;
	ncxss.k3 = ncxss.k31 + ncxss.k32;
	ncxss.k41 = ncxss.h3*var.omega_ca/var.hca;
	ncxss.k42 = ncxss.h2*var.omega_naca;
	ncxss.k4 = ncxss.k41 + ncxss.k42;
	ncxss.k5 = var.kca_off;
	ncxss.k6 = ncxss.h6*var.u[38][site][patch]*var.kca_on;
	ncxss.k7 = ncxss.h5*ncxss.h2*var.omega_na;
	ncxss.k8 = ncxss.h8*ncxss.h11*var.omega_na;

	ncxss.x1 = ncxss.k2*ncxss.k4*(ncxss.k7+ncxss.k6)+ncxss.k5*ncxss.k7*(ncxss.k2+ncxss.k3);
	ncxss.x2 = ncxss.k1*ncxss.k7*(ncxss.k4+ncxss.k5)+ncxss.k4*ncxss.k6*(ncxss.k1+ncxss.k8);
	ncxss.x3 = ncxss.k1*ncxss.k3*(ncxss.k7+ncxss.k6)+ncxss.k8*ncxss.k6*(ncxss.k2+ncxss.k3);
	ncxss.x4 = ncxss.k2*ncxss.k8*(ncxss.k4+ncxss.k5)+ncxss.k3*ncxss.k5*(ncxss.k1+ncxss.k8);

	ncxss.E1 = ncxss.x1/(ncxss.x1+ncxss.x2+ncxss.x3+ncxss.x4);	
	ncxss.E2 = ncxss.x2/(ncxss.x1+ncxss.x2+ncxss.x3+ncxss.x4);	
	ncxss.E3 = ncxss.x3/(ncxss.x1+ncxss.x2+ncxss.x3+ncxss.x4);	
	ncxss.E4 = ncxss.x4/(ncxss.x1+ncxss.x2+ncxss.x3+ncxss.x4);	

	ncxss.allo = 1.0/(1.0+(var.km_ca_act/var.u[38][site][patch])*(var.km_ca_act/var.u[38][site][patch]));

	ncxss.jnaca_na = 3.0*(ncxss.E4*ncxss.k7 - ncxss.E1*ncxss.k8) + ncxss.E3*ncxss.k42 - ncxss.E2*ncxss.k32;
	ncxss.jnaca_ca = ncxss.E2*ncxss.k2 - ncxss.E1*ncxss.k1;
	var.inaca_ss[site][patch] = var.Gnaca*0.2*ncxss.allo*(var.zna*ncxss.jnaca_na + var.zca*ncxss.jnaca_ca);

    var.inaca[site][patch] = var.inaca_i[site][patch]+var.inaca_ss[site][patch];
    
}

// Sodium-Potassium Pump

void comp_inak (int site, int patch)
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (var.u[0][site][patch]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	inak.knai = inak.Tknai[iV]*d2 + inak.Tknai[iV+1]*d1;
	inak.knao = inak.Tknao[iV]*d2 + inak.Tknao[iV+1]*d1;

	inak.P = inak.SigP/(1.0 + (inak.H/inak.HP) + (var.u[33][site][patch]/inak.nap) + (var.u[35][site][patch]/inak.kp));

	inak.a1 = inak.kp1*(var.u[33][site][patch]/inak.knai)*(var.u[33][site][patch]/inak.knai)*(var.u[33][site][patch]/inak.knai)/
				((1.0+var.u[33][site][patch]/inak.knai)*(1.0+var.u[33][site][patch]/inak.knai)*(1.0+var.u[33][site][patch]/inak.knai)+(1+var.u[35][site][patch]/inak.kki)*(1+var.u[35][site][patch]/inak.kki)-1.0);
	inak.b2 = inak.km2*(var.nao/inak.knao)*(var.nao/inak.knao)*(var.nao/inak.knao)/
	 			((1.0+var.nao/inak.knao)*(1.0+var.nao/inak.knao)*(1.0+var.nao/inak.knao)+(1.0+var.ko/inak.kko)*(1.0+var.ko/inak.kko)-1.0);
	inak.a3 = inak.kp3*(var.ko/inak.kko)*(var.ko/inak.kko)/
				((1.0+var.nao/inak.knao)*(1.0+var.nao/inak.knao)*(1.0+var.nao/inak.knao)+(1.0+var.ko/inak.kko)*(1.0+var.ko/inak.kko)-1.0);
	inak.b3 = inak.km3*inak.P*inak.H/(1.0+inak.MgATP/inak.k_MgATP);
	inak.b4 = inak.km4*(var.u[35][site][patch]/inak.kki)*(var.u[35][site][patch]/inak.kki)/
				((1.0+var.u[33][site][patch]/inak.knai)*(1.0+var.u[33][site][patch]/inak.knai)*(1.0+var.u[33][site][patch]/inak.knai)+(1.0+var.u[35][site][patch]/inak.kki)*(1.0+var.u[35][site][patch]/inak.kki)-1.0);

	inak.x1 = inak.a4*inak.a1*inak.a2+inak.b2*inak.b4*inak.b3+inak.a2*inak.b4*inak.b3+inak.b3*inak.a1*inak.a2;
	inak.x2 = inak.b2*inak.b1*inak.b4+inak.a1*inak.a2*inak.a3+inak.a3*inak.b1*inak.b4+inak.a2*inak.a3*inak.b4;
	inak.x3 = inak.a2*inak.a3*inak.a4+inak.b3*inak.b2*inak.b1+inak.b2*inak.b1*inak.a4+inak.a3*inak.a4*inak.b1;
	inak.x4 = inak.b4*inak.b3*inak.b2+inak.a3*inak.a4*inak.a1+inak.b2*inak.a4*inak.a1+inak.b3*inak.b2*inak.a1;

	inak.E1 = inak.x1/(inak.x1+inak.x2+inak.x3+inak.x4);
	inak.E2 = inak.x2/(inak.x1+inak.x2+inak.x3+inak.x4);
	inak.E3 = inak.x3/(inak.x1+inak.x2+inak.x3+inak.x4);
	inak.E4 = inak.x4/(inak.x1+inak.x2+inak.x3+inak.x4);

	inak.jna = 3.0*(inak.E1*inak.a3 - inak.E2*inak.b3);
	inak.jk = 2.0*(inak.E4*inak.b1 - inak.E3*inak.a1);

	inak.inak[site][patch] = inak.G*(var.zna*inak.jna + var.zk*inak.jk);

}

// Sarcolemmal Ca Pump 

void comp_ipca (int site, int patch)
{

	ipca.ca[site][patch] = ipca.G*var.u[37][site][patch]/(ipca.km + var.u[37][site][patch]);

}

// K Background Current
void comp_ikb (int site, int patch)
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (var.u[0][site][patch]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	ikb.xkb = ikb.Txkb[iV]*d2 + ikb.Txkb[iV+1]*d1;

	ikb.k[site][patch] = ikb.G*ikb.xkb*(var.u[0][site][patch] - var.Ek[site][patch]);
}

// Ca Background Current 

void comp_icab (int site, int patch)
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (var.u[0][site][patch]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	icab.exp = icab.Texp[iV]*d2 + icab.Texp[iV+1]*d1;
	
	icab.ca[site][patch] = icab.pcab*4.0*var.u[0][site][patch]*F*F/R/T*(icab.gacai*var.u[37][site][patch]*icab.exp-icab.gacao*var.cao)/(icab.exp-1.0);

}

// Na Background Current 

void comp_inab (int site, int patch)
{
	MKL_INT iV=0;
	double V1,V2,d1,d2;
	V1 = (var.u[0][site][patch]+Emax)*dvm;
	V2 = (int)V1;
	d1 = V1-V2;
	d2 = 1.0-d1;
	iV = (int)V2;

	inab.exp = inab.Texp[iV]*d2 + inab.Texp[iV+1]*d1;

	inab.na[site][patch] = inab.pnab*var.u[0][site][patch]*F*F/R/T*(var.u[33][site][patch]*inab.exp-var.nao)/(inab.exp-1.0);

}

void comp_CaMK (int site, int patch)
{
	CaMK.bound[site][patch] = CaMK.z*(1.0-var.u[30][site][patch])/(1.0+CaM.Km/var.u[38][site][patch]);

	CaMK.active[site][patch] = CaMK.bound[site][patch] + var.u[30][site][patch];

}

void comp_diff (int site, int patch)
{

	jdiff.na[site][patch] = (var.u[34][site][patch]-var.u[33][site][patch])/jdiff.tau_na;
	jdiff.k[site][patch] = (var.u[36][site][patch]-var.u[35][site][patch])/jdiff.tau_k;
	jdiff.ca[site][patch] = (var.u[38][site][patch]-var.u[37][site][patch])/jdiff.tau_ca;

}

void comp_jrel (int site, int patch)
{

	jrel.NPss[site][patch] = jrel.p*(-jrel.a*ical.ica[site][patch]/(1.0+pow((1.5/var.u[40][site][patch]),8.0)));
	jrel.tau_NP[site][patch] = jrel.b_tau/(1.0+(0.0123/var.u[40][site][patch]));
	if (jrel.tau_NP[site][patch] < 0.001){
		jrel.tau_NP[site][patch] = 0.001;
		//printf("jrel tauNP is small\n");
	}

	jrel.CaMKss[site][patch] = jrel.p*(-jrel.a_CaMK*ical.ica[site][patch]/(1.0+pow((1.5/var.u[40][site][patch]),8.0)));
	jrel.tau_CaMK[site][patch] = jrel.b_tau_CaMK/(1.0+(0.0123/var.u[40][site][patch]));
	if (jrel.tau_CaMK[site][patch] < 0.001){
		jrel.tau_CaMK[site][patch] = 0.001;
	}

	jrel.pCaMK[site][patch] = 1.0/(1.0+CaMK.Km/CaMK.active[site][patch]);

	jrel.ca[site][patch] = (1.0-jrel.pCaMK[site][patch])*var.u[31][site][patch]+jrel.pCaMK[site][patch]*var.u[32][site][patch];

}

void comp_jup (int site, int patch)
{

	jup.np[site][patch] = jup.p*(0.004375*var.u[37][site][patch]/(0.00092+var.u[37][site][patch]));
	jup.CaMK[site][patch] = jup.p*((1.0+jup.dCaMK)*0.004375*var.u[37][site][patch]/(0.00092 - jup.dKm_PLB + var.u[37][site][patch]));
	
	jup.pCaMK[site][patch] = 1.0/(1.0 + CaMK.Km/CaMK.active[site][patch]);

	jup.leak[site][patch] = 0.0039375*var.u[39][site][patch]/15.0;

	jup.ca[site][patch] = (1.0-jup.pCaMK[site][patch])*jup.np[site][patch] + jup.pCaMK[site][patch]*jup.CaMK[site][patch] - jup.leak[site][patch];

}

void comp_jtr (int site, int patch)
{

	jtr.ca[site][patch] = (var.u[39][site][patch]-var.u[40][site][patch])/jtr.tau;

}

void comp_concentration (int site, int patch)
{
	var.b_Ca_i[site][patch] = 1.0/(1.0+var.cmdnbar*var.kmcmdn/((var.kmcmdn+var.u[37][site][patch])*(var.kmcmdn+var.u[37][site][patch]))+var.trpnbar*var.kmtrpn/((var.kmtrpn+var.u[37][site][patch])*(var.kmtrpn+var.u[37][site][patch])));
	var.b_Ca_ss[site][patch] = 1.0/(1.0+var.bsrbar*var.kmbsr/((var.kmbsr+var.u[38][site][patch])*(var.kmbsr+var.u[38][site][patch]))+var.bslbar*var.kmbsl/((var.kmbsl+var.u[38][site][patch])*(var.kmbsl+var.u[38][site][patch])));
	var.b_Ca_jsr[site][patch] = 1.0/(1.0+var.csqnbar*var.kmcsqn/((var.kmcsqn+var.u[40][site][patch])*(var.kmcsqn+var.u[40][site][patch])));

}


// Reversal potentials */

void comp_reversal_potential(int site, int patch)
{
	var.Ena[site][patch] = var.RTonF*log(var.nao/var.u[33][site][patch]);
	var.Ek[site][patch] = var.RTonF*log(var.ko/var.u[35][site][patch]);
	var.Eks[site][patch] = var.RTonF*log((var.ko+var.prnak*var.nao)/(var.u[35][site][patch]+var.prnak*var.u[33][site][patch]));
	
	//printf("Ena=%lf, Ek=%lf, Eks=%lf\n",var.Ena,var.Ek,var.Eks);
}

