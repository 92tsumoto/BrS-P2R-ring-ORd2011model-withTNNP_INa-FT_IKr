#include "syspara.h"

typedef double Number;
typedef long long Lint;

void initial_mem()
{
	int i,k,m,l,z;
	int site,patch;

	site=MEDIA_SITE;
	patch=MEDIA_PATCH;

// initialized memory for state variables and arguments for Eular method 
	var.u = (double***)malloc3d(sizeof(double),NN,site,patch);
	var.k1 = (double***)malloc3d(sizeof(double),NN,site,patch);

	var.inner_v = (double**)malloc2d(sizeof(double*),site,patch);
	var.inner_current = (double**)malloc2d(sizeof(double*),site,patch);
	var.dvdt = (double**)malloc2d(sizeof(double*),site,patch);

// reversal potentail
	var.Ena = (double**)malloc2d(sizeof(double*),site,patch);
	var.Ek = (double**)malloc2d(sizeof(double*),site,patch);
	var.Eks = (double**)malloc2d(sizeof(double*),site,patch);

//	Incx
	var.inaca_i = (double**)malloc2d(sizeof(double*),site,patch);
	var.inaca_ss = (double**)malloc2d(sizeof(double*),site,patch);
	var.inaca = (double**)malloc2d(sizeof(double*),site,patch);

// Ionic currents
	var.Ina_i_total = (double**)malloc2d(sizeof(double*),site,patch);
	var.Ina_ss_total = (double**)malloc2d(sizeof(double*),site,patch);
	var.Ik_i_total = (double**)malloc2d(sizeof(double*),site,patch);
	var.Ik_ss_total = (double**)malloc2d(sizeof(double*),site,patch);
	var.Ica_i_total = (double**)malloc2d(sizeof(double*),site,patch);
	var.Ica_ss_total = (double**)malloc2d(sizeof(double*),site,patch);
	var.Itotal = (double**)malloc2d(sizeof(double*),site,patch);

// Ca buffer
	var.b_Ca_i = (double**)malloc2d(sizeof(double*),site,patch);
	var.b_Ca_ss = (double**)malloc2d(sizeof(double*),site,patch);
	var.b_Ca_jsr = (double**)malloc2d(sizeof(double*),site,patch);

// initialized tablization memorys for Exp functions

	// ina_fast
	ina.Tmss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttaum=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Thss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttauh_fast=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttauh_slow=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Tjss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttauj=(Number *)calloc(VNMAX,sizeof(Number));
	ina.ThCaMKss=(Number *)calloc(VNMAX,sizeof(Number));
	if( ina.Tmss==NULL || ina.Ttaum==NULL || ina.Thss==NULL || ina.Ttauh_fast==NULL || ina.Ttauh_slow==NULL 
		|| ina.Tjss==NULL || ina.Ttauj==NULL || ina.ThCaMKss==NULL ) exit(1);
	// ina_late
	ina.Tmlss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Ttauml=(Number *)calloc(VNMAX,sizeof(Number));
	ina.Thlss=(Number *)calloc(VNMAX,sizeof(Number));
	ina.ThlCaMKss=(Number *)calloc(VNMAX,sizeof(Number));
	if( ina.Tmlss==NULL || ina.Ttauml==NULL || ina.Thlss==NULL || ina.ThlCaMKss==NULL ) exit(1);

	ina.fast_pCaMK = (double**)malloc2d(sizeof(double*),site,patch);
	ina.gnaf_local = (double**)malloc2d(sizeof(double*),site,patch);
	ina.late_pCaMK = (double**)malloc2d(sizeof(double*),site,patch);
	ina.gnal_local = (double**)malloc2d(sizeof(double*),site,patch);
	ina.fast = (double**)malloc2d(sizeof(double*),site,patch);
	ina.late = (double**)malloc2d(sizeof(double*),site,patch);
	ina.total = (double**)malloc2d(sizeof(double*),site,patch);

	// ito
	ito.Tass=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Ttaua=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Tiss=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Ttaui_fast=(Number *)calloc(VNMAX,sizeof(Number));
	ito.Ttaui_slow=(Number *)calloc(VNMAX,sizeof(Number));
	ito.TAi_fast=(Number *)calloc(VNMAX,sizeof(Number));
	ito.TaCaMKss=(Number *)calloc(VNMAX,sizeof(Number));
	ito.TdeltaCaMK_dev=(Number *)calloc(VNMAX,sizeof(Number));
	ito.TdeltaCaMK_rec=(Number *)calloc(VNMAX,sizeof(Number));
	if(var.celltype == 3){
		ito.Tdepi=(Number *)calloc(VNMAX,sizeof(Number));
		if(ito.Tdepi == NULL) exit(1);
	}
	if( ito.Tass==NULL || ito.Ttaua==NULL 
		|| ito.Tiss==NULL || ito.Ttaui_fast==NULL || ito.Ttaui_slow==NULL
		|| ito.TAi_fast==NULL || ito.TaCaMKss==NULL
		|| ito.TdeltaCaMK_dev==NULL || ito.TdeltaCaMK_rec==NULL ) exit(1);
	
	ito.pCaMK = (double**)malloc2d(sizeof(double*),site,patch);
	ito.ik = (double**)malloc2d(sizeof(double*),site,patch);

	// ical
	ical.Tdss=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttaud=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Tfss=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttauf_fast=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttauf_slow=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttaufca_fast=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Ttaufca_slow=(Number *)calloc(VNMAX,sizeof(Number));
	ical.TAfca_fast=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Texp_Ca=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Texp_Na=(Number *)calloc(VNMAX,sizeof(Number));
	ical.Texp_K=(Number *)calloc(VNMAX,sizeof(Number));
	if( ical.Tdss==NULL || ical.Ttaud==NULL 
		|| ical.Tfss==NULL || ical.Ttauf_fast==NULL || ical.Ttauf_slow==NULL
		|| ical.Ttaufca_fast == NULL || ical.Ttaufca_slow==NULL || ical.TAfca_fast==NULL 
		|| ical.Texp_Ca==NULL || ical.Texp_Na==NULL || ical.Texp_K == NULL ) exit(1);
	
	ical.phi_ca = (double**)malloc2d(sizeof(double*),site,patch);
	ical.ibarcal = (double**)malloc2d(sizeof(double*),site,patch);
	ical.ibarcal_CaMK = (double**)malloc2d(sizeof(double*),site,patch);
	ical.phi_na = (double**)malloc2d(sizeof(double*),site,patch);
	ical.ibarcana = (double**)malloc2d(sizeof(double*),site,patch);
	ical.ibarcana_CaMK = (double**)malloc2d(sizeof(double*),site,patch);
	ical.phi_k = (double**)malloc2d(sizeof(double*),site,patch);
	ical.ibarcak = (double**)malloc2d(sizeof(double*),site,patch);
	ical.ibarcak_CaMK = (double**)malloc2d(sizeof(double*),site,patch);
	ical.pCaMK = (double**)malloc2d(sizeof(double*),site,patch);
	ical.ica = (double**)malloc2d(sizeof(double*),site,patch);
	ical.icana = (double**)malloc2d(sizeof(double*),site,patch);
	ical.icak = (double**)malloc2d(sizeof(double*),site,patch);

	// ikr
	ikr.Txrss=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.Txrss2=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.Ttauxr_fast=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.Ttauxr_slow=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.TAxr_fast=(Number *)calloc(VNMAX,sizeof(Number));
	ikr.Trkr=(Number *)calloc(VNMAX,sizeof(Number));
	if( ikr.Txrss==NULL || ikr.Txrss2==NULL || ikr.Ttauxr_fast==NULL || ikr.Ttauxr_slow==NULL 
		|| ikr.TAxr_fast==NULL || ikr.Trkr==NULL ) exit(1);
	ikr.without = (double**)malloc2d(sizeof(double*),site,patch);
	ikr.with = (double**)malloc2d(sizeof(double*),site,patch);
	ikr.ik = (double**)malloc2d(sizeof(double*),site,patch);

	// iks
	iks.Txs1ss=(Number *)calloc(VNMAX,sizeof(Number));
	iks.Ttauxs1=(Number *)calloc(VNMAX,sizeof(Number));
	iks.Ttauxs2=(Number *)calloc(VNMAX,sizeof(Number));
	if( iks.Txs1ss==NULL || iks.Ttauxs1==NULL || iks.Ttauxs2==NULL ) exit(1);
	iks.ik = (double**)malloc2d(sizeof(double*),site,patch);

	// ik1
	ik1.Tk1ss=(Number *)calloc(VNMAX,sizeof(Number));
	ik1.Ttauk1=(Number *)calloc(VNMAX,sizeof(Number));
	ik1.Trk1=(Number *)calloc(VNMAX,sizeof(Number));
	if( ik1.Tk1ss == NULL || ik1.Ttauk1==NULL || ik1.Trk1==NULL ) exit(1);
	ik1.ik = (double**)malloc2d(sizeof(double*),site,patch);

	// inaca
	var.Thca=(Number *)calloc(VNMAX,sizeof(Number));
	var.Thna=(Number *)calloc(VNMAX,sizeof(Number));
	if( var.Thca==NULL || var.Thna==NULL ) exit(1);

	// inak
	inak.Tknai=(Number *)calloc(VNMAX,sizeof(Number));
	inak.Tknao=(Number *)calloc(VNMAX,sizeof(Number));
	if( inak.Tknai==NULL || inak.Tknao==NULL ) exit(1);
	inak.inak = (double**)malloc2d(sizeof(double*),site,patch);
	
	// ipca
	ipca.ca = (double**)malloc2d(sizeof(double*),site,patch);

	// ikb	
	ikb.Txkb=(Number *)calloc(VNMAX,sizeof(Number));
	if( ikb.Txkb==NULL ) exit(1);
	ikb.k = (double**)malloc2d(sizeof(double*),site,patch);

	// icab
	icab.Texp=(Number *)calloc(VNMAX,sizeof(Number));
	if( icab.Texp==NULL ) exit(1);
	icab.ca = (double**)malloc2d(sizeof(double*),site,patch);

	// inab
	inab.Texp=(Number *)calloc(VNMAX,sizeof(Number));
	if( inab.Texp==NULL ) exit(1);
	inab.na = (double**)malloc2d(sizeof(double*),site,patch);

	// Ca/CaMK
	CaMK.bound = (double**)malloc2d(sizeof(double*),site,patch);
	CaMK.active = (double**)malloc2d(sizeof(double*),site,patch);
	
	// SR calcium release flux, via RyR (Jrel)
	jrel.NPss = (double**)malloc2d(sizeof(double*),site,patch);
	jrel.tau_NP = (double**)malloc2d(sizeof(double*),site,patch);
	jrel.CaMKss = (double**)malloc2d(sizeof(double*),site,patch);
	jrel.tau_CaMK = (double**)malloc2d(sizeof(double*),site,patch);
	jrel.pCaMK = (double**)malloc2d(sizeof(double*),site,patch);
	jrel.ca = (double**)malloc2d(sizeof(double*),site,patch);

	// Calcium uptake via SERCA pump
	jup.np = (double**)malloc2d(sizeof(double*),site,patch);
	jup.CaMK = (double**)malloc2d(sizeof(double*),site,patch);
	jup.pCaMK = (double**)malloc2d(sizeof(double*),site,patch);
	jup.ca = (double**)malloc2d(sizeof(double*),site,patch);
	jup.leak = (double**)malloc2d(sizeof(double*),site,patch);
	
	// diffusion flux
	jdiff.na = (double**)malloc2d(sizeof(double*),site,patch);
	jdiff.k = (double**)malloc2d(sizeof(double*),site,patch);
	jdiff.ca = (double**)malloc2d(sizeof(double*),site,patch);

	// Translocation of Ca Ions from NSR to JSR
	jtr.ca = (double**)malloc2d(sizeof(double*),site,patch);

}
		
/* -------------------------------------------------------------------- */
/* .. Termination and release of memory. */
/* -------------------------------------------------------------------- */
void close_mem()
{

	int i,j;

		free(var.u);free(var.k1);
		free(var.inner_v);free(var.inner_current);free(var.dvdt);
		free(var.Ena);free(var.Ek);free(var.Eks);
		free(var.inaca_i);free(var.inaca_ss);free(var.inaca);
		free(var.Ina_i_total);free(var.Ina_ss_total);
		free(var.Ik_i_total);free(var.Ik_ss_total);
		free(var.Ica_i_total);free(var.Ica_ss_total);
		free(var.Itotal);
		free(var.b_Ca_i);free(var.b_Ca_ss);free(var.b_Ca_jsr);

	// ina_fast
		free(ina.Tmss); free(ina.Ttaum); free(ina.Thss); free(ina.Ttauh_fast); free(ina.Ttauh_slow); free(ina.Tjss); free(ina.Ttauj); 
		free(ina.ThCaMKss); free(ina.fast_pCaMK); free(ina.late_pCaMK);
		free(ina.gnaf_local); free(ina.gnal_local); free(ina.fast);free(ina.late);free(ina.total);

	// ina_late
		free(ina.Tmlss); free(ina.Ttauml); free(ina.Thlss); free(ina.ThlCaMKss);

	// ito
		free(ito.Tass); free(ito.Ttaua); free(ito.Tiss); free(ito.Ttaui_fast); free(ito.Ttaui_slow);
		free(ito.TAi_fast); free(ito.TaCaMKss); if(var.celltype==3){free(ito.Tdepi);}
		free(ito.TdeltaCaMK_dev); free(ito.TdeltaCaMK_rec);
		free(ito.pCaMK);free(ito.ik);

	// ical
		free(ical.Tdss); free(ical.Ttaud); free(ical.Tfss); free(ical.Ttauf_fast); free(ical.Ttauf_slow);
		free(ical.Ttaufca_fast); free(ical.Ttaufca_slow); free(ical.TAfca_fast);
		free(ical.Texp_Ca); free(ical.Texp_Na); free(ical.Texp_K);
		free(ical.phi_ca);free(ical.ibarcal);free(ical.ibarcal_CaMK);
		free(ical.phi_na);free(ical.ibarcana);free(ical.ibarcana_CaMK);
		free(ical.phi_k);free(ical.ibarcak);free(ical.ibarcak_CaMK);
		free(ical.pCaMK);
		free(ical.ica);free(ical.icana);free(ical.icak);

	// ikr
		free(ikr.Txrss); free(ikr.Ttauxr_fast); free(ikr.Ttauxr_slow); 
		free(ikr.Txrss2);
		free(ikr.TAxr_fast); free(ikr.Trkr);
		free(ikr.ik); free(ikr.without); free(ikr.with);

	// iks
		free(iks.Txs1ss); free(iks.Ttauxs1); free(iks.Ttauxs2);
		free(iks.ik);
	
	// ik1
		free(ik1.Tk1ss); free(ik1.Ttauk1); free(ik1.Trk1);
		free(ik1.ik);

	// inaca
		free(var.Thca); free(var.Thna);

	// inak
		free(inak.Tknai); free(inak.Tknao);
		free(inak.inak);

	// ipca
		free(ipca.ca);

	// ikb	
		free(ikb.Txkb);free(ikb.k);

	// icab	
		free(icab.Texp);free(icab.ca);

	// inab	
		free(inab.Texp);free(inab.na);

	// Ca/CaMK
		free(CaMK.bound);free(CaMK.active);

	// SR calcium release flux, via RyR (Jrel)
		free(jrel.NPss);free(jrel.tau_NP);
		free(jrel.CaMKss);free(jrel.tau_CaMK);
		free(jrel.pCaMK);free(jrel.ca);

	// Calcium uptake via SERCA pump
		free(jup.np);free(jup.CaMK);free(jup.pCaMK);free(jup.ca);free(jup.leak);

	// diffusion flux
		free(jdiff.na);free(jdiff.k);free(jdiff.ca);

	// Translocation of Ca Ions from NSR to JSR
		free(jtr.ca);


}

// メモリ領域が連続な2次元配列
void *malloc2d(size_t size, int row, int col)
{
	char **a, *b;
	int  t = size * col;
	int i;

	// インデックスと要素を一気に確保
	a = (char**)malloc((sizeof(*a) + t) * row);

	if (a) {
		// [インデックス, インデックス, ..., 要素, 要素, 要素, ...]と整列させるため要素の開始位置をずらす
		b = (char*)(a + row);

		// 各行の先頭アドレスを与える
		for (i = 0; i < row; i++) {
			a[i] = b;
			b += t; // 要素のサイズ×列の長さの分だけずらす
		}
		return a;
	}

	return NULL;
}

// メモリ領域が連続な3次元配列
void *malloc3d(size_t size, int i, int j, int k)
{
	char ***a, **b, *c;
	int  t = size * k;
	int idx1,idx2;

	// インデックスと要素を一気に確保
	a = (char***)malloc((sizeof(*a) + sizeof(**a) * j + t * j) * i);

	if (a) {
		b = (char**)(a + i);
		c = (char *)(b + i * j);

		for (idx1 = 0; idx1 < i; idx1++) {
			a[idx1] = b;
			for (idx2 = 0; idx2 < j; idx2++) {
				b[idx2] = c;
				c += t;
			}
			b += j;
		}

		return a;
	}

	return NULL;
}
