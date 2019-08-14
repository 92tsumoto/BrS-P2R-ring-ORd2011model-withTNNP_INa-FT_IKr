//#ifndef __SYSPARA_H_INCLUDE 
//#define __SYSPARA_H_INCLUDE

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "mkl.h"

#define NN 43
#define BUF 200
#define NUM 60
#define MEDIA_SITE 600
#define MEDIA_PATCH 3
#define MAT_SIZE MEDIA_SITE*MEDIA_PATCH
#define beats 30

//#define R 8314.472
//#define F 96485.33771638995
//#define T 310
#define R 8314.0		// J/kmol/K
#define F 96485.0	// C/mol
#define T 310.0		// K
#define Q 3         // Tempurture coefficent

#define dvm 5
#define Emax 2000
#define Emin -2000
#define VNMAX (Emax-Emin)*dvm+1

struct varstruct {

	int datas;
	int line_wid[NUM];
	
	int n;
	double Rmyo,Gj,D,Istim;
	double coef,dist;
	double Ri,Rg,Rd,Rj;
	double gx0,ggx;
	double gj_rate;

	// for EF mechanism
	double gamma[MEDIA_SITE-1],delta[MEDIA_SITE-1];

	// used variable in linear resistance matrix
	double p1,p2,p3;
	double ri,rd,rj;
	double w1,w2,w3,w4;

	// membrane surface area
	double s1,s2[MEDIA_PATCH];

	// An invariant constant
	double RTonF,RTon2F;

	// Drug concentration for Nifekalant
    double drug_conc,ic50,ec50,hill,ehill;
	double normal_block,fraction_facil;
	double block_rate;

	// Model & Cell type & Drug type
	// model switch
	int celltype;
	int couple_type; // 0:Gj alone, 1:EF+Gj
	int model_type; // 0:original, 1:Dutta et al., PBMB 2017
	int drug_type;	// 0: without facilitation, 1: with facilitation
	int drug_mode;	// 0: concentration, 1: rate block

	// Na channel localize switch
	int na_switch,na_switch2;
	double narate,jmrate1,lmrate1,jmrate2,lmrate2;
	//double vshift;

	// state variables
	double x0[NUM][NN];
	//double u[NN][MEDIA_SITE][MEDIA_PATCH];
	double ***u;
	double **inner_v,**inner_current;
	double **dvdt;
	//double inner_v[MEDIA_SITE][MEDIA_PATCH];
	//double inner_current[MEDIA_SITE][MEDIA_PATCH];
	//double dvdt[MEDIA_SITE][MEDIA_PATCH];
	double disk_membrane[MEDIA_SITE][2];

	// quantifies of the cleft potential
	double cleft_width,cleft_R;
	double cleft_axial[MEDIA_SITE-1][2],cleft_potential[MEDIA_SITE-1];

	// Cell Geometry
	double RGC;
	double length,a;
	double ageo[MEDIA_PATCH];
	double acap[MEDIA_PATCH];
	double vcell[MEDIA_PATCH];
	double vmyo[MEDIA_PATCH];
	double vss[MEDIA_PATCH];
	double vsr[MEDIA_PATCH];
	double vnsr[MEDIA_PATCH];
	double vjsr[MEDIA_PATCH];
	double vcleft[MEDIA_PATCH];

	double vr1[MEDIA_PATCH];
	double vr2[MEDIA_PATCH];
	double vr3[MEDIA_PATCH];
	double vr4[MEDIA_PATCH];
	double vr5[MEDIA_PATCH];
	double vr6[MEDIA_PATCH];
	double vr7[MEDIA_PATCH];
	double vr8[MEDIA_PATCH];

// Ion Valences 
	double zna,zk,zca;

// Ion concentrations
	double nao,ko,cao;

// Reversal potentials
	double prnak;
	//double Ena[MEDIA_SITE][MEDIA_PATCH];
	//double Ek[MEDIA_SITE][MEDIA_PATCH];
	//double Eks[MEDIA_SITE][MEDIA_PATCH];
	double **Ena,**Ek,**Eks;

// IKr adding variable
	double Tf,gkr;

// Sodium-Calcium Exchanger V-S
	double *Thca,*Thna;
	
	double hca,hna;
	double kna1,kna2,kna3,kasym;
	double omega_na,omega_ca,omega_naca;
	double kca_on,kca_off,qna,qca;
	double km_ca_act;
	double Gnaca;
	//double inaca_i[MEDIA_SITE][MEDIA_PATCH];
	//double inaca_ss[MEDIA_SITE][MEDIA_PATCH];
	//double inaca[MEDIA_SITE][MEDIA_PATCH];
	double **inaca_i,**inaca_ss,**inaca;

// Total Ion currents 
	//double Ina_i_total[MEDIA_SITE][MEDIA_PATCH], Ina_ss_total[MEDIA_SITE][MEDIA_PATCH];
	//double Ik_i_total[MEDIA_SITE][MEDIA_PATCH], Ik_ss_total[MEDIA_SITE][MEDIA_PATCH];
	//double Ica_i_total[MEDIA_SITE][MEDIA_PATCH], Ica_ss_total[MEDIA_SITE][MEDIA_PATCH];
	//double Itotal[MEDIA_SITE][MEDIA_PATCH];
	double **Ina_i_total,**Ina_ss_total;
	double **Ik_i_total,**Ik_ss_total;
	double **Ica_i_total,**Ica_ss_total;
	double **Itotal;

// for measurement of entire INa current per cell
	double *total_ina;

// Ion concentration and Buffers
	double cmdnbar,kmcmdn;
	double trpnbar,kmtrpn;
	double bsrbar,kmbsr;
	double bslbar,kmbsl;
	double csqnbar,kmcsqn;
	//double b_Ca_i[MEDIA_SITE][MEDIA_PATCH],b_Ca_ss[MEDIA_SITE][MEDIA_PATCH],b_Ca_jsr[MEDIA_SITE][MEDIA_PATCH];
	double **b_Ca_i,**b_Ca_ss,**b_Ca_jsr;

// test variable
	double dt;

// Base Currnt Stimulus
	int stim_sw;
	double Istim_base;

// Stimulus parameters
	double BCL;  // Base cycle length = stimulus period
	int beat; // Number of stimulus

    int m;
    int l;

    double tsign[NUM];
    double tend[NUM];

    int write,out_data,out_data_plus,out_data_plus2;
    int half;

// linear algebra // Do not change!!
	double *AT,*AT2;
// with Pardiso rutine
	double *b, *solv;
	long long *ia, *ja, *stok_ia, *stok_ja;
	int par_k;
	long long mtype,nrhs,iparm[64],maxfct,mnum,phase;
	long long error, msglvl;
	double ddum;
	long long idum;

/* Internal solver memory pointer pt, */
/* 32-bit: int pt[64]; 64-bit: long int pt[64] */
/* or void *pt[64] should be OK on both architectures */
	void *pt[64];

// using Eular
//double k1[NN][MEDIA_SITE][MEDIA_PATCH];
	double ***k1;


} var;


// Fast and Late sodium currnets
struct inastruct {

	double *Tmss,*Ttaum,*Thss,*Ttauh_fast,*Ttauh_slow,*Tjss,*Ttauj;
	double *ThCaMKss;

	double mss,taum,hss,tauh_fast,tauh_slow,jss,tauj;
	double Ah_fast,Ah_slow;
	double hCaMKss,tauh_CaMK_slow,jCaMKss,tauj_CaMK;
	double Ah_CaMK_fast,Ah_CaMK_slow;
	double h,hCaMK;
	
	double Gna_fast,gnaf_all,gnaf_site;
	//double fast_pCaMK[MEDIA_SITE][MEDIA_PATCH];
	//double gnaf_local[MEDIA_SITE][MEDIA_PATCH];
	double **fast_pCaMK,**gnaf_local;

	// Late sodium currnet
	double *Tmlss,*Ttauml,*Thlss,*ThlCaMKss;

	double mlss,tauml,hlss,tauhl,hlCaMKss,tauhl_CaMK;

	double Gna_late,gnal_all,gnal_site;
	//double late_pCaMK[MEDIA_SITE][MEDIA_PATCH];
	//double gnal_local[MEDIA_SITE][MEDIA_PATCH];
	double **late_pCaMK,**gnal_local;

	// ina
	//double fast[MEDIA_SITE][MEDIA_PATCH];
	//double late[MEDIA_SITE][MEDIA_PATCH];
	//double total[MEDIA_SITE][MEDIA_PATCH];
	double **fast,**late,**total;

} ina;

// Transient Outward Current (Ito)
struct itostruct {

	double *Tass,*Ttaua,*Tiss,*Ttaui_fast,*Ttaui_slow,*TAi_fast,*Tdepi;
	double *TaCaMKss,*TdeltaCaMK_dev,*TdeltaCaMK_rec;

	double ass,taua,iss,taui_fast,taui_slow,Ai_fast,Ai_slow,i,depi;
	double aCaMKss,taua_CaMK,iCaMKss,deltaCaMK_dev,deltaCaMK_rec;
	double taui_CaMK_fast,taui_CaMK_slow,Ai_CaMK_fast,Ai_CaMK_slow,iCaMK;
	double Gto;

	//double pCaMK[MEDIA_SITE][MEDIA_PATCH];
	//double ik[MEDIA_SITE][MEDIA_PATCH];
	double **pCaMK,**ik;

} ito;

// L-type Calcium channel current (IcaL)
struct icalstruct {
	
	double rate;
	double *Texp_Ca,*Texp_Na,*Texp_K;
	double *Tdss,*Ttaud,*Tfss,*Ttauf_fast,*Ttauf_slow;
	double *Ttaufca_fast,*Ttaufca_slow,*TAfca_fast; 
	
	double dss,taud,fss,tauf_fast,tauf_slow;
	double f,Af_fast,Af_slow;

	double fcass,taufca_fast,taufca_slow,Afca_fast,Afca_slow,fca; 

	double jcass,taujca; 

	double f_CaMKss,f_CaMK,tauf_CaMK_fast,Af_CaMK_fast,Af_CaMK_slow,f_CaMK_slow; 

	double fca_CaMKss,fca_CaMK,taufca_CaMK_fast,Afca_CaMK_fast,Afca_CaMK_slow,fca_CaMK_slow; 

	double kmn,kp2n,km2n,alpha_n,tmp;
	double pca,gacai,gacao;
	double exp_Ca,exp_Na,exp_K;
	double pcana,ganai,ganao; 
	double pcak,gaki,gako;
	double pca_CaMK,pcana_CaMK,pcak_CaMK;

	//double phi_ca[MEDIA_SITE][MEDIA_PATCH],ibarcal[MEDIA_SITE][MEDIA_PATCH],ibarcal_CaMK[MEDIA_SITE][MEDIA_PATCH];
	//double phi_na[MEDIA_SITE][MEDIA_PATCH],ibarcana[MEDIA_SITE][MEDIA_PATCH],ibarcana_CaMK[MEDIA_SITE][MEDIA_PATCH]; 
	//double phi_k[MEDIA_SITE][MEDIA_PATCH],ibarcak[MEDIA_SITE][MEDIA_PATCH],ibarcak_CaMK[MEDIA_SITE][MEDIA_PATCH];
	//double pCaMK[MEDIA_SITE][MEDIA_PATCH];
	//double ica[MEDIA_SITE][MEDIA_PATCH],icana[MEDIA_SITE][MEDIA_PATCH],icak[MEDIA_SITE][MEDIA_PATCH];
	double **phi_ca,**ibarcal,**ibarcal_CaMK;
	double **phi_na,**ibarcana,**ibarcana_CaMK; 
	double **phi_k,**ibarcak,**ibarcak_CaMK;
	double **pCaMK;
	double **ica,**icana,**icak;

} ical;

// Rapid activating potassium current (Ikr)
struct ikrstruct {

	double *Txrss,*Ttauxr_fast,*Ttauxr_slow,*TAxr_fast;
	double *Trkr;
	double *Txrss2;

	double xr,xrss,tauxr_fast,tauxr_slow,Axr_fast,Axr_slow,rkr;
	double xr2,xrss2,tauxr_fast2,tauxr_slow2,Axr_fast2,Axr_slow2,rkr2;
	double Gkr,rategkr;

	//double ik[MEDIA_SITE][MEDIA_PATCH];
	double **ik,**without,**with;

} ikr;

// Slowlactivating potassium current (Iks)
struct iksstruct {

	double *Txs1ss,*Ttauxs1,*Ttauxs2;
		
	double xs1ss,tauxs1,xs2ss,tauxs2;
	double Gks,KsCa;

	//double ik[MEDIA_SITE][MEDIA_PATCH];
	double **ik;

} iks;

// Inward rectifier potassium current (Ik1)
struct ik1struct {

	double *Tk1ss,*Ttauk1,*Trk1;
	
	double k1ss,tauk1,rk1;
	double Gk1,rategk1;

	//double ik[MEDIA_SITE][MEDIA_PATCH];
	double **ik;

} ik1;

struct ncxistruct {

	double h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12;
	double k1,k2,k3,k31,k32,k4,k41,k42;
	double k5,k6,k7,k8;
	double x1,x2,x3,x4;
	double E1,E2,E3,E4;
	double allo;
	double jnaca_na,jnaca_ca;

} ncxi;

struct ncxssstruct {

	double h1,h2,h3,h4,h5,h6,h7,h8,h9,h10,h11,h12;
	double k1,k2,k3,k31,k32,k4,k41,k42;
	double k5,k6,k7,k8;
	double x1,x2,x3,x4;
	double E1,E2,E3,E4;
	double allo;
	double jnaca_na,jnaca_ca;

} ncxss;

// Sodium-Potassium Pump
struct inakstruct {

// Sodium-Potassium Pump
//	double fnak,sigma,ibarnak,kmnai,kmko;
//	double inak[MEDIA_SITE][MEDIA_PATCH];
	
	double *Tknai,*Tknao;

	double kp1,km1,kp2,km2;
	double kp3,km3,kp4,km4;
	double ko_nai,ko_nao,delta;
	double knai,knao;
	double kki,kko,MgADP,MgATP,k_MgATP,H,SigP;
	double HP,nap,kp,P;
	double a1,a2,a3,a4;
	double b1,b2,b3,b4;
	double x1,x2,x3,x4;
	double E1,E2,E3,E4;
	double G,jna,jk;

	//double inak[MEDIA_SITE][MEDIA_PATCH];
	double **inak;

} inak;

// Sarcolemmal Ca Pump
struct ipcastruct {

	double G,km;

	//double ca[MEDIA_SITE][MEDIA_PATCH];
	double **ca;

} ipca;

// K Background Current
struct ikbstruct {

	double *Txkb,**k;

	double xkb;
	double G;
	//double k[MEDIA_SITE][MEDIA_PATCH];

} ikb;

// Ca Background Current
struct icabstruct {

	double *Texp,**ca;

	double exp;
	double pcab,gacai,gacao;
	
	//double ca[MEDIA_SITE][MEDIA_PATCH];

} icab;

// Na Background Current
struct inabstruct {

	double *Texp,**na;

	double exp;
	double pnab;
	
	//double na[MEDIA_SITE][MEDIA_PATCH];


} inab;

// Ca/Calmodulin dependent protein kinese (CaMK)
struct CaMKstruct{
	
	double a,b,z;
	double Km;

	//double bound[MEDIA_SITE][MEDIA_PATCH],active[MEDIA_SITE][MEDIA_PATCH];
	double **bound,**active;

} CaMK;

// Ca/Calmodulin dependent protein kinese (CaMK)
struct CaMstruct{
	
	double Km;

} CaM;

// SR calcium release flux, via RyR (Jrel)
struct jrelstruct {

	double b_tau,a,b_tau_CaMK,a_CaMK;
	double p;

	//double NPss[MEDIA_SITE][MEDIA_PATCH],tau_NP[MEDIA_SITE][MEDIA_PATCH];
	//double CaMKss[MEDIA_SITE][MEDIA_PATCH],tau_CaMK[MEDIA_SITE][MEDIA_PATCH];
	//double pCaMK[MEDIA_SITE][MEDIA_PATCH],ca[MEDIA_SITE][MEDIA_PATCH];
	double **NPss,**tau_NP,**CaMKss,**tau_CaMK;
	double **pCaMK,**ca;

} jrel;

// Calcium uptake via SERCA pump
struct jupstruct {

	double p;
	double dKm_PLB,dCaMK;
	//double np[MEDIA_SITE][MEDIA_PATCH],CaMK[MEDIA_SITE][MEDIA_PATCH];
	//double pCaMK[MEDIA_SITE][MEDIA_PATCH],ca[MEDIA_SITE][MEDIA_PATCH],leak[MEDIA_SITE][MEDIA_PATCH];
	double **np,**CaMK,**pCaMK,**ca,**leak;

} jup;

// diffusion flux
struct jdiffstruct {

	double tau_na,tau_k,tau_ca;
	//double na[MEDIA_SITE][MEDIA_PATCH],k[MEDIA_SITE][MEDIA_PATCH],ca[MEDIA_SITE][MEDIA_PATCH];
	double **na,**k,**ca;

} jdiff;

// Translocation of Ca Ions from NSR to JSR
struct jtrstruct {

	double tau;
	
	//double ca[MEDIA_SITE][MEDIA_PATCH];
	double **ca;

} jtr;

// for main
void data_out(FILE *, FILE *, FILE *, FILE *, FILE *, FILE *, FILE *);
void make_ExPTable();
void eular(int n, double h, double t);
void function(int site, int patch, double t);
void input_para(FILE *);
void static_paras(FILE *);
void initial_mem();
void close_mem();
void *malloc2d(size_t size, int row, int col);
void *malloc3d(size_t size, int i, int j, int k);

void comp_reversal_potential(int site, int patch);
void comp_ina(int site, int patch);
void comp_ito(int site, int patch);
void comp_ical(int site, int patch);
void comp_ikr(int site, int patch);
void comp_iks(int site, int patch);
void comp_ik1(int site, int patch);
void comp_inaca(int site, int patch);
void comp_inak(int site, int patch);
void comp_ipca(int site, int patch);
void comp_ikb(int site, int patch);
void comp_icab(int site, int patch);
void comp_inab(int site, int patch);
void comp_CaMK(int site, int patch);
void comp_diff(int site, int patch);
void comp_jrel(int site, int patch);
void comp_jup(int site, int patch);
void comp_jtr (int site, int patch);
void comp_concentration (int site, int patch);

// for dataout
void vm_data(FILE *, double time);
void ina_data(FILE *, double time);
void ical_data(FILE *, double time);
//void cai_data(FILE *, double time);
//void ik_data(FILE *, FILE *, FILE *, FILE *, double time);
//void ik_data(FILE *, double time);
//void ilca_data(FILE *, double time);
//void itotal_data(FILE *, double time);
//void intra_i_data(FILE *, double time);
//void intra_v_data(FILE *, double time);
//void cleft_v_data(FILE *, FILE *, double time);
//void dvdt_data(FILE *, double time);
void vm1_data(FILE *, double time);
void vm2_data(FILE *, double time);

MKL_INT init_pardiso();
void linear_coeff();
void linear_vec();
void linear_vec3();
void linear_solve();
void release_of_memory();
//void dgesv_(int *N, int *NRHS, double AT[], int *LDA, int pivot[], double b[], int *LDB, int *ok);

//void main(int argc, char **argv);

