
#include "syspara.h"

void make_ExpTable()
{

	int vindex;
	double v;
	double fact,fdeact,sact,sdeact;
	double am,bm,ah,bh,aj,bj;
    
	for(vindex=0;vindex<VNMAX;vindex++){

        v = (double)vindex/dvm-(double)Emax;
		
        /** for fast INa & Late INa **/
		if(var.model_type==1){ // O'Hara original model
			ina.Tmss[vindex] = 1.0/(1.0+exp(-(v+39.57)/9.871));
			ina.Ttaum[vindex] = 1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
			ina.Thss[vindex] = 1.0/(1.0+exp((v+82.90)/6.086));
			ina.Ttauh_fast[vindex] = 1.0/(1.432E-5*exp(-(v+1.196)/6.285)+6.149*exp((v+0.5096)/20.27));
			ina.Ttauh_slow[vindex] = 1.0/(0.009764*exp(-(v+17.95)/28.05)+0.3343*exp((v+5.730)/56.66));
			ina.Ttauj[vindex] = 2.038+1.0/(0.02136*exp(-(v+100.6)/8.281)+0.3052*exp((v+0.9941)/38.45));
			ina.ThCaMKss[vindex] = 1.0/(1.0+exp((v+89.1)/6.086));
			// for inal 
			ina.Tmlss[vindex] = 1.0/(1.0 + exp(-(v+42.85)/5.264));
			ina.Ttauml[vindex] = 1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
			ina.Thlss[vindex] = 1.0/(1.0 + exp((v + 87.61)/7.488));
			ina.ThlCaMKss[vindex] = 1.0/(1.0 + exp((v+93.81)/7.488));
		} else if(var.model_type==2){ // Dutta's modification (Passini et al. JMCC, 2016)
			ina.Tmss[vindex] = 1.0/(1.0+exp(-(v+39.57)/9.871));
			ina.Ttaum[vindex] = 1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
			ina.Thss[vindex] = 1.0/(1.0+exp((v+78.5)/6.22));
			ina.Ttauh_fast[vindex] = 1.0/(3.686E-6*exp(-(v+3.8875)/7.8579)+16.0*exp((v-0.4963)/9.1843));
			ina.Ttauh_slow[vindex] = 1.0/(0.009764*exp(-(v+17.95)/28.05)+0.3343*exp((v+5.730)/56.66));
			ina.Ttauj[vindex] = 4.8590+1.0/(0.8628*exp(-(v+116.7258)/7.6005)+1.1096*exp((v+6.2719)/9.0358));
			ina.ThCaMKss[vindex] = 1.0/(1.0+exp((v+84.7)/6.22));
			// for inal 
			ina.Tmlss[vindex] = 1.0/(1.0 + exp(-(v+42.85)/5.264));
			ina.Ttauml[vindex] = 1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
			ina.Thlss[vindex] = 1.0/(1.0 + exp((v + 87.61)/7.488));
			ina.ThlCaMKss[vindex] = 1.0/(1.0 + exp((v+93.81)/7.488));
		} else if(var.model_type==3){ // TP06 INa substitution (Elshef&Cherry, PONE, 2014)
			// TP INa activation
			ina.Tmss[vindex] = 1.0/((1.0+exp(-(v+56.86)/9.03))*(1.0+exp(-(v+56.86)/9.03)));
			am=1.0/(1.+exp((-60.0-v)/5.0));
			bm=0.1/(1.0+exp((v+35.0)/5.0))+0.10/(1.0+exp((v-50.0)/200.0));
			ina.Ttaum[vindex] = am*bm;
			// TP INa inactivation
			ina.Thss[vindex] = 1.0/((1.0+exp((v+71.55)/7.43))*(1.0+exp((v+71.55)/7.43)));
			if(v>=-40.0){
				ah = 0.0;
				bh = 0.77/(0.13*(1.0+exp(-(v+10.66)/11.1)));
			} else {
				ah = 0.057*exp(-(v+80.0)/6.8);
				bh = 2.7*exp(0.079*v)+3.1E+5*exp(0.3485*v);
			}
			// TP INa recovery 
			ina.Ttauh_fast[vindex] = 1.0/(ah+bh);
			ina.Ttauh_slow[vindex] = 1E+99;
			ina.Tjss[vindex] = 1.0/((1.0+exp((v+71.55)/7.43))*(1.0+exp((v+71.55)/7.43)));
			if(v>=-40.0){
				aj = 0.0;
				bj = 0.6*exp(0.057*v)/(1.0+exp(-0.1*(v+32.0)));
			} else {
				aj = (-2.5428E+4*exp(0.2444*v)-6.948E-6*exp(-0.04391*v))*(v+37.78)/(1.0+exp(0.311*(v+79.23)));
				bj = 0.02424*exp(-0.01052*v)/(1.0+exp(-0.1378*(v+40.14)));
			}
			ina.Ttauj[vindex] = 1.0/(aj+bj);
			ina.ThCaMKss[vindex] = 1.0/(1.0+exp((v+89.1)/6.086));
			// for inal 
			ina.Tmlss[vindex] = 1.0/(1.0 + exp(-(v+42.85)/5.264));
			ina.Ttauml[vindex] = 1.0/(6.765*exp((v+11.64)/34.77)+8.552*exp(-(v+77.42)/5.955));
			ina.Thlss[vindex] = 1.0/(1.0 + exp((v + 87.61)/7.488));
			ina.ThlCaMKss[vindex] = 1.0/(1.0 + exp((v+93.81)/7.488));
		}

		// ito
		ito.Tass[vindex] = 1.0/(1.0+exp(-(v-14.34)/14.82));
		ito.Ttaua[vindex] = 1.0515/(1.0/(1.2089*(1.0+exp(-(v-18.4099)/29.3814)))+3.5/(1.0+exp((v+100.0)/29.3814)));

		ito.Tiss[vindex] = 1.0/(1.0+exp((v+43.94)/5.711));
		if(var.celltype==3){
			ito.Tdepi[vindex] = 1.0-(0.95/(1.0+exp((v+70.0)/5.0)));
		} 
		ito.Ttaui_fast[vindex] = 4.562+1.0/(0.3933*exp(-(v+100.0)/100.0)+0.08004*exp((v+50.0)/16.59));
		ito.Ttaui_slow[vindex] = 23.62+1.0/(0.001416*exp(-(v+96.52)/59.05)+1.7808E-8*exp((v+114.1)/8.079));
		ito.TAi_fast[vindex] = 1.0/(1.0+exp((v-213.6)/151.2));

		ito.TaCaMKss[vindex] = 1.0/(1.0+exp(-(v-24.34)/14.82));
		ito.TdeltaCaMK_dev[vindex] = 1.354+1.0E-4/(exp((v-167.4)/15.89)+exp(-(v-12.23)/0.2154));
		ito.TdeltaCaMK_rec[vindex] = 1.0-0.5/(1.0+exp((v+70.0)/20.0));

		// for ical
		ical.Tdss[vindex] = 1.0/(1.0+exp(-(v+3.940)/4.230));
		ical.Ttaud[vindex] = 0.6+1.0/(exp(-0.05*(v+6.0))+exp(0.09*(v+14.0)));

		ical.Tfss[vindex] = 1.0/(1.0+exp((v+19.58)/3.696));
		ical.Ttauf_fast[vindex] = 7.0+1.0/(0.0045*exp(-(v+20.0)/10.0)+0.0045*exp((v+20.0)/10.0));
		ical.Ttauf_slow[vindex] = 1000.0+1.0/(0.000035*exp(-(v+5.0)/4.0)+0.000035*exp((v+5.0)/6.0));

		ical.Ttaufca_fast[vindex] = 7.0+1.0/(0.04*exp(-(v-4.0)/7.0)+0.04*exp((v-4.0)/7.0));
		ical.Ttaufca_slow[vindex] = 100.0+1.0/(0.00012*exp(-v/3.0)+0.00012*exp(v/7.0));
		ical.TAfca_fast[vindex] = 0.3+0.6/(1.0+exp((v-10.0)/10.0));

		ical.Texp_Ca[vindex] = exp(v/var.RTon2F);
		ical.Texp_Na[vindex] = exp(v/var.RTonF);
		ical.Texp_K[vindex] = exp(v/var.RTonF);
		
		// for ikr 
		// ORd original expressions (non-facilitation component for ikr)
		//ikr.Txrss[vindex] = 1.0/(1.0+exp(-(v+8.337)/6.789));
		//ikr.Ttauxr_fast[vindex] = 12.98+1.0/(0.3652*exp((v-31.66)/3.869)+4.123E-5*exp(-(v-47.78)/20.38));
		//ikr.Ttauxr_slow[vindex] = 1.865+1.0/(0.06629*exp((v-34.70)/7.355)+1.128E-5*exp(-(v-29.74)/25.94));
		//ikr.TAxr_fast[vindex] = 1.0/(1.0+exp((v+54.81)/38.21));
		//ikr.Trkr[vindex] = 1.0/(1.0+exp((v+55.0)/75.0))*1.0/(1.0+exp((v-10.0)/30.0));

		// A modified ikr model using our HEK data
		ikr.Txrss[vindex] = 1.0/(1.0+exp(-(v+10.75)/7.4));
		fact = 1400.0/(1.0+exp((v+10.0)/13.8));
		fdeact = 1300.0/(1.0+exp(-(v+55.0)/12.3));
		ikr.Ttauxr_fast[vindex] = (fact+fdeact)-1300.0;
		sact = 5000.0/(1.0+exp((v+22.8)/11.5))+2000/(1.0+exp((v+25)/2.5));
		sdeact = 4000.0/(1.0+exp(-(v+61.5)/5.0))+3000/(1.0+exp(-(v+54)/2.5));
		if(v > -80.0){
			ikr.Ttauxr_slow[vindex] = (sact+sdeact)-7000.0;
		} else {
			ikr.Ttauxr_slow[vindex] = 2000.0+8000.0/(1.0+exp(-(v+74.9)/10.2));
		}
		ikr.TAxr_fast[vindex] = 0.05+0.95/(1.0+exp((v+73.4)/16.0));
		ikr.Trkr[vindex] = (0.12+0.88/(1.0+exp((v+60.3)/27.0)))*(0.32+0.68/(1.0+exp((v-7.0)/20.0)));
		// for ikr (facilitated fraction)
		ikr.Txrss2[vindex] = 1.0/(1.0+exp(-(v+37.20)/7.7));


		// for iks 
		iks.Txs1ss[vindex] = 1.0/(1.0+exp(-(v+11.60)/8.932));
		iks.Ttauxs1[vindex] = 817.3+1.0/(2.326E-4*exp((v+48.28)/17.80)+0.001292*exp(-(v+210.0)/230.0));
		iks.Ttauxs2[vindex] = 1.0/(0.01*exp((v-50.0)/20.0)+0.0193*exp(-(v+66.54)/31.0));

		// ik1 
		ik1.Tk1ss[vindex] = 1.0/(1.0+exp(-(v+2.5538*var.ko+144.59)/(1.5692*var.ko+3.8115)));
		ik1.Ttauk1[vindex] = 122.2/(exp(-(v+127.2)/20.36)+exp((v+236.8)/69.33));
		ik1.Trk1[vindex] = 1.0/(1.0+exp((v+105.8-2.6*var.ko)/9.493));

		// inaca
		var.Thca[vindex] = exp(var.qca*v/var.RTonF);
		var.Thna[vindex] = exp(var.qna*v/var.RTonF);

		// inak 
		inak.Tknai[vindex] = inak.ko_nai*exp((inak.delta*v*F)/(3.0*R*T));
		inak.Tknao[vindex] = inak.ko_nao*exp(((1.0-inak.delta)*v*F)/(3.0*R*T));

		// ikb
		ikb.Txkb[vindex] = 1.0/(1.0+exp(-(v-14.48)/18.34));

		// icab
		icab.Texp[vindex] = exp(v/var.RTon2F);

		// inab
		inab.Texp[vindex] = exp(v/var.RTonF);
	}

}
