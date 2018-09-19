/*−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−
*
* This code is based on the rabbit ventricular myocyte model of Mahajan et al (Biophysical Journal, 2008) and modified for simulating the mouse ventricular myocyte action potential and calcium transient and finding combinations of ionic conductances that produce a normal electrophysiological phenotype.
*
* Contact Information:
*
* Center for interdisciplinary research on complex systems
* Departments of Physics, Northeastern University
*
* Alain Karma a.karma (at) northeastern.edu
*
* The code was used to produce the results published in
*Title: The Ca2+ transient as a feedback sensor controlling cardiomyocyte ionic conductances in mouse populations
*Authors: Colin Rees, Jun-Hai Yang, Marc Santolini, Aldons J. Lusis, James N. Weiss, Alain Karma
*Journal: eLife
*Year: 2018
*−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−− */

// Information of the original Mahajan rabbit ventricular myocyte model:
/*−−−−−−−−−−−−−−−−−−−− UCLA Model ver 1.00 −−−−−−−−−−−−−−−−−−−−−−
*
* Contact Information
*
* Departments of Medicine (Cardiology)
* David Geffen School of Medicine at UCLA
*
* Daisuke Sato
* Yohannes Shiferaw
* James N Weiss
*
* The code was used to produce simulations in
* A. Mahajan, Y. Shiferaw, D. Sato, A. Baher, R. Olcese, L.−H. Xie,
* M.−J. Yang, P.−S. Chen, J. G. Restrepo, A. Karma, A. Garfinkel,
* Z. Qu, and J. N. Weiss, A rabbit ventricular action potential model
* replicating cardiac dynamics at rapid heart rates, Biophysical Journal, vol 94 (2008), pp. 392-410. 
*−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−−− */ 

#ifndef ___CELL
#define ___CELL
#include "cell.h"
#define ccm 1.0//25

CCell::CCell(void) : y(new double[N]),
	xm(y[0]), xh(y[1]), xj(y[2]), xr(y[3]), 
        nks(y[4]), /**/ /*xs1(y[4]), xs2(y[5]),*/ atos(y[6]), itos(y[7]),
	v(y[8]), cp(y[9]), cs(y[10]), ci(y[11]), 
	cj(y[12]), cjp(y[13]), xir(y[14]), c1(y[15]), 
	c2(y[16]), xi1ca(y[17]), xi1ba(y[18]), xi2ca(y[19]), 
	xi2ba(y[20]), xnai(y[21]), atof(y[22]), itof(y[23]), 
  tropi(y[24]), trops(y[25]), jrel(y[26]), ass(y[27]), iss(y[28]), aur(y[29]), iur(y[30])
{
// initial conditions
	jrel=0;
	xm=0.001145222753;// sodium m-gate
	xh=0.9898351676;// sodium h-gate
	xj=0.9930817518;// soiumj-gate

	xr=0.008709989976;// ikr gate variable 
	
	nks=0.000262753;// iks gate variable
	  /*these two no longe needed
	  xs1=0.08433669901;// iks gate variable
	  xs2=0.1412866149;// iks gate varaible */

	atos=0.000417069;//0.003757746357;// ito slow activation
	itos=0.998543;//0.1553336368;// ito slow inactivation

	ass=0.0004147069;// iss activation
	iss=1;// iss inactivation
	
	aur=0.2;//0.000417069;// iur activation
	iur=0.923;//0.998543;// iur inactivation

	v=-82.402;//-86.79545769; // voltage

	cp=1.682601371;// averaged dyadic space con.
	cs=0.3937;//0.133511;//0.115001;//0.3205609256;// averaged submembrane conc.
	ci=0.334;//0.127428;//0.115001;//0.3863687451;// myoplasm conc.

	cj=105;//88.3;//98.94;// NSR load
	cjp=105;//87;//300 // average JSR load

	xir=0.064;//0.006462569526;// SR current flux

	// Markov gate variables 
	c1=2.0192e-05;//1.925580885e-05;// C1
	c2=0.935194;//0.9535940241;// C2
	xi1ca=0.0109354;//0.007052299702;// I1_Ca
	xi1ba=4.25749e-05;//3.629261123e-05;// I1_Ba
	xi2ca=0.0353134;//0.02316349806;// I2_Ca
	xi2ba=0.018492;//0.01613268649;// I2_Ba
	//TAKE ANOTHER GOOD LOOK AT THIS!!!!!
	xnai=14.2371;//14.01807252;// internal Na conc.

	atof=0.0018;//0.00265563;//0.003737842131;// ito fast activation
	itof=0.996;//0.999977;//0.9823715315;// ito slow inactivation

	tropi=27.43;//29.64807803;// time dependent buffers in myplasm (troponin)
	trops=31.25;//26.37726416;// time dependent buffers in submembrane (troponin)

	hpde=0.1;
	vold = v;
	jparam=1;

#ifdef ___USE_VAR_FOR_CONST
	xnao=140.0;//136.0;//mM	  external Na
	xki=143.5;// mM	 internal K
	xko=5.40;//mM	 external K
	cao=1.8/ccm;// mM	 external Ca

	gca=172.9*1.5*0.68; //16*6;// /1.5 //182;// ica conductance //Relevent
	gtos=0.0;//0.04;// ito slow conductance 
	gtof=0.8*0.2;//0.4067;//0.11;// ito fast conductance 
	gnaca=292.8;//0.84;// exchanger strength 
	gkr=0.078;//0.0125;// Ikr conductance 
	gks=0.0575;//0.32;
	gkss=0.05;
	gkur=0.16;
	gkix=0;//0.3;// Ik1 conductance
	gnak=0.88;//1.5; //RELEVENT /sodiumpotassium
	vup=0.45;//0.3;// uptake strength
	taus=0.75;// diffusional delay (ms)				// ------------CHANGED
	gna=13.0;//  /2  12.0;// sodium conductance (mS/micro F) //RELEVENT (Changed to reproduce peak value, chosen to match experimental values for dV/dt )
	taur=10;//30.0;// spark lifetime (ms)
	taua=20;//100.0;// NSR-JSR diffusional delay (ms)
	av=4;//11.3;
	cstar=90.0;
#endif
}
CCell::~CCell()
{
	delete[] y;
}
void CCell::Prepare(double BCL, int Iter)
{
	if (Iter==0)
	{
		double dciold=0;
		double dciold2=0;
		bool first=false;
		int Tn=BCL*10000/hpde, BCLn=BCL/hpde, Durn=stimduration/hpde;
		for (int tn=0;tn<Tn;tn++)
		{
			double t=tn*hpde;
			if (tn%BCLn < Durn)
			{
				if (first)
				{
					if (fabs(ci-dciold2)<0.00001 && t>BCL*300)
					{
						break;
					}
					dciold2=dciold;
					dciold=ci;
					first=false;
				}
				Pace(stim);
			}
			else
			{
				first=true;
				Pace();
			}
		}
	}
	else
	{
		int Tn=BCL*Iter/hpde, BCLn=BCL/hpde, Durn=stimduration/hpde;
		for (int tn=0;tn<Tn;tn++)
		{
			if (tn%BCLn < Durn)
				Pace(stim);
			else
				Pace();
		}
	}
}
CCell& CCell::operator=(const CCell& cell)
{
	if (&cell!=this)
	{
		for (int i=0;i<N;i++)
		{
			y[i]=cell.y[i];
		}
		jparam=cell.jparam;
		vold=cell.vold;
		hpde=cell.hpde;
	#ifdef ___USE_VAR_FOR_CONST
		xnao=cell.xnao;
		xki=cell.xki;
		xko=cell.xko;
		cao=cell.cao;

		gca=cell.gca;
		gtos=cell.gtos;
		gtof=cell.gtof;
		gnaca=cell.gnaca;
		gkr=cell.gkr;
		gks=cell.gks;
		gkix=cell.gkix;
		gnak=cell.gnak;
		vup=cell.vup;
		taus=cell.taus;
		gna=cell.gna;
		taur=cell.taur;
		taua=cell.taua;
		av=cell.av;
		cstar=cell.cstar;
	#endif
	}
	return(*this);
}
void CCell::ClampAP(double t, double T, double APD)
{
	const double Vmin=-80;//-80mV
	const double Vmax=30;//30mV
	double clampv;
	if (APD==0)
	{
		const double a=2.0/3.0*1000;
		double x=a/(a+T);
		int m=(int)(t/T);
		if (m*T+x*T>t)
		{
			clampv=Vmin+(Vmax-Vmin)*sqrt(1-((t-m*T)/x/T)*((t-m*T)/x/T));
		}
		else
		{
			clampv=Vmin;
		}
	}
	else
	{
		double x=APD/T;
		int m=(int)(t/T);
		if (m*T+x*T>t)
		{
			clampv=Vmin+(Vmax-Vmin)*sqrt(1-((t-m*T)/x/T)*((t-m*T)/x/T));
		}
		else
		{
			clampv=Vmin;
		}
	}

	double dv=(vold-v)/hpde;
	vold=v;
	double Itotal;
	if(fabs(dv)>25.0)// then finer time step when dv/dt large
	{
		hode=hpde/10;
		for (int iii=0;iii<10;iii++)
		{
			v=clampv;
			Itotal=PaceX(0);
		}
	}
	else
	{
		hode=hpde;
		v=clampv;
		Itotal=PaceX(0);
	}
}
double CCell::Pace(double Istim)
{
// -------------time step adjustment ------------------------
	double dv=(vold-v)/hpde;
	vold=v;
	double Itotal;
	if(fabs(dv)>25.0)// then finer time step when dv/dt large
	{
		hode=hpde/10;
		//cout << v << endl;
		for (int iii=0;iii<10;iii++)
		{
		  //v=Vinput[iii];
			Itotal=PaceX(Istim);
		}
	}
	else
	{
		hode=hpde;
		//cout << v << endl;
		//v=Vinput[0];
		Itotal=PaceX(Istim);
	}
	return Itotal;
}
double CCell::PaceVClamp(double clampv)
{
// -------------time step adjustment ------------------------
	double dv=(vold-v)/hpde;
	vold=v;
	double Itotal;
	if(fabs(dv)>25.0)// then finer time step when dv/dt large
	{
		hode=hpde/10;
		for (int iii=0;iii<10;iii++)
		{
			v=clampv;
			Itotal=PaceX(0);
		}
	}
	else
	{
		hode=hpde;
		v=clampv;
		Itotal=PaceX(0);
	}
	return Itotal;
}

double CCell::PaceX(double Istim)
{
	double xik1=comp_ik1();
	double xito=comp_ito();//itos and itof
	double xinak=comp_inak();
	double csm=cs/1000.0;// convert micro M to mM
	double xinacaq=comp_inaca(csm); //RELEVENT
//----------- Equations for Ca cycling -------------------------
	double xdif=(cs-ci)/taus;//diffusion from submembrane to myoplasm
	// Troponin kinetics
	const double xkon=0.0327/ccm;
	const double xkoff=0.0196;
	const double btrop=70.0;
	double xbi=xkon*ci*(btrop-tropi)-xkoff*tropi;
	double xbs=xkon*cs*(btrop-trops)-xkoff*trops;

	double xiup=comp_iup();
	double xileak=comp_ileak();

	double po=comp_icalpo();
	double rxa=comp_rxa(csm);
	double xicaq=gca*po*rxa/ccm;// Ca current in micro M/ms

	double xsvipca=comp_svipca();// PMCA //Not PMCA!! CaB

	//double dcs=comp_inst_buffer(cs)*(50.0*(xir-xdif-xicaq+xinacaq)-xbs);
	double dcs=comp_inst_buffer(cs)*(50.0*(xir-xdif-xicaq+xinacaq)-xbs);  //----with PMCA----
	double dci=comp_inst_buffer(ci)*(xdif-xiup+xileak-xbi-xsvipca/32.0);
	//cout << xicaq << " " << xinacaq << " " << xir << " " << xdif << " " << xiup << " " << xileak << " " << xbi << " " << xbs << " " << xsvipca/16. <<  endl;
	double dcj=-xir+xiup-xileak;// SR load dynamics 
	double dcjp=(cj-cjp)/taua;// NSR-JSR relaxation dynamics 
	double Qr=comp_Q();
	double dir=comp_dir(po, Qr, rxa, dcj);
	double dcp=comp_dcp(po, Qr, rxa);

	double xina=comp_ina();
	double xikr=comp_ikr();
	double xiks=comp_iks();
	
	cp+=dcp*hode;
	cs+=dcs*hode;
	ci+=dci*hode;
	cj+=dcj*hode;
	xir+=dir*hode;
	cjp+=dcjp*hode;

	tropi+=xbi*hode;
	trops+=xbs*hode;

	jrel=xicaq;//-------------------detect--------------------;

//-------convert ion flow to current---------
	const double wca=16.0;//conversion factor between micro molar/ms to micro amps/ micro farads
	double xinaca=wca*xinacaq;
	double xica=2.0*wca*xicaq;
//--------sodium dynamics -------------------------
	const double xrr=(1.0/wca)/1000.0;// note: sodium is in m molar so need to divide by 1000
	xnai+=200*(-xrr*(xina+3.0*xinak+3.0*xinaca))*hode;//this 100 speeds up relaxation to steady state.!!!!
	//cout << xina << " " << xinak << " " << xinaca << " " << xnai << endl;
// --------	dV/dt ------------------------------------
	//double Itotal=(-(xina+xik1+xikr+xiks+xito+xinaca+xica+xinak)+ Istim);
	double Itotal=(-(xina+xik1+xikr+xiks+xito+xinaca+xica+xinak+xsvipca)+ Istim); //----with PMCA----
	v+=Itotal*hode; 


	//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



	#ifdef ___REC_CURRENTS
	_inaca=xinaca/16.;_ica=xica/32.;_iks=xiks;_ikr=/*xikr*/xir;_ik1=xik1;/*_ina=xina;*/_inak=xinak;_iup=xiup;_svipca=xsvipca/32.;//_itof=xito;
	#endif
	return Itotal;
}
//---------------  PMCA-------------------------
double CCell::comp_svipca(void)
{//ICaB + Ipump
  //double xsvipca=svpmca*1.15*ci/(3.6+ci);
  //double xsvipca=svpmca*0.18*1.15*cs/(3.6+cs); 
  double ecan = (1.0/2/frt)*log(cao/ci);
  double xsvipca=(0.000367*(v-ecan)+svpmca*ci*ci/(0.5*0.5+ci*ci) );
  return xsvipca;
}
//-----------	sodium current following Hund-Rudy -------------------
double CCell::comp_ina(void)
{
	double ena = (1.0/frt)*log(xnao/xnai);
	double am;
	if (fabs(v+47.13)<0.001/0.1)
		am=3.2;
	else
		am = 0.32*(v+47.13)/(1.0-exp(-0.1*(v+47.13)));
	double bm = 0.08*exp(-v/11.0);
	double ah,bh,aj,bj;
	if(v<(-40.0))
	{
		 ah=0.135*exp((80.0+v)/(-6.8));
		 bh=3.56*exp(0.079*v)+310000.0*exp(0.35*v);
		 aj=((-127140.0*exp(0.2444*v)-0.00003474*exp(-0.04391*v))*(v+37.78))/(1.0+exp(0.311*(v+79.23)));
		 bj=(0.1212*exp(-0.01052*v))/(1.0+exp(-0.1378*(v+40.14)));
	}
	else
	{
		 ah=0.0;
		 bh=1.0/(0.13*(1.0+exp((v+10.66)/(-11.1))));
		 aj=0.0;
		 bj=(0.3*exp(-0.0000002535*v))/(1.0+exp(-0.1*(v+32.0)));
	}
	
	double tauh=1.0/(ah+bh);
	double tauj=1.0/(aj+bj)*jparam;
	double taum=1.0/(am+bm);
	double xina= svna*gna*xh*xj*xm*xm*xm*(v-ena);
	double inab = 0.0026*(v-ena);/*NaB*/
#ifdef REC_CURRENTS
	_ina=xina;
	_inab = inab;
#endif      

	xh = ah/(ah+bh)-((ah/(ah+bh))-xh)*exp(-hode/tauh);
	xj = aj/(aj+bj)-((aj/(aj+bj))-xj)*exp(-hode/tauj);
	xm = am/(am+bm)-((am/(am+bm))-xm)*exp(-hode/taum);
	return xina+inab;
}


//-------------- Ikr following Shannon------------------ 
double CCell::comp_ikr(void)
{
  /*
        #ifdef ___USE_VAR_FOR_CONST
	double ek = (1.0/frt)*log(xko/xki);// K reversal potential
	#endif
	const double gss=sqrt(xko/5.4);
	double xkrv1;
	if (fabs(v+7.0)<0.001/0.123)
		xkrv1=0.00138/0.123;
	else
		xkrv1=0.00138*(v+7.0)/( 1.-exp(-0.123*(v+7.0)));
	double xkrv2;
	if (fabs(v+10.0)<0.001/0.145)
		xkrv2=0.00061/0.145;
	else
		xkrv2=0.00061*(v+10.0)/(exp( 0.145*(v+10.0))-1.0);
	double taukr=1.0/(xkrv1+xkrv2);
	double xkrinf=1.0/(1.0+exp(-(v+50.0)/7.5));
	double rg=1.0/(1.0+exp((v+33.0)/22.4));
	double xikr=svkr*gkr*gss*xr*rg*(v-ek);
	xr=xkrinf-(xkrinf-xr)*exp(-hode/taukr);
  */
  double ekna = (1.0/frt)*log( (0.98*xko+0.02*xnao)/(0.98*xki+0.02*xnai) );
	
	return svkr*0.005*gkr*(v-ekna);
	//return xikr;
	
}

/*
// ----- Iks modified from Shannon, with new Ca dependence------------
double CCell::comp_iks(void)
{
	const double prnak=0.018330;
	double eks=(1.0/frt)*log((xko+prnak*xnao)/(xki+prnak*xnai));
	double xs1ss=1.0/(1.0+exp(-(v-1.50)/16.70));
	double xs2ss=xs1ss;
	double tauxs1;
	if (fabs(v+30.0)<0.001/0.0687)
		tauxs1=1/(0.0000719/0.148+0.000131/0.0687);
	else
		tauxs1=1.0/(0.0000719*(v+30.0)/(1.0-exp(-0.148*(v+30.0)))+0.000131*(v+30.0)/(exp(0.0687*(v+30.0))-1.0));
	double tauxs2=4*tauxs1;
	double gksx=0.433*(1+0.8/(1+pow((0.5/ci),3)));
	double xiks=svks*gks*gksx*xs1*xs2*(v-eks);
	xs1=xs1ss-(xs1ss-xs1)*exp(double(-hode/tauxs1));
	xs2=xs2ss-(xs2ss-xs2)*exp(double(-hode/tauxs2));
	return xiks;
}


double CCell::comp_ikr(void)
{
	#ifdef ___USE_VAR_FOR_CONST
	double ek = (1.0/frt)*log(xko/xki);// K reversal potential
	#endif
	return 0;
}

*/
double CCell::comp_iks(void)
{
#ifdef ___USE_VAR_FOR_CONST
  double ek = (1.0/frt)*log(xko/xki);// K reversal potential
#endif
  double r1 = -0.12*(v+26.5);
  double r2 = -0.038*(v+26.5);
  double an = 0.00000481333*(v+26.5)/(1-exp(r1) );
  double bn = 0.0000953333*exp(r2);

  double taun=1/(an+bn);
  double n_inf=an*taun;
  double xiks = svks*gks*nks*nks*(v-ek);
  nks = n_inf-(n_inf-nks)*exp(-hode/taun);
  return 0;//xiks;
}

//------Ik1 following Luo-Rudy formulation (from Shannon model) ------
double CCell::comp_ik1(void)
{
	#ifdef ___USE_VAR_FOR_CONST
	double ek = (1.0/frt)*log(xko/xki);// K reversal potential
	#endif
	/*const double gki=(sqrt(xko/5.4));
	double aki=1.02/(1.0+exp(0.2385*(v-ek-59.215)));
	double bki=(0.49124*exp(0.08032*(v-ek+5.476))+exp(0.061750*(v-ek-594.31)))/(1.0+exp(-0.5143*(v-ek+4.753)));
	double xkin=aki/(aki+bki);
	double xik1=svk1*gkix*gki*xkin*(v-ek);*/
	double xik1=0.2938*xko/(xko+0.210)*(v-ek)/(1+exp(0.0896*(v-ek) ) );
	return xik1;
}
//------- Ito slow following Shannon et. al. 2005 -----------
//------- Ito fast following Shannon et. al. 2005 -----------
double CCell::comp_ito(void)
{
	#ifdef ___USE_VAR_FOR_CONST
	double ek = (1.0/frt)*log(xko/xki);// K reversal potential
	#endif
	double rt1=-(v+22.5)/7.7;
	double rt2=(v+45.2)/5.7;
	double atos_inf=1.0/(1.0+exp(rt1));
	double itos_inf=1.0/(1.0+exp(rt2));
	double tas=0.493*exp(-0.0629*v)+2.058;
	double tis=270+1050/(1+exp(rt2));
	double xitos=svtos*gtos*atos*itos*(v-ek);// ito slow
	atos = atos_inf-(atos_inf-atos)*exp(-hode/tas);
	itos = itos_inf-(itos_inf-itos)*exp(-hode/tis);

	double rt3=0.03577*(v+30.0+15.);  //+15
	double rt4=-0.06237*(v+30.0+15.);   //+15
	double rt5=(v+13.5)/7.0;
	double rt6=(v+33.5)/7.0;
	double aatof = 0.18264*exp(rt3);
	double batof = 0.3956*exp(rt4);
	double aitof = 10.*0.000152*exp(-1*rt5)/(0.067083*exp(-1*rt6) + 1);
	double bitof = 10.*0.00095*exp(rt6)/(0.051335*exp(rt6) + 1);
	double taf=1/(aatof+batof);
	double tif=1/(aitof+bitof);
	double atof_inf=aatof*taf;
	double itof_inf=aitof*tif*0.63+0.37;
	double xitof=svtof*0.8*1.3*gtof*atof*atof*atof*itof*(v-ek);// ito fast
	atof = atof_inf-(atof_inf-atof)*exp(-hode/taf);
	itof = itof_inf-(itof_inf-itof)*exp(-hode/tif);
	
	double taur= 0.493*exp(-0.0629*v)+2.058;
	double tiur= 1200.0-170.0/(1.0+exp(rt2) );
	double xiur = svkur*0.8*gkur*aur*iur*(v-ek);
	aur = atos_inf-(atos_inf-aur)*exp(-hode/taur);
	iur = itos_inf-(itos_inf-iur)*exp(-hode/tiur);

	double tass= 39.3*exp(-0.0862*v)+13.17;
	double xiss = svkss*gkss*0.5*ass*iss*(v-ek);
	ass = atos_inf-(atos_inf-ass)*exp(-hode/tass);

	


	#ifdef ___REC_CURRENTS
	//_itof=xitof;_itos=xitos;!!!!
	_itof=xitof;
	_itos=xitos;
	_ikss=xiss;
	_ikur=xiur;
	//cerr << xiur << endl;
	#endif
	return xitos+xitof+xiur+xiss;
}
// -------Inak (sodium-potassium exchanger) following Shannon --------------
double CCell::comp_inak(void)
{
  const double xkmko=1.5;//1.5;	//these are Inak constants adjusted to fit
							//the experimentally measured dynamic restitution curve //Relevent
  const double xkmnai=21.0;//12.0;
	const double sigma = 0.424706;//(exp(xnao/67.3)-1.0)/7.0;
	double fnak = 1.0/(1+0.1245*exp(-0.1*v*frt)+0.0365*sigma*exp(-v*frt));
	double xinak = svnak*gnak*1.95*fnak*(1./(1.+pow(xkmnai/xnai,1.5) ))*xko/(xko+xkmko);
	return xinak;
}
// --- Inaca (sodium-calcium exchange) following Shannon and Hund-Rudy------
//	Note: all concentrations are in mM
double CCell::comp_inaca(double csm)
{
  
  double eta = 0.35;
  double xmnao = 87.5*87.5*87.5;
  double xmcao = 1380;
  double nao3 = xnao*xnao*xnao;
  double nai3 = xnai*xnai*xnai;
  double ex1 = exp( (eta-1)*v*frt );
  double ex2 = exp( eta*v*frt );
  double xinacaq = svncx*gnaca*(ex2*nai3*cao*1000-ex1*nao3*ci)/(xmnao+nao3)/(xmcao+cao*1000)/(1+0.1*ex1)/8.;
  /*
	double zw3=pow(xnai,3)*cao*exp(v*0.35*frt)-pow(xnao,3)*csm*exp(v*(0.35-1.)*frt);
	double zw4=1.0+0.2*exp(v*(0.35-1.0)*frt);
	const double xkdna=0.3;// micro M
	double aloss=1.0/(1.0+pow((xkdna/cs),3));
	const double xmcao=1.3;
	const double xmnao=87.5;
	const double xmnai=12.3;
	const double xmcai=0.0036;
	double yz1=xmcao*pow(xnai,3)+pow(xmnao,3)*csm;
	double yz2=pow(xmnai,3)*cao*(1.0+csm/xmcai);
	double yz3=xmcai*pow(xnao,3)*(1.0+pow((xnai/xmnai),3));
	double yz4=pow(xnai,3)*cao+pow(xnao,3)*csm;
	double zw8=yz1+yz2+yz3+yz4;
	double xinacaq=svncx*gnaca*aloss*zw3/(zw4*zw8);*/
	return xinacaq;
}
//	compute driving force
double CCell::comp_rxa(double csm)
{
	const double pca=0.00054;
	double za=v*2.0*frt;
	double factor1=4.0*pca*xf*xf/(xxr*temp);
	double factor=v*factor1;
	double rxa;
	if(fabs(za)<0.001)
	{
		rxa=factor1*(csm*exp(za)-0.341*(cao))/(2.0*frt);
	}
	else
	{
		rxa=factor*(csm*exp(za)-0.341*(cao))/(exp(za)-1.0);
	}
	return rxa*svca*0.7*2.7;
}
// ------ Markovian Ca current --------------------------------
//	Markov model:All parameters have been fitted directly to 
//	experimental current traces using a multidimensional current fitting
//	routine.	
double CCell::comp_icalpo(void)
{
  /*
	const double vth=0.0;
	const double s6=8.0;

	const double taupo=1.0;
	double poinf=1.0/(1.0+exp(-(v-vth)/s6));
	
	double alpha=poinf/taupo;
	double beta=(1.0-poinf)/taupo;

	const double r1=0.30;
	const double r2=3.0;

	const double cat=3.0/ccm;
	double fca=1.0/(1.0+pow(double(cat/cp),3));

	double s1=0.0182688*fca*svca;
	
	const double s1t=0.00195;
	 
	double xk1=0.024168*fca*svca;
	const double xk2=1.03615e-4;

	const double xk1t=0.00413;
	const double xk2t=0.00224;

	double s2=s1*(r1/r2)*(xk2/xk1);
	const double s2t=s1t*(r1/r2)*(xk2t/xk1t);

	const double vx=-40;
	const double sx=3.0;
	double poi=1.0/(1.0+exp(-(v-vx)/sx));
	const double tau3=3.0;
	 
	double xk3=(1.0-poi)/tau3;
	double xk3t=xk3;
			
	const double vy=-40.0;
	const double sy=4.0;
	double prv=1.0-1.0/(1.0+exp(-(v-vy)/sy));

	double recov=10.0+4954.0*exp(v/15.6);

	const double tca=78.0329;
	const double cpt=6.09365/ccm;
	double tau_ca=tca/(1.0+pow((cp/cpt),4));

#ifdef ___FORTHREED
	double tauca=(recov-tau_ca)*prv+tau_ca+1;
#else
	double tauca=(recov-tau_ca)*prv+tau_ca;
#endif
	double tauba=(recov-450.0)*prv+450.0;

	const double vyr=-40.0;
	const double syr=11.32;
	double poix=1.0/(1.0+exp(-(v-vyr)/syr));

	double xk6=fca*poix/tauca;
	double xk5=(1.0-poix)/tauca;
		 
	double xk6t=poix/tauba;
	double xk5t=(1.0-poix)/tauba;

	double xk4=xk3*(alpha/beta)*(xk1/xk2)*(xk5/xk6);
	double xk4t=xk3t*(alpha/beta)*(xk1t/xk2t)*(xk5t/xk6t);

	double po=1.0-xi1ca-xi2ca-xi1ba-xi2ba-c1-c2;

	double dc2= beta*c1+xk5*xi2ca+xk5t*xi2ba-(xk6+xk6t+alpha)*c2;
	double dc1=alpha*c2+xk2*xi1ca+xk2t*xi1ba+r2*po-(beta+r1+xk1t+xk1)*c1;

	double dxi1ca=xk1*c1+xk4*xi2ca+s1*po-(xk3+xk2+s2)*xi1ca;
	double dxi2ca=xk3*xi1ca+xk6*c2-(xk5+xk4)*xi2ca;

	double dxi1ba=xk1t*c1+xk4t*xi2ba+s1t*po-(xk3t+xk2t+s2t)*xi1ba;
	double dxi2ba=xk3t*xi1ba+xk6t*c2-(xk5t+xk4t)*xi2ba;

	c1+=dc1*hode;
	c2+=dc2*hode;
	xi1ca+=dxi1ca*hode;
	xi1ba+=dxi1ba*hode;
	xi2ca+=dxi2ca*hode;
	xi2ba+=dxi2ba*hode;
	//jrel=po;
	return po;
  */
  

  
double po;
  double dv5 = -12.6+8;
  double dvk = 6.3;

  double fv5 = -22.8;
  double fvk = 6.1;
  
  //#ifdef iso2
  //dv5 = 0;//-12.6;
  //fv5 = -28;//-36; -23; -28
  //fvk = 8.5;
  //#endif

  double dinf = 1.0/(1.0+exp(-(v-dv5)/dvk));
  double taudin = (0.035*(v-dv5))/dinf/(1.0-exp(-(v-dv5)/dvk));
  if( v > dv5-0.0001 && v < dv5+0.0001)
    taudin = 0.035*dvk/dinf;
  if( v < -80)
    dinf = 0.00000001;
  double finf = 1.-1.0/(1.0+exp(-(v-fv5)/fvk));///(1.+exp((v-60)/12.));
  double ptau = 0.0337*(v+10.5);
  double taufin = (0.02-0.007*exp(-ptau*ptau));


  c1+=(dinf-c1)*taudin*hode;
  xi1ba+=(finf-xi1ba)*taufin*hode;
  //xi2ca+=(  0.06*(1-xi2ca)-0.275/(1+(5.5/cp)*(5.5/cp))*xi2ca )*hode;
  //xi2ca+=(  0.06*2*(1-xi2ca)-0.275*3.*1.5/(1+(30./cp)*(30./cp))*xi2ca )*hode;
  xi2ca+=(  0.06*2*(1-xi2ca)-0.275*3.*1.7/(1+(30./cp)*(30./cp)*(30./cp)*(30./cp))*xi2ca )*hode;
  po=0.083*c1*xi1ba*xi2ca;


	//jrel=po;
	return po;
  

}
//----- SERCA2a uptake current ------------------------------------
double CCell::comp_iup(void)
{
	const double xup=0.5;// uptake threshold
	double xiup=svserca*2.*1.3*vup/ccm*ci*ci/(ci*ci + xup*xup);		//----CHANGED
	return xiup;
}
// ---------leak from the SR--------------------------
double CCell::comp_ileak(void)
{
  
  const double gleak=1.74e-5;//0.00002069;
  return gleak*(cj*12.4-ci);//12.4 is the size difference between the SR and whole cytoplasm
  //return gleak*(cj*cj/(cj*cj+50.0/ccm*50.0/ccm))*(cj*16.667-ci);//vsr/vcell=0.06
  
  /*
  const double gleak=1.74e-5;
  return gleak*(cj-ci);*/
}
// ---------- buffer dynamics in the myoplasm -----------------------
//buffering to calmodulin and SR are instantaneous, while buffering to
//Troponin C is time dependent.These are important to have reasonable
//Ca transient.Note: we have buffering in the submembrane space and 
//the myoplasm.
double CCell::comp_inst_buffer(double c)
{/*
	const double bcal=24.0;
	const double xkcal=7.0;
	const double srmax=47.0;
	const double srkd=0.6;
	const double bmem=15.0;
	const double kmem=0.3;
	const double bsar=42.0;
	const double ksar=13.0;*/
  const double bcal=24.0/ccm;//cmrmodification
  const double xkcal=7.0/ccm;//cmrmodification
  const double srmax=47.0/ccm;//cmrmodification
  const double srkd=0.6/ccm;//cmrmodification
  const double bmem=15.0/ccm;//cmrmodification
  const double kmem=0.3/ccm;//cmrmodification
  const double bsar=42.0/ccm;//cmrmodification
  const double ksar=13.0/ccm;//cmrmodification
	double bpx=bcal*xkcal/((xkcal+c)*(xkcal+c));
	double spx=srmax*srkd/((srkd+c)*(srkd+c));
	double mempx=bmem*kmem/((kmem+c)*(kmem+c));
	double sarpx=bsar*ksar/((ksar+c)*(ksar+c));
	//_inaca=1.0/(1.0+bpx+spx+mempx+sarpx);
	//return 1.0/(1.0+bpx+spx+mempx+sarpx);
	double beta = 1.0/(1.0+(110*0.538/(c+0.538)/(c+0.538)));
	//_inaca=beta;
	return 1/(1.+bpx+spx+mempx+sarpx);
	  //return beta;
}
// --------- release-load functional dependence ----------------
double CCell::comp_Q(void)
{
	double bv=(cstar-50.)-av*cstar;
	double Qr;
	if (cjp<50)
	{
		Qr=0.0;
	}
	else if (cjp>50.0 && cjp<cstar)
	{
		Qr=cjp-50.0;
	}
	else
	{
		Qr=av*cjp+bv;
	}
	//jrel=xiup;
	return cj*Qr/cstar;
}
double CCell::comp_dir(double po, double Qr, double rxa, double dcj)
{
	const double ay=0.05;
	double sparkV=exp(-ay*(v+30))/(1.+exp(-ay*(v+30)));
	const double g=svryr*2.58079 /*!!!!!!*/ /0.2;
	double spark_rate=g*po*fabs(rxa)*sparkV;
	//jrel=spark_rate*Qr-xir*(1-taur*dcj/cj)/taur;
	return spark_rate*Qr-xir*(1-taur*dcj/cj)/taur;
}
// ----------- dyadic junction dynamics ------------------------
double CCell::comp_dcp(double po, double Qr, double rxa)
{
  const double grel=2.5808*50*5; //26841.8;// m mol/(cm C) 
  const double ax=0.05;//0.3576;
	const double gdyad=9000.0;// m mol/(cm C) 
	double ssr=exp(-ax*(v+30))/(1.0+exp(-ax*(v+30)));
	double gain=po*Qr*fabs(rxa)*ssr;
	double xirp=svryr*grel*gain;

	double xicap=po*gdyad*fabs(rxa);
	const double taups=0.5;
	//cout << xirp << " " << xicap << " " << -(cp-cs)/taups << endl;
	return xirp+xicap-(cp-cs)/taups;
}
#endif /* ___CELL */
