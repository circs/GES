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
*Title: The Ca2+ transient as a feedback sensor controlling ionic conductance variability and compensation in mouse cardiomyocytes
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


#ifndef ___CELL_H
#define ___CELL_H

/* ---------------- UCLA Model ver 1.00 ---------------- **
*
* Contact Information
*
* Departments of Medicine (Cardiology)
* David Geffen School of Medicine at UCLA
*
* Daisuke Sato		 dasato (at) mednet.ucla.edu
* Yohannes Shiferaw  yshiferaw (at) csun.edu
* James N Weiss 	 JWeiss (at) mednet.ucla.edu
*
** ---------------- ------------------- ---------------- */


// #define ___REC_CURRENTS //record currents (more memory)
// #define ___USE_VAR_FOR_CONST //use variables for Gto Gks Gkr etc. instead of constants (more memory)

#include <iostream>
using namespace std;
#define _USE_MATH_DEFINES
#include <cmath>

double svca,svtof,svtos,svkr,svks,svk1,svncx,svryr,svserca,svpmca,svna,svnak,svkur,svkss;

class CCell{
private:
	double jparam;//tauj*jparam
	double PaceX(double stim=0);
	static const int N=31;
	static const double Vc;
	static const double stim;
	static const double stimduration;
	static const double temp;// temperature (K)
	static const double xxr;//
	static const double xf;// Faraday's constant
	static const double frt;

	#ifndef ___USE_VAR_FOR_CONST
	static const double xnao;//mM	 external Na
	static const double xki;// mM	 internal K
	static const double xko;// mM	 external K
	static const double cao;// mM	 external Ca
	static const double ek;

	static const double gca;// ica conductance
	static const double gtos;// ito slow conductance 
	static const double gtof;//ito fast conductance 
	static const double gnaca;// exchanger strength 
	static const double gks;
	static const double gkr; 
	static const double gkss;
	static const double gkur;
	static const double vup;// uptake strength
	static const double gna;// sodium conductance (mS/micro F) 
	static const double gkix;// Ik1 conductance
	static const double gnak;

	static const double taur;// spark lifetime (ms)
	static const double taus;// diffusional delay (ms)
	static const double taua;// NSR-JSR diffusional delay (ms)
	static const double av;
	static const double cstar;
	#endif

	double comp_svipca(void);//-------------- PMCA-------------
	double comp_ina (void);
	double comp_ikr(void);
	double comp_iks(void);
	double comp_ik1(void);
	double comp_ito(void);
	double comp_inak(void);
	double comp_inaca(double csm);
	double comp_icalpo(void);
	double comp_iup(void);
	double comp_ileak(void);
	double comp_inst_buffer(double c);

	double comp_rxa(double csm);
	double comp_Q(void);
	double comp_dir(double po, double Qr, double rxa, double dcj);
	double comp_dcp(double po, double Qr, double rxa);
	double vold;
	double hode,hpde;

public:
	
	double Vinput[10];


	double Pace(double stim=0);
	double PaceVClamp(double clampv);
	double setJparam(double newjp){jparam=newjp;return newjp;}
	double setdt(double dtt){hpde=dtt;return hpde;}
	double getdt(void){return hpde;}
	int getDim(void){return N;}
	double getVc(void){return Vc;}
	double getstim(void){return stim;}
	double getstimduration(void){return stimduration;}
	void ClampAP(double t, double BCL, double APD=0);//BCL ms
	CCell(void);
	virtual ~CCell();
	CCell& operator=(const CCell& cell);
	void Prepare(double BCL=300, int Iter=0);
	double *y;
	double &xm,&xh,&xj,&xr,&nks,/*&xs1,&xs2,*/&atos,&itos,&v,&ci,&cs,&cj,&cjp,&cp, &aur, &iur, &ass, &iss;
	double &xir,&c1,&c2,&xi1ca,&xi1ba,&xi2ca,&xi2ba,&xnai,&atof,&itof,&tropi,&trops,&jrel;

	#ifdef ___USE_VAR_FOR_CONST
	double gca;//ica conductance
	double gtos;// ito slow conductance 
	double gtof;// ito fast conductance 
	double gnaca;// exchanger strength 
	double gks;
	double gkr; 
	double gkss;
	double gkur;
	double vup;
	double gna;// sodium conductance (mS/micro F) 
	double gkix;// Ik1 conductance
	double gnak;

	double xnao;//mM external Na
	double xki;//mM internal K
	double xko;//mM external K
	double cao;//mM external Ca

	double taus;//	diffusional delay (ms)
	double taur;// spark lifetime (ms)
	double taua;// NSR-JSR diffusional delay (ms)
	double av;
	double cstar;
	#endif

	#ifdef ___REC_CURRENTS
	double _inaca,_ica,_iks,_ikr,_itof,_itos,_ik1,_ina,_inak,_iup,_svipca,_ikss, _ikur, _inab;
	#endif
};

const double CCell::Vc=-80;
const double CCell::stim=80;
const double CCell::stimduration=0.5;//2;


// ---------------constant parameters ------------------------------
const double CCell::temp=298.0; //308.0;// temperature (K)
const double CCell::xxr=8.314;//***// Ideal Gas Constant
const double CCell::xf=96.485;//***// Faraday's constant 
const double CCell::frt=xf/(xxr*temp);

#ifndef ___USE_VAR_FOR_CONST
const double CCell::xnao=140.0;//136.0//mM	  external Na
const double CCell::xki=140.0;// mM	 internal K
const double CCell::xko=5.40;//***//mM	 external K
const double CCell::cao=1.8;//***/ mM	 external Ca
const double CCell::ek = (1.0/frt)*log(xko/xki);// K reversal potential

const double CCell::gca=172.9;//182// ica conductance SF^-1
const double CCell::gtos=0.0;//0.04// ito slow conductance 
const double CCell::gtof=0.4067;//0.11// ito fast conductance 
const double CCell::gnaca=0.84;//18015;//0.84 exchanger strength (uM/s)
/* gnaca = knaca *Cm/F/V, knaca = 292.8AF^-1, Cm = 1.534E-4 uF, V = 25.84E-6uL*/

const double CCell::gkr=0.078;//0.0125;// Ikr conductance (replace w/ m Ikr)
const double CCell::gks=0.0575;//0.32;// Iks conductance (replace w/ m Ik )
const double CCell::gkss=0.05;
const double CCell::gkur=0.16;
const double CCell::gkix=0;//0.3;// Ik1 conductance (unused )


const double CCell::gnak=0.88;//1.5;
const double CCell::vup=0.4;//0.3;// uptake strength
const double CCell::taus=8//4.0;// diffusional delay (ms)
const double CCell::gna=13.0;//12.0;// sodium conductance (mS/micro F) 
const double CCell::taur=30.0;// spark lifetime (ms)
const double CCell::taua=20;//100.0;// NSR-JSR diffusional delay (ms)
const double CCell::av=11.3;//realease slope
const double CCell::cstar=90.0;//Threshold for steep release function (uM/cytostol)
#endif

#endif /* ___CELL_H */
