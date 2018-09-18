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

//icpc uc.cpp -o ucc -O3
#define ___REC_CURRENTS
#define ___USE_VAR_FOR_CONST
#define dim 6 //----------------DIMENSION-------------
#define ROUND atoi(argv[1])


#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <time.h>
using namespace std;
#include "cell.h"
#include "cell2.h"

#define devres deva + devp + devna;// + devcj + devna;// + devapd + devcj + devdci;

CCell *cell;
double bcl=250.0;
double bcl2=0.0;
int beats=100;
double vnor[400],cainor[400],simplex[dim+1][dim], vvt[400], caca[400];
double maxnor, minnor, dnor, nordcidt, norAPD, nordiascj,nanor;
double finaldca = 0;
double finaldcab = 0;
double finaldca2 = 0;
double finaldca2b = 0;
double finalcj, finalna;
double finalapd;
double dnor2;
double tol=0.0001;
double ac = 0.0;

double getdeviation(double cond[dim]);
double getdeviationfinal(double cond[dim]);
double getdevA(double cond[dim]);
double getdevP(double cond[dim]);
double Prediction(double cond[dim]);
void sort (double *input, int num);
void swap(double *aa, double *bb);

int main( int argc, char** argv)
{
  int rrdd=ROUND;
  srand(rrdd+1);
  char filename[30],filenum[8];
  double dev[dim+1],cm[dim],ref[dim],epn[dim],contr[dim];
  int index,i,j,n,kk,limit;
  double startat[dim];
  double sum,refdev,epndev,contrdev;
  double temp;
  
  //	double bank[dim][5]={{0.8,0.9,1.0,1.1,1.2},{0.5,0.8,1.2,1.5,1.8},
  //	{0.5,0.8,1.2,1.5,1.8},{0.1,0.5,1,1.5,2},{0.5,0.75,1,1.25,1.5}};
  //Ica, Serca, NCX, Ryr, Itof
  
  //sprintf(filenum,"%d",rrdd);
  sprintf(filename,"AMB5DRWF_5SP400_part%04i.txt",rrdd);
  ofstream outfile(filename,ios::out);
  //---------------Calculate the normal value---------------
  cell= new CCell;
  svca=1;svtof=1;svtos=1;svkr=1;svks=1;svk1=1;svncx=1;svryr=1;svserca=1;svpmca=1;svna=1;svnak=1;svkss=1;svkur=1;
  maxnor = 0.0; nordcidt = 0;
  int Tn=static_cast<int>(bcl*(beats+1)/cell->getdt()+0.1), bcln=static_cast<int>(bcl/cell->getdt()+0.1), durn=static_cast<int>(1/cell->getdt()+0.1);
  double dcitemp, APDs, APDe;
  for (int tn=0;tn<Tn;tn++)
    {
      double t=tn*cell->getdt();
      if (tn%10==0 && t>=bcl*beats) 
	{
	  index=int(tn-bcln*beats)/10;
	  vnor[index]=cell->v;
	  cainor[index]=cell->ci;
	  ac += cell->ci;
	  if(cell->ci > maxnor)
	    maxnor = cell->ci;
	  if(index > 0)
	    {
	      //find action potential starting point and ending point
	      if(vnor[index - 1] < -60 && vnor[index] > -60)
		APDs = t-bcl*beats;
	      if(vnor[index - 1] > -60 && vnor[index] < -60)
		APDe = t-bcl*beats;
	      dcitemp = cainor[index] - cainor[index - 1];
	      if(dcitemp > nordcidt)
		nordcidt = dcitemp;
	    }
	  //outfile<<t-bcl*beats<<"\t"<<cell->ci<<"\t"<<cell->v<<"\t"<<cell->_svipca/16<<"\n";
	}
      if (tn%bcln < durn)
	cell->Pace(50.0);
      else
	cell->Pace();
    }
  norAPD = APDe - APDs;
  nordiascj = cell->cj;
  minnor = cell->ci;
  dnor = (maxnor - minnor);
  nanor = cell->xnai;
  delete cell;


  cell= new CCell;
  Tn=static_cast<int>(bcl2*(beats+1)/cell->getdt()+0.1), bcln=static_cast<int>(bcl2/cell->getdt()+0.1), durn=static_cast<int>(1/cell->getdt()+0.1);
  for (int tn=0;tn<Tn;tn++)
    {
      double t=tn*cell->getdt();
      if (tn%10==0 && t>=bcl2*beats) 
	{
	  //ac += cell->ci;
	  if(cell->ci > maxnor)
	    maxnor = cell->ci;
	}
      if (tn%bcln < durn)
	cell->Pace(50.0);
      else
	cell->Pace();
    }
  minnor = cell->ci;
  dnor2 = (maxnor - minnor)/3;
  delete cell;


  //trials start now
  //bcl = 102;

  for (int kk1=0; kk1<100; kk1++)
    {
      cout<<kk1<<"\n";
      /*
      temp = log(0.6)/log(2.0) + (log(2.4) - log(0.6))/log(2.0)* (rand()%1001)/1000;		//--- Ica
      *(startat+0) = pow(2.0, temp);
      temp = log(0.6)/log(2.0) + (log(1.9) - log(0.6))/log(2.0) * (rand()%1001)/1000;		//--- SERCA
      *(startat+1) = pow(2.0, temp);
      temp = log(0.1)/log(2.0) + (log(3.0) - log(0.1))/log(2.0) * (rand()%1001)/1000;		//--- NCX
      *(startat+2) = pow(2.0, temp);
      temp = log(0.05)/log(2.0) + (log(4.45) - log(0.05))/log(2.0) * (rand()%1001)/1000;		//--- RyR
      *(startat+3) = pow(2.0, temp);
      //		temp = log(0.1)/log(2.0) + (log(2.0) - log(0.1))/log(2.0) * (rand()%1001)/1000;		//--- Itof*/

      for (i=0;i<dim;++i)
	*(startat+i) = 0.2 + (2.5 - 0.2) * (rand()%1001)/1000;
      
      //------------Initialize simplex---------------------
      for (i=0;i<=dim;i++)
	{
	  for (j=0;j<dim;j++) *(*(simplex+i)+j)=*(startat+j);
	  if (i!=0)	*(*(simplex+i)+i-1)+=0.2;
	}
      for (i=0;i<dim;i++) outfile<<*(*(simplex)+i)<<"\t";
      outfile<<"\t";
      //dev[0]=getdeviation(simplex[0]);
      //cout<<dev[0]<<"\n";
      for (i=0;i<=dim;i++)	dev[i]=getdeviation(simplex[i]);
      //cout<<"end";
      limit=1;
    again:	for(int trail=0;trail<20;trail++)
	{
	  //for (i=0;i<=dim;i++) cout<<simplex[i][0]<<"\t"<<simplex[i][1]<<"\t"<<simplex[i][2]<<"\t"<<simplex[i][3]<<"\t"<<dev[i]<<"\n";
	  
	  cout<<limit;
	  sort(dev,dim+1);
	  //cout<<dev[0]<<"\t"<<getdeviation(simplex[0])<<"\t"<<Prediction(simplex[0])<<"\n";
	  
	  refdev=0;epndev=0;contrdev=0;
	  //for (i=0;i<=dim;i++) cout<<simplex[i][0]<<"\t"<<simplex[i][1]<<"\t"<<simplex[i][2]<<"\t"<<simplex[i][3]<<"\t"<<dev[i]<<"\n";
	  for (j=0;j<=dim-1;j++)
	    {
	      sum=0.0;
	      for (i=0;i<=dim;i++) sum+=*(*(simplex+i)+j);
	      cm[j]=sum/(dim+1);
	      ref[j]=cm[j]+(cm[j]-*(*(simplex+dim)+j));
			}
	  refdev=getdeviation(ref);
	  if ((refdev<=dev[dim])&&(refdev>=dev[0]))
	    {
	      for (j=0;j<=dim-1;j++) *(*(simplex+dim)+j)=ref[j];
	      dev[dim]=refdev;
	      //cout<<"ref"<<"\t"<<refdev<<"\n";
	    }
	  else if (refdev<dev[0])
	    {
	      for (j=0;j<=dim-1;j++) epn[j]=cm[j]+2*(cm[j]-*(*(simplex+dim)+j));
	      epndev=getdeviation(epn);
	      if (epndev<=refdev) {for (j=0;j<=dim-1;j++) *(*(simplex+dim)+j)=epn[j];dev[dim]=epndev;}
	      else {for (j=0;j<=dim-1;j++) *(*(simplex+dim)+j)=ref[j];dev[dim]=refdev;}
	      //cout<<"exp"<<"\t"<<epndev<<"\n";
	    }
	  else
	    {
	      
	      for (j=0;j<=dim-1;j++) contr[j]=cm[j]+0.5*(cm[j]-*(*(simplex+dim)+j));
	      contrdev=getdeviation(contr);
	      //cout<<"contr"<<"\t"<<contrdev<<"\n";
	      if (contrdev<dev[dim]) {
		for (j=0;j<=dim-1;j++) 
		  *(*(simplex+dim)+j)=contr[j];dev[dim]=contrdev;
	      }
	      else
		{
		  for (i=1;i<=dim;i++)
		    {
		      for (j=0;j<=dim-1;j++) 
			*(*(simplex+i)+j)=*(*simplex+j)+(*(*simplex+j)-*(*(simplex+i)+j))*0.5;
		      dev[i]=getdeviation(simplex[i]);
		    }
		}
	    }
	}
      if (dev[0]>tol && limit<=10)
	{
	  limit++;
	  for (i=1;i<=dim;i++)
	    {
	      for (j=0;j<dim;j++) 
		*(*(simplex+i)+j)=*(*(simplex+0)+j);
	      *(*(simplex+i)+i-1)*=(1 + 0.2 * pow(-1.0, i+limit));
	      dev[i]=getdeviation(simplex[i]);
	    }
	  goto again;
	}
      for (j=0;j<dim;j++) outfile<<simplex[0][j]<<"\t";
      //outfile<<getdevA(simplex[0])<<"\t"<<getdevP(simplex[0])<<"\t"<<dev[0]<<"\n";
      outfile<<getdeviationfinal(simplex[0])<<"\t";
      outfile<<finaldca<<"\t"<<finaldcab<<"\t"<<finaldca2<<"\t"<<finaldca2b<<"\t"<<finalna<<"\t"<<finalcj<<"\t"<<finalapd<<"\n";
	//<<"\t"<<Prediction(simplex[0])<<"\n";
    }
  outfile.close();
  return 0;
}



double getdeviation(double cond[dim])
{
  cell= new CCell;
  int Tn=static_cast<int>(bcl*(beats+1)/cell->getdt()+0.1), bcln=static_cast<int>(bcl/cell->getdt()+0.1), durn=static_cast<int>(1/cell->getdt()+0.1);
  double dev, deva, devp, devcj, devdci, devapd, time, devna, na;
  int index;
  double dcitemp, APD, diascj, dcidt, APDs, APDe, minca, maxca, dca, sumca;
  svca=cond[0];svserca=cond[1];svncx=cond[2];svryr=cond[3];svtof=cond[4];svkur=cond[5];
  //svk1=cond[5];svkur=cond[6];svkss=cond[7];
  sumca = 0; maxca = 0; dcidt = 0;
  for (int tg=0;tg<Tn;tg++)
    {
      time=tg*cell->getdt();
      if (tg%10==0 && time>=bcl*beats)
	{
	  index=int(tg-bcln*beats)/10;
	  vvt[index]=cell->v;
	  caca[index]=cell->ci;
	  sumca += cell->ci;
	  //outfile<<t-bcl*beats<<"\t"<<cell->ci<<"\t"<<cell->_inaca/8*(-1)<<"\t"<<cell->_svipca/16<<"\n";
	  if(cell->ci > maxca)
	    maxca = cell->ci;
	  if(index > 0)
	    {
	      if(vvt[index - 1] < -60 && vvt[index] > -60)
		APDs = time-bcl*beats;
	      if(vvt[index - 1] > -60 && vvt[index] < -60)
		APDe = time-bcl*beats;
	      dcitemp = caca[index] - caca[index - 1];
	      if(dcitemp > dcidt)
		dcidt = dcitemp;
	    }
	}
      if (tg%bcln < durn)
	cell->Pace(50.0);
      else
	cell->Pace();
    }
  diascj = cell->cj;
  APD = APDe - APDs;
  minca = cell->ci;
  na = cell->xnai;
  
  devapd = pow((APD - norAPD)/norAPD, 2);
  devcj = pow((diascj - nordiascj)/diascj, 2);
  devdci = pow((dcidt - nordcidt)/dcidt, 2);
  deva = pow((sumca - ac)/ac, 2);
  dca = maxca - minca;
  devp = pow((dca - dnor)/dnor, 2);
  devna = pow((na - nanor)/nanor, 2);
  delete cell;

  cell= new CCell;
  Tn=static_cast<int>(bcl2*(beats+1)/cell->getdt()+0.1), bcln=static_cast<int>(bcl2/cell->getdt()+0.1), durn=static_cast<int>(1/cell->getdt()+0.1);
  maxca = 0;
  for (int tn=0;tn<Tn;tn++)
    {
      double t=tn*cell->getdt();
      if (tn%10==0 && t>=bcl2*beats) 
	{
	  //ac += cell->ci;
	  if(cell->ci > maxca)
	    maxca = cell->ci;
	}
      if (tn%bcln < durn)
	cell->Pace(50.0);
      else
	cell->Pace();
    }
  minca = cell->ci;
  double dca2 =  maxca - minca;
  double devp2 = pow((dca2-dnor2)/dnor2,2);
  delete cell;
  
  dev = devres;//deva + devp + devcj + devna;// + devapd + devcj + devdci;


  return dev;
}


double getdeviationfinal(double cond[dim])
{
  cell= new CCell;
  int Tn=static_cast<int>(bcl*(beats+2)/cell->getdt()+0.1), bcln=static_cast<int>(bcl/cell->getdt()+0.1), durn=static_cast<int>(1/cell->getdt()+0.1);
  double dev, deva, devp, devcj, devdci, devapd, time, devna, na;
  int index;
  double dcitemp, APD, diascj, dcidt, APDs, APDe, minca, maxca, dca, sumca;
  double mincab;
  svca=cond[0];svserca=cond[1];svncx=cond[2];svryr=cond[3];svtof=cond[4];svkur=cond[5];
  //svk1=cond[5];svkur=cond[6];svkss=cond[7];
  sumca = 0; maxca = 0; dcidt = 0;
  finaldcab = 0;
  finaldca2b = 0;
  for (int tg=0;tg<Tn;tg++)
    {
      time=tg*cell->getdt();
      if (tg%10==0 && time>=bcl*beats && time<bcl*(beats+1))
	{
	  index=int(tg-bcl*beats/cell->getdt())/10;
	  vvt[index]=cell->v;
	  caca[index]=cell->ci;
	  sumca += cell->ci;
	  //outfile<<t-bcl*beats<<"\t"<<cell->ci<<"\t"<<cell->_inaca/8*(-1)<<"\t"<<cell->_svipca/16<<"\n";
	  if(cell->ci > maxca)
	    maxca = cell->ci;
	  if(index > 0)
	    {
	      if(vvt[index - 1] < -60 && vvt[index] > -60)
		APDs = time-bcl*beats;
	      if(vvt[index - 1] > -60 && vvt[index] < -60)
		APDe = time-bcl*beats;
	      dcitemp = caca[index] - caca[index - 1];
	      if(dcitemp > dcidt)
		dcidt = dcitemp;
	    }
	  minca = cell->ci;
	}
      if (tg%10==0 && time>=bcl*(beats+1) )
	{
	  if(cell->ci > finaldcab)
	    finaldcab = cell->ci;
	  mincab = cell->ci;
	}
      if (tg%bcln < durn)
	cell->Pace(50.0);
      else
	cell->Pace();
    }
  finaldcab-=mincab;

  APD = APDe - APDs;
  na = cell->xnai;
  finalna = na;
  diascj = cell->cj;
  finalcj = diascj;
  
  devapd = (APD)/norAPD;
  finalapd = devapd;
  devcj = pow((diascj - nordiascj)/diascj, 2);
  devdci = pow((dcidt - nordcidt)/dcidt, 2);
  deva = pow((sumca - ac)/ac, 2);
  dca = maxca - minca;
  finaldca = dca;
  devp = pow((dca - dnor)/dnor, 2);
  devna = pow((na - nanor)/nanor, 2);
  delete cell;

  cell= new CCell;
  maxca = 0;
  Tn=static_cast<int>(bcl2*(beats+2)/cell->getdt()+0.1), bcln=static_cast<int>(bcl2/cell->getdt()+0.1), durn=static_cast<int>(1/cell->getdt()+0.1);
  cout << endl;
  for (int tn=0;tn<Tn;tn++)
    {
      double t=tn*cell->getdt();
      if (tn%10==0 && t>=bcl2*beats && t<bcl2*(beats+1)) 
	{
	  if(cell->ci > maxca)
	    maxca = cell->ci;
	  minca = cell->ci;
	}
      if (tn%10==0 && t>=bcl2*(beats+1) )
	{
	  if(cell->ci > finaldca2b)
	    finaldca2b = cell->ci;
	  mincab = cell->ci;
	}
	
      if (tn%bcln < durn)
	cell->Pace(50.0);
      else
	cell->Pace();
    }
  finaldca2b-=mincab;
  double dca2 =  maxca - minca;
  finaldca2 = dca2;
  double devp2 = pow((dca2-dnor2)/dnor2,2);
  delete cell;
  
  dev = devres;//deva + devp + devcj + devna;// + devapd + devcj + devdci;


  return dev;
}




double Prediction(double cond[dim]){
  double fa_ca, fa_ser, fa_ncx, fa_ryr, fa_tof;
  double fp_ca, fp_ser, fp_ncx, fp_ryr, fp_tof;
  double predict;
  
  fa_ca = 1.225*(cond[0] - 1) + 0.47383 * pow(cond[0] - 1, 2);
  fp_ca = 1.8195*(cond[0] - 1) + 0.93463 * pow(cond[0] -1, 2);
  
  fa_tof = 0.011546*(cond[4] - 1) - 0.76475*pow(cond[4] - 1, 2);
  fp_tof = 0.16047*(cond[4] - 1) - 1.2947*pow(cond[4] - 1, 2);
  
  fa_ser = -0.16924*(cond[1] -1) - 0.028285*pow(cond[1]-1, 2);
  fp_ser = 0.33497*(cond[1] -1) - 0.40756*pow(cond[1]-1, 2) + 0.11344*pow(cond[1]-1, 3);

  fa_ncx = 1.0905*exp(-1.8365*cond[2]) - 0.17345;
  fp_ncx = 1.6258*exp(-1.9848*cond[2]) - 0.22354;

  fa_ryr = -0.30844*exp(-12.4785*cond[3]) + 0.093432*exp(-0.25693*cond[3]) - 0.07396;
  fp_ryr = -0.57065*exp(-3.4085*cond[3]) + 1.661*exp(-0.027482*cond[3]) - 1.5997;

  predict = pow(fa_ca+fa_tof+fa_ser+fa_ncx+fa_ryr, 2) + pow(fp_ca+fp_tof+fp_ser+fp_ncx+fp_ryr, 2);
  return predict;
}

void sort(double *input, int num)
{
  int i,j,k;
  double t;
  for (i=0;i<num-1;i++)
    {
      k=i;
      for (j=i+1;j<num;j++) 
	if (*(input+j)<*(input+k)) k=j;
      t=*(input+k);*(input+k)=*(input+i);*(input+i)=t;
      for (j=0;j<=dim-1;j++) swap(*(simplex+i)+j,*(simplex+k)+j);
    }
}
void swap(double *aa,double *bb)
{
  double temp=*aa;
  *aa=*bb;
  *bb=temp;
}
/*
double getdevA(double cond[dim])
{
  cell= new CCell;
  int Tn=bcl*(beats+1)/cell->getdt(), bcln=bcl/cell->getdt(), durn=1/cell->getdt();
  double deva, sumca, time;
  int index;
  svca=cond[0];svserca=cond[1];svncx=cond[2];svryr=cond[3];svtof=cond[4];
  sumca = 0;
  for (int tn=0;tn<Tn;tn++)
    {
      time=tn*cell->getdt();
      if (tn%10==0 && time>=bcl*beats)
	{
	  sumca += cell->ci;
	}
      if (tn%bcln < durn)
	cell->Pace(50.0);
      else
	cell->Pace();
    }
  deva = (sumca - ac)/ac;
  delete cell;
  return deva;
}

double getdevP(double cond[dim])
{
  cell= new CCell;
  int Tn=bcl*(beats+1)/cell->getdt(), bcln=bcl/cell->getdt(), durn=1/cell->getdt();
  double devp, minca, maxca, dca, time;
  int index;
  svca=cond[0];svserca=cond[1];svncx=cond[2];svryr=cond[3];svtof=cond[4];
  maxca = 0;
  for (int tn=0;tn<Tn;tn++)
    {
      time=tn*cell->getdt();
      if (tn%10==0 && time>=bcl*beats)
	{
	  index=int(tn-bcl*beats/cell->getdt())/10;
	  if(cell->ci > maxca)
	    maxca = cell->ci;
	}
      if (tn%bcln < durn)
	cell->Pace(50.0);
      else
	cell->Pace();
    }
  minca = cell->ci;
  dca = maxca - minca;
  devp = (dca - dnor)/dnor;
  delete cell;
  return devp;
}
*/
