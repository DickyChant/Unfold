/*************************************************************************
  > File Name: UnfoldAMPT.C
 ************************************************************************/


#if !defined(__CINT__) || defined(__MAKECINT__)
#include <iostream>
#include "TRandom.h"
#include "TRandom2.h"
#include "TH1.h"
#include "TProfile.h"
#include "TF1.h"
#include "TH3.h"
#include "TH2.h"
#include "TMath.h"
#include "TFile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TComplex.h"
#include "TSystem.h"
#include "/Users/dicky/codes/RooUnfold/src/RooUnfoldResponse.h"
#include "/Users/dicky/codes/RooUnfold/src/RooUnfoldBayes.h"
#include "/Users/dicky/codes/RooUnfold/src/RooUnfoldSvd.h"
#include "/Users/dicky/codes/RooUnfold/src/RooUnfoldTUnfold.h"
#endif

#define NBIN 22050
int NITERS[9]={0,1,2,4,8,16,32,64,128};
using namespace std;


const int NITR= 9;
//const int NITR= 4;
const int NHAR= 6;
//const int NHAR= 6;
const int Ncent=14;
//const int Ncent=16;
//const int Nmax=90;     //(NHAR-1)*Ncent;



double xran[NHAR][Ncent]={1.0};
double sig2[NHAR][Ncent]={0.};
int numbins[NHAR][Ncent] = {0};

double ran;
int nbins;
double sig;

char name[200];
const double pi2 = 2*TMath::Pi();
const double pi = TMath::Pi();
const double bias=pi/20;


TF1* efit;
TH1D* hfun;
TH1*hMeas,*hReco,*htrs,*htrain;

TH1* h_vnc[NHAR][Ncent];
TH1* h_qn[NHAR][Ncent];
TH1* h_dqxqy[NHAR][Ncent];
TH2* h_qndxy[NHAR][Ncent];

TH1* hUnfold[NHAR][Ncent][NITR];
TH1* hRefold[NHAR][Ncent][NITR];

TFile* fout;

const int ih=  0 ;  //warning: ih-[0,5], icent-[0,15]
const int icent= 0  ;



double bg1(double x,double p0,double p1){//BG function
    double ret=0;
    double aver=p0/p1;
    double r=x/p1;
    if(aver<10){
        ret = r*exp(-(r*r+aver*aver))*TMath::BesselI0(2*aver*r);
    }else{
        double zed = 16*aver*r;
        ret =0.282094792*exp(-(r-aver)*(r-aver))*sqrt(r/aver)*(1+1./zed+9./(2*zed*zed)+225./6/pow(zed,3)+11025./(24*pow(zed,4)));
    }
    return ret;
}
double bg2(double *x,double *p){//BG function
    double ret=0;
    double aver=p[0]/p[1];
    double r=x[0]/p[1];
    if(aver<10){
        ret = r*exp(-(r*r+aver*aver))*TMath::BesselI0(2*aver*r);
    }else{
        double zed = 16*aver*r;
        ret =0.282094792*exp(-(r-aver)*(r-aver))*sqrt(r/aver)*(1+1./zed+9./(2*zed*zed)+225./6/pow(zed,3)+11025./(24*pow(zed,4)));
    }
    return ret;
}

//double myfun2(double* x,double* par){
//    double ret=0;
//    ret = par[0]*x[0]*exp(-(par[1]*par[1]+x[0]*x[0])/(par[2]*par[2]*2))*TMath::BesselI0(par[1]*x[0]/(par[2]*par[2]));
//    return ret;
//}

int    Rebin[]={20,  50, 100};
int    NI[NITR]={0,1,2,4,8,16,32,64,128};
//int    NI[NITR]={0,1,2,4};
//double ran=0;

//void unfoldonce(TH1*hin,TH1*hin2, int nitr){
//void unfoldonce(TH1*hin,int nitr,double nbins, double ran,double sig){

void myunfold_1()
{
    float v=0.08,vs=0.035,sig=0.025;
    
    TFile *f = new TFile("myunfold3d.root","recreate");
    
    TF1*fun = new TF1("fun",bg2,0,0.3,2);
    fun->SetParameters(v,vs);
    
    TRandom*gran = new TRandom();
    
    
    double wd=vs; if(sig>wd) wd=sig;
    double ran1 = v+wd*6;
    double ran2 = v+wd*4;
    
    
    
    
    
    
    TH1D *P=new TH1D("v/q_p","v/q_p",50,0,ran1);
    TH1D *M=new TH1D("v/q_m","v/q_m",7,-2*vs,2*vs);
    TH1D *T=new TH1D("v/q_the","v/q_the",21,-pi,pi);
    
    TH1D *Uv=new TH1D("v_for_Unfold","v_for_Unfold",NBIN,0,NBIN);
    TH1D *Uq=new TH1D("q_for_Unfold","q_for_Unfold",NBIN,0,NBIN);
    
    
    for(int ir=0;ir<10000000;ir++){
        double v1 = fun->GetRandom();
        double the1=gran->Rndm()*pi2; //angle for v1
        
        double x1 = gran->Gaus(0,sig);
        double y1 = gran->Gaus(0,sig);
        double x2 = gran->Gaus(0,sig);
        double y2 = gran->Gaus(0,sig);
        
        
        
        double x3 = gran->Gaus(0,sig);
        
        TComplex V1(v1+x1,y1),V2(v1+x2,y2);
        double q1 = V1.Rho();
        TComplex V3=(V2*V2.Conjugate(V1))/q1;
        double q2x=V3.Re(),q2y=V3.Im();
        double q2 = sqrt(q2x*q2x+q2y*q2y);
        double q4 = V2.Rho();
        TComplex q1q2=V2*V2.Conjugate(V1);
        double The=q1q2.Theta();
        double the = V3.Theta();
        TComplex V4(v1,0);
        TComplex V5=(V4*V4.Conjugate(V4))/v1;
        double the2=V5.Theta();
        if(the>pi) the-=pi2; if(the<-pi) the+=pi2;
        if(the2>pi) the2-=pi2; if(the2<-pi) the2+=pi2;
        
        
        
        
        //h3[1]->Fill(q1,q2,the*(q1+q2)/2);
        //    Int_t number_dec = detectorDistribution->GetGlobalBinNumber(q1,q2,the);
        //    Int_t number_gen = bgrBinning->GetGlobalBinNumber(v1,v1,the2);
        
        double qp=(q1+q4)/2;
        double qm=(q1-q4)/2;
        int v_index=P->FindBin(v1)*21*21+M->FindBin(0.0)*21+T->FindBin(0.000000001);
        int q_index=P->FindBin(qp)*21*21+M->FindBin(qm)*21+T->FindBin(the);
        
        
        //h3[1]->Fill(q1,q2,the*(q1+q2)/2);
        //    Int_t number_dec = detectorDistribution->GetGlobalBinNumber(q1,q2,the);
        //    Int_t number_gen = bgrBinning->GetGlobalBinNumber(v1,v1,the2);
        
        Uv->Fill(v_index);
        Uq->Fill(q_index);
        
        
        
        
        
        
    }
   
    Uv->Write();
    Uq->Write();
    
    RooUnfoldResponse response (NBIN,0,NBIN);
    TH2D* histM=new TH2D("mig_matrix","mig_matrix",NBIN,0,NBIN,NBIN,0,NBIN);
    
    
    
    for(int ir=0;ir<1000000000;ir++){
        double v1 = fun->GetRandom();
        double v2 = fun->GetRandom();
        double the1=gran->Rndm()*pi2-pi; //angle for v1
        
        double x1 = gran->Gaus(0,sig);
        double y1 = gran->Gaus(0,sig);
        double x2 = gran->Gaus(0,sig);
        double y2 = gran->Gaus(0,sig);
        
        
        
        double x3 = gran->Gaus(0,sig);
        
        TComplex V1(v1+x1,y1),V2(v2*cos(the1)+x2,v2*sin(the1)+y2);
        double q1 = V1.Rho();
        TComplex V3=(V2*V2.Conjugate(V1))/q1;
        double q2x=V3.Re(),q2y=V3.Im();
        double q2 = sqrt(q2x*q2x+q2y*q2y);
        double q4 = V2.Rho();
        TComplex q1q2=V2*V2.Conjugate(V1);
        double The=q1q2.Theta();
        double the = V3.Theta();
        TComplex V4(v1,0);
        TComplex V5=(V4*V4.Conjugate(V4))/v1;
        double the2=V5.Theta();
        if(the>pi) the-=pi2; if(the<-pi) the+=pi2;
        if(the2>pi) the2-=pi2; if(the2<-pi) the2+=pi2;
        
        
        
        double vp=(v1+v2)/2;
        double vm=(v1-v2)/2;
        double qp=(q1+q4)/2;
        double qm=(q1-q4)/2;
        int v_index=P->FindBin(vp)*21*21+M->FindBin(vm)*21+T->FindBin(the1);
        int q_index=P->FindBin(qp)*21*21+M->FindBin(qm)*21+T->FindBin(the);
        
        response.Fill(q_index,v_index);
        histM->Fill(q_index,v_index);
        
        
        
        
        
    }
    response.Write();
     histM->Write();
    
    // char name[100];
    // TCanvas* canRefold = new TCanvas("unfold","hunfold", 800, 600);

    // for(int ii=0;ii<9;ii++)
    // {
    //     sprintf(name,"itr_%d",NITERS[ii]);
    //     if(NITERS[ii]==0)
    //     {
    //         Uv->Draw();
    //         Uq->Draw("same");
    //     }
    //     else{
    //         RooUnfoldBayes unfold (&response, Uq, NITERS[ii]);
            
    //         TH1D* hReco;
    //         hReco = (TH1D*) unfold.Hreco();
    //         hReco->Write();
    //         TH1D*h=(TH1D*)hReco->Clone(name);
    //         h->SetTitle(name);
    //         h->Write();
    //         h->Draw("same");
    //         h->SetMarkerStyle(20);
    //         h->SetMarkerColor(ii);
    //     }
    //     canRefold->Write();
    // }
}


