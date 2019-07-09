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
#include "TMath.h"
#include "TFile.h"
#include "TLine.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TSystem.h"
#include "/Users/dickie/codes/RooUnfold/src/RooUnfoldResponse.h"
#include "/Users/dickie/codes/RooUnfold/src/RooUnfoldBayes.h"
#include "/Users/dickie/codes/RooUnfold/src/RooUnfoldSvd.h"
#include "/Users/dickie/codes/RooUnfold/src/RooUnfoldTUnfold.h"
#endif
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

char name[200];


double myfun(double x,double p0,double p1){
    
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



double myfun2(double *x,double *par){
    
    double ret=0;
    double aver=par[0]/par[1];
    double r=x[0]/par[1];
    if(aver<10){
        ret = par[2]*r*exp(-(r*r+aver*aver))*TMath::BesselI0(2*aver*r);
    }else{
        double zed = 16*aver*r;
        ret =par[2]*0.282094792*exp(-(r-aver)*(r-aver))*sqrt(r/aver)*(1+1./zed+9./(2*zed*zed)+225./6/pow(zed,3)+11025./(24*pow(zed,4)));
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
void unfoldonce(TH1*hin,TH1*hin2,int nitr,double nbins, double ran,double sig){
//    //need response function which is hmeasure
//    //need prior c^0, which is htrain
//    TRandom*gran = new TRandom();
//    cout<<"our sig is "<<sig<<endl;
//    TH1* htmp  = hin; //obs
//    //    htrain = hin;  //true
//    htrain = hin2;  //true
//
//    if(nitr==0){
//        int b1 = htmp->FindBin(0.000001);
//        int b2 = htmp->FindBin(ran-0.00001);
//        cout<<b1<<" "<<b2<<endl;
//        hMeas = new TH1D("hMeas","hMeas",nbins,0,ran);
//        for(int bin =b1;bin<=b2;bin++){
//            if(htmp->GetBinContent(bin)!=0){
//                hMeas->SetBinContent(bin-b1+1,htmp->GetBinContent(bin));
//                hMeas->SetBinError(bin-b1+1,htmp->GetBinError(bin));
//            }
//        }
//        hMeas->SetBinContent(nbins+1,0);  hMeas->SetBinError(nbins+1,1e-7);
//        hMeas->SetBinContent(0,0);  hMeas->SetBinError(0,1e-7);
//        hReco = hMeas;
//        //        hReco=(TH1D*)hMeas->Clone("hReco");
//    }else{//does unfolding
//        //cout << "========================= TRAIN ========================" << endl;
//        RooUnfoldResponse response (nbins, 0, ran);
//
//
//    }
//    //cout << "=========================== UNFOLD =======================" << endl;
//        RooUnfoldBayes   unfold (&response, hMeas, nitr);    // OR
//        hReco= (TH1D*) unfold.Hreco();
//        //unfold.PrintTable (cout, hReco);
//    }
//    //do refold
//    htrs  = (TH1*)hReco->Clone("htrs");
//   if(nitr!=0){
//    htrs->Reset();
//    int nb = hReco->GetNbinsX();
//    for (Int_t i=1; i<nb; i++) {
//        double xt = hReco->GetBinCenter(i);
//        if(xt>ran) continue;
//
//        int  N = hReco->GetBinContent(i);
//        double from = xt-6*sig;
//        double to = xt+6*sig;
//        if(from<0) from =0;
//        hfun->Reset();
//        int bin1 = hfun->FindBin(from);
//        int bin2 = hfun->FindBin(to);
//        for(int j=bin1;j<bin2;j++){
//            hfun->SetBinContent(j,myfun(hfun->GetBinCenter(j),xt,sig));
//        }
//        for(int j=0;j<N;j++){
//            Double_t x=  hfun->GetRandom();
//            htrs->Fill(x);
//        }
//     }
//    }
}

void readin(){
    
    TFile* fin = TFile::Open("/Users/dickie/codes/UUdataUnfolding/output.root");
//    TFile* fin1 = TFile::Open("histo_singlevnscal.root");
    
    for(int ih=0; ih<NHAR-1; ih++){
        int fb=2;
        for(int icent=0; icent!=Ncent; icent++){
        sprintf(name,"h_qn_fb%d_ih%d_cent%d", fb, ih,icent);
        h_qn[ih][icent] =(TH1D*)fin->Get(name);

        sprintf(name,"h_qndxy_ih%d_cent%d",ih,icent);
        h_qndxy[ih][icent] =(TH2D*)fin->Get(name);
    
        sprintf(name,"hdeltaqnxqny_ih%d_cent%d",ih,icent);
        h_dqxqy[ih][icent] = (TH1D*)fin->Get(name);

//        sprintf(name,"VnCorScal3_w0_har%d_eta0_cent%d", ih+2,icent);
//        h_vnc[ih][icent]  =(TH1D*)fin1->Get(name);
        }
        
    }
    
    cout<<"Finish readin"<<endl;

//    fout = new TFile("./myUnfold/unfoldout.root","recreate");

}


void UnfoldbayesUU_quicklook(){

    hfun = new TH1D("hfun","hfun",10000,0,1);
    efit = new TF1("efit","pol0",0,0.5);
    efit->SetLineWidth(1);
    efit->SetLineColor(2);

    readin();
    
//    double mysig=h_dqxqy[ih][icent]->GetRMS();
    // h_vnc[0]->Draw();
    h_qn[ih][icent]->Rebin(4);
    h_dqxqy[ih][icent]->Rebin(4);

//    for(int ih=0; ih<NHAR-1; ih++){
//        for(int icent=0; icent!=Ncent; icent++){
//    h_qn[ih]->Rebin(4);
//    nbins = h_qn[ih]->FindBin(ran);
//    nbins =nbins-1;
    
    int NmaxBin = h_qn[ih][icent]->FindLastBinAbove(0,1);
    numbins[ih][icent]=NmaxBin;
    double Binwidth= h_qn[ih][icent]->GetBinWidth(NmaxBin);
    
    xran[ih][icent] = h_qn[ih][icent]->GetBinCenter(NmaxBin);
    xran[ih][icent]+= 0.5*Binwidth;

    cout<<"**************************************"<<endl;
    cout<<" nbins= "<< numbins[ih][icent]<<"  ran= "<<xran[ih][icent]<<endl;
//    cout<<" nbins= "<< numbins<<"  ran= "<<xran <<endl;
//    cout<<" Binwidth= "<< Binwidth<<endl;
    cout<<"**************************************"<<endl;
    
    double sigx = h_qndxy[ih][icent]->GetRMS(1);
    double sigy = h_qndxy[ih][icent]->GetRMS(2);
 
//    std::cout<<"sigx= "<<sigx<<" sigy= "<<sigy<<std::endl;
 
//    double sigy = 0.02102 ;
    sig2[ih][icent] = 0.5*(sigx+sigy);
    sig2[ih][icent] /= sqrt(2);
//        }
//    }
    sprintf(name,"outUnfolding_ih%d_cent%d.root",ih,icent);
    TFile* fout = new TFile(name,"recreate");
    fout->cd();
    
//    for(int ih=0; ih<NHAR-1; ih++){
//        for(int icent=0; icent!=Ncent; icent++){
           h_qn[ih][icent]->Write();
           h_dqxqy[ih][icent]->Write();
//        }
//    }
   
    
//    for(int ih=0; ih<NHAR-1; ih++){
//        for(int icent=0; icent!=Ncent; icent++){
            for(int it=0; it<NITR; it++){
//        unfoldonce(h_qn[ih], h_qn[ih],  NI[it]);
                 
            nbins = numbins[ih][icent];
            ran = xran[ih][icent];
            sig = sig2[ih][icent];
                                                        
//        unfoldonce(h_qn[ih][icent], NI[it],nbins,ran,sig);
//        unfoldonce(h_qn[ih][icent], NI[it],nbins,ran,sig);
        unfoldonce(h_qn[ih][icent],h_qn[ih][icent+1], NI[it],nbins,ran,sig);

        sprintf(name, "hUnfold_ih%d_cent%d_itr%d", ih,icent,it);
        hUnfold[ih][icent][it] = (TH1D*)hReco->Clone(name);
        hUnfold[ih][icent][it] ->SetTitle(name);
        sprintf(name, "hRefold_ih%d_cent%d_itr%d", ih,icent,it);
        hRefold[ih][icent][it] = (TH1D*)htrs->Clone(name);
        hRefold[ih][icent][it] ->SetTitle(name);
        
        hUnfold[ih][icent][it]->Write();
        hRefold[ih][icent][it]->Write();
            }
//        }
//    }
    
    TCanvas* canUnfold = new TCanvas("canUnfold","hUnfold", 800, 600);

    for(int it=0; it<NITR; it++){
        if(it==0) hUnfold[ih][icent][it]->Draw();
        else hUnfold[ih][icent][it]->Draw("same");
        hUnfold[ih][icent][it]->SetLineColor(it+1);
        hUnfold[ih][icent][it]->SetMarkerStyle(20);
        hUnfold[ih][icent][it]->SetMarkerColor(it+1);
    }

    TLegend *l = new TLegend(0.6002506,0.6191304,0.8922306,0.8730435,NULL,"brNDC");   
 // TLegend*l = new TLegend(0.6,0.7,1,0.95);
    for(int iter=0;iter<NITR;iter++){
        sprintf(name, "Iter %d",NI[iter]);
        l->AddEntry(hUnfold[ih][icent][iter],name,"pl");
    }
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(43);
    l->SetTextSize(19);
    l->Draw("same");


    TCanvas* canRefold = new TCanvas("canRefold","hRefold", 800, 600);

    for(int it=0; it<NITR; it++){
        if(it==0) hRefold[ih][icent][it]->Draw();
        else hRefold[ih][icent][it]->Draw("same");
        hRefold[ih][icent][it]->SetLineColor(it+1);
        hRefold[ih][icent][it]->SetMarkerStyle(20);
        hRefold[ih][icent][it]->SetMarkerColor(it+1);
	
    }
    
//    gStyle->SetOptStat(0);
//    gStyle->SetOptTitle(0);

    TLegend *l1 = new TLegend(0.6002506,0.6191304,0.8922306,0.8730435,NULL,"brNDC");   
 // TLegend*l = new TLegend(0.6,0.7,1,0.95);
    for(int iter=0;iter<NITR;iter++){
        sprintf(name, "Iter %d",NI[iter]);
        l1->AddEntry(hRefold[ih][icent][iter],name,"pl");
    }
    l1->SetBorderSize(0);
    l1->SetFillStyle(0);
    l1->SetTextFont(43);
    l1->SetTextSize(19);
    l1->Draw("same");

    canUnfold->Write();
    canRefold->Write();
 
    
    TH1* hres  = h_dqxqy[ih][icent]; //obs
    TH1D *hResfun= new TH1D("hResfun","",nbins,0,ran);
    int b1 = hres->FindBin(0.000001);
    int b2 = hres->FindBin(ran-0.00001);
        
    for(int bin =b1;bin<=b2;bin++){
        if(hres->GetBinContent(bin)!=0){
            hResfun->SetBinContent(bin-b1+1,hres->GetBinContent(bin));
            hResfun->SetBinError(bin-b1+1,hres->GetBinError(bin));
        }
    }
    hResfun->SetBinContent(nbins+1,0);  hResfun->SetBinError(nbins+1,1e-7);
    hResfun->SetBinContent(0,0);  hResfun->SetBinError(0,1e-7);
    hResfun->Draw();
    hResfun->Write();
        
/*    for(int bin=1; bin!=nbins;bin++){
    double ResX = h_qn[ih][icent]->GetBinCenter(bin);
  
    int  N = h_qn[ih][icent]->GetBinContent(bin);
    
    for(int j=0;j<N;j++){
        
    double y=myfun(ResX,ResX,sig);
        
        Resfun->Fill(ResX,y);
    }

    }
    Resfun->Draw();
    Resfun->Write();

    double nentries=h_qn[ih][icent]->GetEntries();
    
    double mean1=h_qn[ih][icent]->GetMean();
    double mean2=h_vnc[ih][icent]->GetMean();
    
    TF1* Resfun = new TF1("Resfun",myfun2,0,ran,3);
    Resfun->SetParameters( mean2, sig,nentries);
    Resfun->Draw();
    Resfun->Write();
*/
    
    
//    fun1->Write();
}

void myunfold()
{
    TH1D* V;
    TH1D* Q;
    TFile* fin = TFile::Open("./unfold.root");
    V = (TH1D*)fin->Get("v_for_unfold");
    Q = (TH1D*)fin->Get("q_for_unfold");
    TH2D*tmp = (TH2D*)fin->Get("mig");
    TH1D* tv=(TH1D*)fin->Get("v_for_train");
    TH1D* tq=(TH1D*)fin->Get("q_for_train");
    
    RooUnfoldResponse response (tq,tv,tmp,"res","res");
    
    RooUnfoldBayes unfold (&response, Q, 8);
    
    TH1D* hReco;
    hReco = (TH1D*) unfold.Hreco();
    
    
    
    TCanvas* canRefold = new TCanvas("unfold","hunfold", 800, 600);
    
    Q->Draw();
    V->Draw("same");
    hReco->Draw("same");
    
    
}


