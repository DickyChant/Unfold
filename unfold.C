#include <iostream>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TH2.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TLine.h>
#include <TNtuple.h>
#include <TLegend.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TView3D.h>
#include <TTree.h>
#include <TLatex.h>
#include <TROOT.h>
#include <stdio.h>
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TRandom.h"
#include "TComplex.h"

#include "TUnfoldDensity.h"
#include "TUnfoldBinning.h"

#define BIN_NUM 50

void Divide(TCanvas* can,int x,int y,float marx, float mary){
  double xcoor1[10][10], xcoor2[10][10], ycoor1[10][10], ycoor2[10][10];
  Double_t xlow,ylow,xup,yup;
  double ratx[]={1,1,1,1,1,1,1,1};
  double raty[]={1,1,1,1,1,1,1,1};
  double fracx[10], fracy[10];

  double xsli=0,ysli=0;//for boundary
  //define the slice size;
  ratx[0]   +=marx;
  raty[y-1] +=mary;
  for(int i=0;i<x;i++){
    xsli+=ratx[i];
  }
  for(int i=0;i<y;i++){
    ysli+=raty[i];
  }
  fracx[0]=0;  fracy[0]=1;
  for(int i=1;i<=x;i++){
    fracx[i]=fracx[i-1]+ratx[i-1]/xsli;
  }
  for(int i=1;i<=y;i++){
    fracy[i]=fracy[i-1]-raty[i-1]/ysli;
  }
  //rescale
  double scal=0.995;
  for(int i=0;i<=x;i++){
    fracx[i]= fracx[i]*scal+(1-scal)*(0.5-fracx[i]);
  }
  for(int i=0;i<=y;i++){
    fracy[i]= fracy[i]*scal+(1-scal)*(0.5-fracy[i]);
  }
  can->cd();
  can->Divide(x,y);
  int count=1;
  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      can->cd(count);
      count++;
      xlow = fracx[j];      xup = fracx[j+1];
      ylow = fracy[i+1];    yup = fracy[i];
      xcoor1[i][j] = xlow;      xcoor2[i][j] = xup;
      ycoor1[i][j] = ylow;      ycoor2[i][j] = yup;
      //cout<<xlow<<" "<<ylow<<" "<<xup<<" "<<yup<<endl;
      gPad->SetPad(xlow,ylow,xup,yup);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.01);
      gPad->SetTopMargin(0.06);
      gPad->SetBottomMargin(0.01);
      if(j==0){
        gPad->SetLeftMargin((marx+0.15)/ratx[0]);
      }
      if(i==y-1){
        gPad->SetBottomMargin(mary/raty[y-1]);
      }
    }
  }
}

void Divide1(TCanvas* can,int x,int y,float marx, float mary){
  double xcoor1[10][10], xcoor2[10][10], ycoor1[10][10], ycoor2[10][10];
  Double_t xlow,ylow,xup,yup;
  double ratx[]={1,1,1,1,1,1,1,1};
  double raty[]={1,1,1,1,1,1,1,1};
  // double raty[]={1,0.5,1,0.5,1,0.5,1,0.3};
  double fracx[10], fracy[10];

  double xsli=0,ysli=0;//for boundary
  //define the slice size;
  ratx[0]   +=marx;
  raty[y-1] +=mary;
  for(int i=0;i<x;i++){
    xsli+=ratx[i];
  }
  for(int i=0;i<y;i++){
    ysli+=raty[i];
  }
  fracx[0]=0;  fracy[0]=1;
  for(int i=1;i<=x;i++){
    fracx[i]=fracx[i-1]+ratx[i-1]/xsli;
  }
  for(int i=1;i<=y;i++){
    fracy[i]=fracy[i-1]-raty[i-1]/ysli;
  }
  //rescale
  double scal=0.995;
  for(int i=0;i<=x;i++){
    fracx[i]= fracx[i]*scal+(1-scal)*(0.5-fracx[i]);
  }
  for(int i=0;i<=y;i++){
    fracy[i]= fracy[i]*scal+(1-scal)*(0.5-fracy[i]);
  }
  can->cd();
  can->Divide(x,y);
  int count=1;
  for(int i=0;i<y;i++){
    for(int j=0;j<x;j++){
      can->cd(count);
      count++;
      xlow = fracx[j];      xup = fracx[j+1];
      ylow = fracy[i+1];    yup = fracy[i];
      xcoor1[i][j] = xlow;      xcoor2[i][j] = xup;
      ycoor1[i][j] = ylow;      ycoor2[i][j] = yup;
      //cout<<xlow<<" "<<ylow<<" "<<xup<<" "<<yup<<endl;
      gPad->SetPad(xlow,ylow,xup,yup);
      gPad->SetLeftMargin(0.15);
      gPad->SetRightMargin(0.01);
      gPad->SetTopMargin(0.01);
      gPad->SetBottomMargin(0.15);
    }
  }
}
void setstyle(TGraph*h){
  h->GetYaxis()->SetNdivisions(505);  h->GetXaxis()->SetNdivisions(505);
  h->GetYaxis()->SetTitleOffset(1.2); h->GetXaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(17);
  h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(15);
  h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(17);
  h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(15);
  //h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
  //h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(13);
}
void setstyle1(TGraph*h){
  h->GetYaxis()->SetNdivisions(505);  h->GetXaxis()->SetNdivisions(505);
  h->GetYaxis()->SetTitleOffset(1.2); h->GetXaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(15);
  h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(13);
  h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(15);
  h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(13);
  //h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
  //h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(13);
}


void setstyle(TH1*h){
  h->GetYaxis()->SetNdivisions(505);  h->GetXaxis()->SetNdivisions(505);
  h->GetYaxis()->SetTitleOffset(1.5); h->GetXaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(17);
  h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(15);
  h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(17);
  h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(15);
  //h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
  //h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(13);
}
void setstyle1(TH1*h){
  h->GetYaxis()->SetNdivisions(505);  h->GetXaxis()->SetNdivisions(505);
  h->GetYaxis()->SetTitleOffset(1.2); h->GetXaxis()->SetTitleOffset(1.2);
  h->GetYaxis()->SetTitleFont(43);    h->GetYaxis()->SetTitleSize(15);
  h->GetYaxis()->SetLabelFont(43);    h->GetYaxis()->SetLabelSize(13);
  h->GetXaxis()->SetTitleFont(43);    h->GetXaxis()->SetTitleSize(15);
  h->GetXaxis()->SetLabelFont(43);    h->GetXaxis()->SetLabelSize(13);
  //h->GetZaxis()->SetTitleFont(43);    h->GetZaxis()->SetTitleSize(17);
  //h->GetZaxis()->SetLabelFont(43);    h->GetZaxis()->SetLabelSize(13);
}
enum{
  NBIN1=50,
  NBIN2=18,
  NBIN=NBIN1+NBIN2,

  NT=14,
  NM=10,
};

char name[200];
const double pi2 = 2*TMath::Pi();
const double pi = TMath::Pi();
//p1 is smearing width in 1D
//p0 is mean
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
void unfold(){
  float v=0.08,vs=0.035,sig=0.025;
    
    TFile *f = new TFile("unfold.root","recreate");

  TF1*fun = new TF1("fun",bg2,0,0.3,2);
  fun->SetParameters(v,vs);

  TRandom*gran = new TRandom();
  TH1*h[10];
  TH1*h2[10][500],*hh[10];
  TH3D*h3[10];

  double wd=vs; if(sig>wd) wd=sig;
  double ran1 = v+wd*6;
  double ran2 = v+wd*4;
  for(int i=0;i<4;i++){
    sprintf(name,"h_%d",i); 
    if(i==0)   h[i] = new TH1D(name,name,500,-ran1+v,ran1); 
    else if(i==1)   h[i] = new TH1D(name,name,500,-ran2,ran2); 
    else  h[i] = new TH1D(name,name,500,0,ran1); 
    h[i]->SetMarkerStyle(1);
    h[i]->SetMarkerSize(0.5);
    h[i]->Sumw2();
  }
  for(int i=0;i<4;i++){
    sprintf(name,"hh_%d",i); 
    if(i==0) hh[i] = new TH2D(name,name,200,0,ran1,200,0,ran1);
    else if(i==1)  hh[i] = new TH2D(name,name,200,0,ran1,200,-ran1+v,ran1);
    else if (i==2)  hh[i] = new TH2D(name,name,200,0,ran1,200,-ran2,ran2); 
    else hh[i] = new TH2D(name,name,200,0,ran1,200,0,ran1);
  }
  hh[0]->SetYTitle("q1");     hh[0]->SetXTitle("v");
  hh[1]->SetYTitle("q2x");    hh[1]->SetXTitle("v");
  hh[2]->SetYTitle("q2y");    hh[2]->SetXTitle("v");
  hh[3]->SetYTitle("q1");     hh[3]->SetXTitle("q2");

  for(int i=0;i<2;i++){
    sprintf(name,"h3_%d",i); 
    if(i==0) h3[i] = new TH3D(name,name,BIN_NUM,-ran1+v,ran1,BIN_NUM,-ran2,ran2,BIN_NUM,0,ran1);
    else h3[i] = new TH3D(name,name,BIN_NUM,0,ran1,BIN_NUM,0,ran1,BIN_NUM,-pi,pi);
    //h3[i] = new TH3D(name,name,BIN_NUM,0,ran1,BIN_NUM,0,ran1,BIN_NUM,-0.5*ran1,0.5*ran1);
    h3[i]->SetMarkerStyle(6);
    h3[i]->SetMarkerColor(1);
  } 

  h3[0]->SetXTitle("q2x");   h3[0]->SetYTitle("q2y");    h3[0]->SetZTitle("q1");
  h3[1]->SetXTitle("q1");    h3[1]->SetYTitle("q2");     h3[1]->SetZTitle("#Delta#Psi (q1+q2)/2");
  
  h3[2]=new TH3D("h3_2","h3_2",BIN_NUM,0,ran1,BIN_NUM,0,ran1,BIN_NUM,-pi,pi);
  h3[2]->SetXTitle("v1");    h3[2]->SetYTitle("v2");     h3[2]->SetZTitle("#Delta#Psi (v1+v2)/2");
  
  h3[3]=new TH3D("h3_3","h3_3",BIN_NUM,0,ran1,BIN_NUM,0,ran1,BIN_NUM,-pi,pi);
  h3[3]->SetXTitle("v1_data");    h3[3]->SetYTitle("v2_data");     h3[3]->SetZTitle("#Delta#Psi (v1+v2)/2_data");
    
    h3[4]=new TH3D("h3_4","h3_4",BIN_NUM,0,ran1,BIN_NUM,-5*vs,5*vs,BIN_NUM,-vs,ran1*ran1);
    h3[4]->SetXTitle("(q1+q2)/2");    h3[4]->SetYTitle("(q1-q2)/2");     h3[4]->SetZTitle("q1.q2*");
    
    h3[5]=new TH3D("h3_5","h3_5",BIN_NUM,0,ran1,BIN_NUM,-5*vs,5*vs,BIN_NUM,-vs,ran1*ran1);
    h3[5]->SetXTitle("(q1+q2)/2_data");    h3[5]->SetYTitle("(q1-q2)/2_data");     h3[5]->SetZTitle("q1.q2*_data");
    
    h3[6]=new TH3D("h3_6","h3_6",BIN_NUM,-ran1,ran1,BIN_NUM,-ran1,ran1,BIN_NUM,-5*vs,5*vs);
    h3[6]->SetXTitle("(q1+q2)/2_cos");    h3[6]->SetYTitle("(q1+q2)/2_sin");     h3[6]->SetZTitle("(q1-q2)/2");
    
    h3[7]=new TH3D("h3_7","h3_7",BIN_NUM,0,ran1,BIN_NUM,-5*vs,5*vs,BIN_NUM,-pi,pi);
    h3[7]->SetXTitle("(q1+q2)/2_data");    h3[7]->SetYTitle("(q1-q2)/2_data");     h3[7]->SetZTitle("theta");
    
    h3[8]=new TH3D("h3_8","h3_8",BIN_NUM,0,ran1,10,-2*vs,2*vs,BIN_NUM,-pi,pi);
    h3[8]->SetXTitle("(q1+q2)/2_data");    h3[8]->SetYTitle("(q1-q2)/2_data");     h3[8]->SetZTitle("theta");
    
  
  


//  TUnfoldBinning *detectorBinning=new TUnfoldBinning("detector");
//  // highest discriminator bin has fine binning
//  TUnfoldBinning *detectorDistribution=
//     detectorBinning->AddBinning("q1q2thedistribution");
//  detectorDistribution->AddAxis("q1",BIN_NUM,0,ran1,
//                                false, // no underflow bin (not reconstructed)
//                                false // no overflow bin
//                                );
//  detectorDistribution->AddAxis("q2",BIN_NUM,0,ran1,
//                                false, // no underflow bin (not reconstructed)
//                                false // no overflow bin (not reconstructed)
//                                );
//  detectorDistribution->AddAxis("qthe",BIN_NUM,-pi,pi,
//                                false, // no underflow bin (empty)
//                                false // no overflow bin (empty)
//                                );
//  /* TUnfoldBinning *detectorExtra=
//     detectorBinning->AddBinning("detectorextra",7,"one;zwei;three");
//     detectorBinning->PrintStream(cout); */
//
//  //=======================================================================
//  // generator level binning
//  TUnfoldBinning *generatorBinning=new TUnfoldBinning("generator");
//
//  TUnfoldBinning *bgrBinning = generatorBinning->AddBinning("v1v2thedistribution");
//  bgrBinning->AddAxis("v1",BIN_NUM,0,ran1,
//                                false, // no underflow bin (not reconstructed)
//                                false // no overflow bin
//                                );
//  bgrBinning->AddAxis("v2",BIN_NUM,0,ran1,
//                                false, // no underflow bin (not reconstructed)
//                                false // no overflow bin (not reconstructed)
//                                );
//  bgrBinning->AddAxis("vthe",BIN_NUM,-pi,pi,
//                                false, // no underflow bin (empty)
//                                false // no overflow bin (empty)
//                                );
//  TH1 *histDataReco=detectorBinning->CreateHistogram("histDataReco");
//  TH1 *histDataTruth=generatorBinning->CreateHistogram("histDataTruth");
//
    TH1D *h1 = new TH1D("q1_data","q1_data",BIN_NUM,0,ran1);
    TH1D *h6 = new TH1D("v1_data","v1_data",BIN_NUM,0,ran1);

  

  for(int ir=0;ir<100000;ir++){
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
      double the = V3.Theta();
      TComplex V4(v1,0);
      TComplex V5=(V4*V4.Conjugate(V4))/v1;
      double the2=V5.Theta();
      if(the>pi) the-=pi2; if(the<-pi) the+=pi2;
      if(the2>pi) the2-=pi2; if(the2<-pi) the2+=pi2;
    
    h3[3]->Fill(q1,q2,the);
    h3[5]->Fill((q1+q4)/2,(q1-q4)/2,q1q2.Re());
    //h3[1]->Fill(q1,q2,the*(q1+q2)/2);
//    Int_t number_dec = detectorDistribution->GetGlobalBinNumber(q1,q2,the);
//    Int_t number_gen = bgrBinning->GetGlobalBinNumber(v1,v1,the2);
    h1->Fill(q1);
    h6->Fill(v1);
    
      
      

    


  }
    h1->Write();
    h6->Write();
    h3[3]->Write();
    h3[5]->Write();
  
//  TH2 *histMCGenRec=TUnfoldBinning::CreateHistogramOfMigrations
//     (generatorBinning,detectorBinning,"histMCGenRec");
    TH2D *histMCGenRec = new TH2D("mig_matrix","mig_matrix",BIN_NUM,0,ran1,BIN_NUM,0,ran1);
    TH2D *histM = new TH2D("mig","mig",10000,0,10000,10000,0,10000);
    TH1D *h5 = new TH1D("q1_simulate","q1_simulate",BIN_NUM,0,ran1);
    TH1D *h4 = new TH1D("v1_simulate","v1_simulate",BIN_NUM,0,ran1);
    
    TH1D *P=new TH1D("v/q_p","v/q_p",50,0,ran1);
    TH1D *M=new TH1D("v/q_m","v/q_m",10,-2*vs,2*vs);
    TH1D *T=new TH1D("v/q_the","v/q_the",20,-pi,pi);
    
    
 
  
  
  for(int ir=0;ir<10000000;ir++){
    double v1 = fun->GetRandom();
      double v2 = fun->GetRandom();
    double the1=gran->Rndm()*pi2; //angle for v1

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

    
    h3[0]->Fill(q2x,q2y,q1);

    h3[1]->Fill(q1,q2,the);
    h3[2]->Fill(v1,v1,the2);
      h3[4]->Fill((q1+q4)/2,(q1-q4)/2,q1q2.Re());
      h3[6]->Fill((q1+q4)/2*cos(The),(q1+q4)/2*sin(The),(q1-q4)/2);
      h3[7]->Fill((q1+q4)/2,(q1-q4)/2,the);
      h3[8]->Fill((q1+q4)/2,(q1-q4)/2,the);
    //h3[1]->Fill(q1,q2,the*(q1+q2)/2);
//    Int_t number_dec = detectorDistribution->GetGlobalBinNumber(q1,q2,the);
//    Int_t number_gen = bgrBinning->GetGlobalBinNumber(v1,v1,the2);
      h5->Fill(v1);
      h4->Fill(q1);
    histMCGenRec->Fill(q1,v1);
      double vp=(v1+v2)/2;
      double vm=(v1-v2)/2;
      double qp=(q1+q4)/2;
      double qm=(q1-q4)/2;
      int v_index=P->FindBin(vp)*200+M->FindBin(vm)*20+T->FindBin(the1);
      int q_index=P->FindBin(qp)*200+M->FindBin(qm)*20+T->FindBin(the);
      
      histM->Fill(v_index,q_index);
      

    


  }
    h3[0]->Write();
    
    h3[1]->Write();
    h3[2]->Write();
    h3[6]->Write();
    h3[7]->Write();
    h3[8]->Write();
    h5->Write();
    h4->Write();
    h3[4]->Write();
    histMCGenRec->Write();
    histM->Write();
    new TCanvas();
    h3[7]->Draw("iso");
    new TCanvas();
    h3[8]->Draw("iso");
    TH1*ht[10],*htmp;
    htmp =  h3[7]->Project3D("zx"); htmp->Write();ht[0] = htmp;
    htmp =  h3[7]->Project3D("zy"); htmp->Write();ht[1] = htmp;
    htmp =  h3[7]->Project3D("yx"); htmp->Write();ht[2] = htmp;
    htmp =  h3[6]->Project3D("zx"); htmp->Write();ht[3] = htmp;
    htmp =  h3[6]->Project3D("zy"); htmp->Write();ht[4] = htmp;
    htmp =  h3[6]->Project3D("yx"); htmp->Write();ht[5] = htmp;
    
    TCanvas *c1 = new TCanvas("c2","c2",1200,900);
    c1->Divide(3,2);
    
    
    
    c1->cd(1);  ht[0]->Draw("colz");
    c1->cd(2);  ht[1]->Draw("colz");
    c1->cd(3);  ht[2]->Draw("colz");
    c1->cd(4);  ht[3]->Draw("colz");
    c1->cd(5);  ht[4]->Draw("colz");
    c1->cd(6);  ht[5]->Draw("colz");
    c1->Write();
    
    
    }
