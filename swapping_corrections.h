class Swap_corr{
 private:
  double           radius;
  dataType         type;
  int              ncent;
  static const int cent_min[5];
  static const int cent_max[5];
  TFile            *f_swap;
  TF1              *fit_p12[5];
  TF1              *fit_p23[5];
  
 public:
  Swap_corr(dataType type, int radius){
   this->radius=radius;
   this->type=type;
   
   if(!(radius==2 || radius==3 || radius==4 || radius==5)){
    cout<<"second parameter radius has to be 2,3,4 or 5, R=radius/10"<<endl;
   }
   
   switch(type){
    case kPbPb_data: case kPbPb_MC:
     f_swap = new TFile(Form("SwappingCorrections/swapPlot_v5_R%d.root",radius));
	   ncent=5;
     break;
    case kPP_data: case kPP_MC:
	   f_swap = new TFile(Form("SwappingCorrections/swapPlot_PP_v5_R%d.root",radius));
     ncent=1;
 	   break;
    default:
	   cout<<"type can take the values kPbPb_data, kPbPb_MC, kPP_data or kPP_MC"<<endl;
	   break;
    }
    for(int icent=0; icent<ncent; icent++){
     fit_p12[icent]=(TF1*)f_swap->Get(Form("fPLeading_%d_%d",cent_min[icent],cent_max[icent]));
     fit_p23[icent]=(TF1*)f_swap->Get(Form("fPSubleading_%d_%d",cent_min[icent],cent_max[icent]));
    }
   }
 
   double get_swap_corr_12(double Aj, int cent=0){
    int cent_bin = 0;
  
    if(type==kPP_data || type==kPP_MC) cent=0;
  
    for(int icent=0; icent<ncent; icent++){
     if(cent>=cent_min[icent] && cent<cent_max[icent]) cent_bin=icent;
    }
  
    double p = fit_p12[cent_bin]->Eval(Aj);
    return p;
   } 
   
   double get_swap_corr_23(double Aj, int cent=0){
    int cent_bin = 0;
  
    if(type==kPP_data || type==kPP_MC) cent=0;
  
    for(int icent=0; icent<ncent; icent++){
     if(cent>=cent_min[icent] && cent<cent_max[icent]) cent_bin=icent;
    }
  
    double p = fit_p23[cent_bin]->Eval(Aj);
    return p;
   } 
};


 const int Swap_corr::cent_min[5]={ 0,20, 60,100,140};
 const int Swap_corr::cent_max[5]={20,60,100,140,200};