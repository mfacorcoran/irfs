/** @file Dispersion.cxx
@brief implementation of class Dispersion

$Header$

*/

#include "Dispersion.h"
#include "TH1F.h"
#include "TPad.h"
#include "TPaveStats.h"
#include "TList.h"
#include "TF2.h"
#include "TMath.h"
#include "TFitResultPtr.h"

#include <cmath>
#include <iomanip>
#include <stdexcept>
#include <map>

namespace {
    // histogram parameters
    static double xmin=-15, xmax=15; 
    static int nbins=200;

  static const char* names[]={"norm","ls1", "rs1", "bias", "ls2",  "rs2"};
  static double pindex[]={1.6,0.6};
  static double psplit(1.5); //tuned for a better fit
  static double fitrange[]={-15, 15};
  static int min_entries(10);
  static int fit_tries=3; // try a fit this many times before giving up



  double edisp_func(double * x, double * par)
  {
    
    double t=fabs(x[0]-par[3]);
    double s1(par[1]);
    double s2(par[4]);
    
    //left or right sigma's are required?
    if (x[0]>par[3]) {
      //right side
      s1=par[2];
      s2=par[5];
    }
   
    //for core
    double g1=exp(-0.5*pow(t/s1,pindex[0]));
    //for tails
    double g2=exp(-0.5*pow(t/s2,pindex[1]));
    
    //scale norm
    double nscale=exp( 0.5*(pow(psplit/s2,pindex[1])-
			    pow(psplit/s1,pindex[0])));
    if (t>psplit)
      return par[0]*nscale*g2;
    //else
    return par[0]*g1;
  }

  //Asymetric Generalized Gaussian
  //the function is smooth for a>1 and normalized to 1
  //It is an asymetric gaussian if a=2 and a strict gaussian if a=2 and k=1, with standard deviation equal to s/sqrt(2)
  double g(double  x,double  s,double  a,double  k,double  b)
  {
    double prefactor1  = a/(s*TMath::Gamma(1./a));
    double prefactor2  = k/(1+pow(k,2));
    double value=x-b;
    if(value>=0){
      value*=k/s;
    }else{
      value*=-1/(s*k);
    }
    return prefactor1*prefactor2*std::exp(-pow(value,a));
  }
  //edisp function version 2 : uses the asymetric generalized gaussian
  double edisp_func2(double * x, double * par)
  {
    double BIAS    = par[0];
    double BIAS2    = par[1];
    double F       = par[2];
    double K1      = par[3];
    double K2      = par[4];
    double PINDEX1 = par[5];
    double PINDEX2 = par[6];
    double S1      = par[7];
    double S2      = par[8];
    
    double g1 = g(x[0],S1,PINDEX1,K1,BIAS);
    double g2 = g(x[0],S2,PINDEX2,K2,BIAS2);
    double result = 0.;
    if(F==1){result=g1;}
    else if(F==0){result=g2;}
    else {result=F*g1+(1.-F)*g2;}
    float  bin_width = (xmax-xmin)/nbins;
    
    return bin_width*result;
  }

  
}// anon namespace


double Dispersion::scaleFactor(double energy,double  zdir,std::vector<double> edisp_scaling_pars)
{
  // following numbers determined empirically to roughly 
  // make the 68% containment radius be 1.0 independent of energy
  //use a TFunction so that all this stuff could be read from fits file
   
  //x is McLogEnergy, y is fabs(McZDir)
  static const char funcdef[]="[0]*x*x+[1]*y*y + [2]*x + [3]*y + [4]*x*y + [5]";

  static TF2 edisp_scale_func("edisp_scale_func",funcdef);

  double vars[2]; vars[0]=::log10(energy); vars[1]=::fabs(zdir);

  //reltying on contiguous storing of std::vector here....
  return edisp_scale_func(vars,&edisp_scaling_pars[0]); 

}

const char* Dispersion::parname(int i){return names[i];}

int Dispersion::npars(){return sizeof(names)/sizeof(void*);}

std::vector<std::string>
      Dispersion::pnames(names, names+npars());

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Dispersion::Dispersion(std::string histname, 
		       std::string title,
		       embed_python::Module & py)
                                         
: m_hist( new TH1F(histname.c_str(),  title.c_str(),  nbins, xmin, xmax))
, m_count(0)
{
  hist().GetXaxis()->SetTitle("scaled deviation");

  // choose what edisp version to use -----------------------------
  int edisp_version=1;
  try{
    py.getValue("Edisp.Version", edisp_version);
  } catch(std::invalid_argument &){;}
  m_edisp_version=edisp_version;

  try{
    py.getDict("Edisp.fit_pars", m_parmap);
  } catch(std::invalid_argument &){;}

  m_fitfunc=TF1("edisp-fit", *this, fitrange[0], fitrange[1], m_parmap.size());
  biaslimit = 3;
  bias2limit = 3;
  if(m_edisp_version==1) m_fitfunc.SetRange(-7,7);

  std::map<std::string,std::vector<double> >::const_iterator 
    it = m_parmap.begin();
  for (unsigned int i=0;i<m_parmap.size();i++){
    m_fitfunc.SetParName(i, ((*it).first).c_str());
    std::vector<double> par_values = (*it).second;
    m_fitfunc.SetParameter(i, par_values[0]);
    m_fitfunc.SetParLimits(i, par_values[1], par_values[2]);
    if(i==0) biaslimit = par_values[2];
    if(i==1) bias2limit = par_values[2];
    it++;
  }

  m_fitfunc.SetLineWidth(1);
}

std::vector<std::string> Dispersion::getFitParNames() {
  std::vector<std::string> names;
  std::map<std::string,std::vector<double> >::const_iterator 
    it = m_parmap.begin();
  for(;it!=m_parmap.end();it++){
    names.push_back((*it).first);}
  return names;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Dispersion::~Dispersion()
{}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::fill(double scaled_delta, double weight)
{
    hist().Fill( scaled_delta, weight );
    m_count++;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::summary_title(std::ostream & out)
{
  out << "\n\t\t Dispersion fit summary\n"
      <<"Nothing here yet! (Dispersion.cxx)" << std::endl;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::summarize(std::ostream & out)
{

  out<<"Nothing here yet! (Dispersion.cxx)" << std::endl;
}

void Dispersion::
setScaleFactorParameters(const std::vector<double> & scaling_pars) {
   unsigned int scaling_par_size = scaling_pars.size();
   m_scaling_pars.resize(scaling_par_size);
   m_scaling_pars.clear();
   if ( scaling_par_size != 6) {
      throw std::runtime_error("Dispersion::setScaleFactorParameters:\n"
                               "the scaling parameters must have "
                               "exactly 6 values.");
   }
   for (size_t i(0); i < scaling_par_size; i++) {
      m_scaling_pars[i] = scaling_pars[i];
   }
}

void Dispersion::
getScaleFactorParameters(std::vector<double> & scaling_pars) {
   scaling_pars.clear();
   for (size_t i(0); i < 6; i++) {
      scaling_pars.push_back(m_scaling_pars[i]);
   }
}


double Dispersion::function(double* delta, double* par) {
  //static function does not work with P8 for now....
    return edisp_func(delta,par);
}

double Dispersion::operator()(double* delta, double* par) {
  if(m_edisp_version==1) {
     return edisp_func(delta,par);
  } else { 
    return edisp_func2(delta,par);
  }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::fit(std::string opts)
{
  
  std::cout << "\rProcessing " << hist().GetTitle()<<std::endl;
    TH1F & h = hist(); 

    // normalize the distribution
    double scale = h.Integral();
    if (scale > 0 && hist().GetEntries()>50) { 
        h.Sumw2(); // needed to preserve errors
        h.Scale(1./scale);
    }

    double mybias = h.GetMean();
    if(m_edisp_version==2)
      {
        m_fitfunc.SetParameter(0,mybias);
        m_fitfunc.SetParameter(1,mybias);
        m_fitfunc.SetParLimits(0,mybias-biaslimit,mybias+biaslimit);
        m_fitfunc.SetParLimits(1,mybias-bias2limit,mybias+bias2limit);
      }

    if( m_count > min_entries ) {
      int fitcount=0;
      int fitStatus=-1;
      // fit till convergence or for a set number of times
      // Minuit at times gives up for no reason
      while ((fitStatus!=0)&&(fitcount<fit_tries)) {
	fitcount++;
        TFitResultPtr fitRes=h.Fit(&m_fitfunc,opts.c_str()); // fit only specified range
	fitStatus = fitRes;
      }
    }
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::draw(double ymin, double ymax, bool ylog)
{
//
    // set up appearance, and draw to current pad
    if( ylog) gPad->SetLogy();
    TH1F & h = hist(); 
    h.SetMaximum(ymax);
    h.SetMinimum(ymin);
    h.SetStats(true);
    h.SetLineColor(kRed);
    h.GetXaxis()->CenterTitle(true);

    /// @todo: do I need this? why is it here
    TList * list = h.GetListOfFunctions();
    TPaveStats  *s = (TPaveStats*)list->FindObject("stats");
    if( s!= 0 ){
        s->SetY1NDC(0.6);
    }
    h.Draw();

    h.Write();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void Dispersion::getFitPars(std::vector<double>& params)const
{
    params.resize(m_fitfunc.GetNpar());
    // weird, but it seems to work 
    const_cast<TF1*>(&m_fitfunc)->GetParameters(&params[0]);
}
