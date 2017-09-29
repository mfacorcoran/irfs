/** 
 * @file IrfAnalysis.cxx
 * @brief Analysis of instrument response function distributions
 * $Header$
*/

#include "IrfAnalysis.h"
#include "PsfPlots.h"
#include "FisheyePlots.h"
#include "DispPlots.h"
#include "EffectiveArea.h"
#include "AeffPhiDep.h"
#include "TreeWrapper.h"
#include "Dispersion.h"
#include "PointSpreadFunction.h"
#include "embed_python/Module.h"

#include "CLHEP/Geometry/Vector3D.h"

#include "TFile.h"
#include "TTree.h"
#include "TEventList.h"

#include <iomanip>
#include <fstream>
#include <ios>
#include <stdexcept>

namespace {
   inline static double sqr(double x) {
      return x*x;
   }
}

IrfAnalysis::IrfAnalysis(std::string output_folder,
                         embed_python::Module & py) 
   : MyAnalysis(py),
     m_binner(py),  // initilize the binner
     m_make_plots(true),
     m_output_type("pdf"),
     m_output_dir(output_folder),
     m_bestXDir("CTBBestXDir"),
     m_bestYDir("CTBBestYDir"),
     m_bestZDir("CTBBestZDir"),
     m_bestEnergy("CTBBestEnergy"),
     m_front_only_psf_scaling(false) {
   std::string logfile;
   std::string selectionName;

   py.getValue("className", m_classname);    
   py.getValue("logFile", logfile);
   py.getValue("selectionName", selectionName);
   
   try {
     int make_plots = 1;
     py.getValue("makePlots", make_plots);
     m_make_plots = bool(make_plots);
     py.getValue("outputType", m_output_type);
   } catch (std::invalid_argument &) { }

   m_filename_root = selectionName;
   m_outputfile = selectionName + ".root";

   if (logfile.empty()) {
      m_log = &std::cout;
   } else {
      m_log = new std::ofstream(logfile.c_str(), std::ios_base::app);
   }

   py.getValue("parameterFile", m_parameterFile);

   py.getValue("Data.generate_area", m_generate_area);

   try {
      py.getValue("Data.var_xdir", m_bestXDir);
      py.getValue("Data.var_ydir", m_bestYDir);
      py.getValue("Data.var_zdir", m_bestZDir);
      py.getValue("Data.var_energy", m_bestEnergy);
   } catch (std::invalid_argument &) {
      /// These are missing, so use defaults.
   }
   std::cout << "Using variables "
             << m_bestXDir << ", "
             << m_bestYDir << ", "
             << m_bestZDir << ", "
             << m_bestEnergy << std::endl;


   std::vector<double> generated, logemins, logemaxes;
   py.getList("Data.generated", generated);
   py.getList("Data.logemin", logemins);
   py.getList("Data.logemax", logemaxes);
   for (size_t i=0; i<generated.size(); ++i) {
      normalization().push_back(Normalization(generated[i],
                                              logemins[i], logemaxes[i]));
   }

   setName(m_classname);
   current_time(out());
   out() << "Event class is " << m_classname << std::endl;

   project(py);
}

void IrfAnalysis::project(embed_python::Module & py) {

  // If skim file is undefined load events directly from the TChain
  if(skim_filename().empty())
    makeCutTree();
  else
    open_input_file();

   // for the histograms
   TFile * m_hist_file= new TFile(summary_filename().c_str(), "recreate");
   std::cout << " writing irf summary plots to " 
             << summary_filename() << std::endl;
   out() << " writing irf summary plots to " 
         << summary_filename() << std::endl;
   std::cout << std::endl;

   //---- declare the IRF plots------
   //---------------------------
   m_psf = new PsfPlots(*this, out(), py);

   m_disp = new DispPlots(*this, out(), py);

   m_aeff = new EffectiveArea(*this, out());
   m_phi_dep = new AeffPhiDep(*this);

   m_fisheye = new FisheyePlots(*this,out(),py);

   std::cout << "Selecting columns in tree " << tree().GetName() << std::endl;
   unsigned EvtRun(0);
   float McEnergy(0.0), Tkr1FirstLayer(0.0);
   double BestEnergy(0.0), McXDir(0.0), McYDir(0.0), McZDir(0.0),
     BestXDir(0.0), BestYDir(0.0), BestZDir(0.0);

   tree().SetBranchAddress("EvtRun", &EvtRun);
   tree().SetBranchAddress("McEnergy", &McEnergy);
   tree().SetBranchAddress("Tkr1FirstLayer", &Tkr1FirstLayer);
   tree().SetBranchAddress(m_bestEnergy.c_str(), &BestEnergy);
   tree().SetBranchAddress("McXDir", &McXDir);
   tree().SetBranchAddress("McYDir", &McYDir);
   tree().SetBranchAddress("McZDir", &McZDir);
   tree().SetBranchAddress(m_bestXDir.c_str(), &BestXDir);
   tree().SetBranchAddress(m_bestYDir.c_str(), &BestYDir);
   tree().SetBranchAddress(m_bestZDir.c_str(), &BestZDir);
   
   int lastrun(0), selected(0), total(0), nruns(0);
   double minlogE(1e6), maxlogE(0);
   double minzdir(1), maxzdir(-1);

   int nentries(0);
   if(list() != NULL) 
     nentries = list()->GetN();
   else
     nentries = tree().GetEntries();

   for(int i = 0; i < nentries; i++) {

      int ii = i;
      if(list() != NULL) {
	ii = list()->GetEntry(i);
      }
      tree().GetEvent(ii);
      if (EvtRun != lastrun) {
         ++nruns;
         lastrun = EvtRun;
      }
      double logE(log10(McEnergy));
      if (logE < minlogE) {
         minlogE = logE;
      }
      if (logE > maxlogE) {
         maxlogE = logE;
      }
      if (McZDir < minzdir) {
         minzdir = McZDir;
      }
      if (McZDir > maxzdir) {
         maxzdir = McZDir;
      }
      
      bool front = Tkr1FirstLayer > 5;
      ++selected;

      // calculate theta, phi components of the error
      HepGeom::Vector3D<double>
         mc_dir(McXDir, McYDir, McZDir), 
   	fit_dir(BestXDir, BestYDir, BestZDir);

      //Approximated version to check theta- and phi- projections
      HepGeom::Vector3D<double> 
   	mc_error(mc_dir - fit_dir), zhat(0, 0, 1),
   	phi_hat = zhat.cross(mc_dir).unit(),
   	theta_hat = phi_hat.cross(mc_dir).unit();

      double phi_err = mc_error*phi_hat,
   	theta_err = mc_error*theta_hat;
   	//	diff = sqrt(sqr(theta_err) + sqr(phi_err));

      //exact version
      double diff=mc_dir.angle(fit_dir);
      
      // choose one of the following
      double measured_energy = BestEnergy, 
         mc_energy = McEnergy,
         //ratio =measured_energy/mc_energy,
         dsp = measured_energy/mc_energy - 1;

      m_fisheye->fill(theta_err, McEnergy, McZDir);
      m_psf->fill(diff, McEnergy, McZDir);
      m_disp->fill(dsp, McEnergy, McZDir);
      m_aeff->fill(mc_energy, McZDir, front, total);
      m_phi_dep->fill(McXDir, McYDir, McEnergy, McZDir);
   }
     
   out() << "\nFound " << nruns <<" run numbers" 
         << " and " << selected<< "/" 
         <<  tree().GetEntries() << " events" <<  std::endl;
   out() << "Log energy range: " << minlogE << " to " << maxlogE << std::endl;
   out() << "McZDir range: " << minzdir << " to " << maxzdir << std::endl;
   
   m_hist_file->Write();}

void IrfAnalysis::fit(bool make_plots) {

   make_plots = make_plots || m_make_plots;
   const std::string& output_type = m_output_type;
   m_psf->fit(); 
   m_fisheye->fit(); 
   m_disp->fit();
   m_phi_dep->fit();
   
   m_psf->summarize();
   m_disp->summarize();
   m_aeff->summarize();
   if (make_plots) {
      m_psf->draw(m_filename_root + "_psf." + output_type);
      m_fisheye->draw(m_filename_root + "_fisheye." + output_type);
      m_disp->draw(m_filename_root + "_disp." + output_type);
      m_aeff->draw(m_filename_root + "_aeff." + output_type);
      m_phi_dep->draw(m_filename_root + "_phi_dep." + output_type);
   }
   if (!m_parameterFile.empty()) {
      writeFitParameters(m_output_dir + "/" + m_parameterFile);
   }
}

void IrfAnalysis::writeFitParameters(std::string outputFile) { 
   // Open the output file, overwriting any existing file.
    TFile* file = new TFile(outputFile.c_str(), "RECREATE");

    // add directory if not already there
    TDirectory* td = file->mkdir((m_classname).c_str());

    std::string classpath=std::string("/")+m_classname;
    bool check = file->cd(classpath.c_str());
    if( ! check) throw("could not cd to "+classpath);

    out() << "Writing summaries to " << outputFile << "/"
          << m_classname << std::endl;

    // now make 2-D histograms of the values
    TDirectory* psfdir=td->mkdir("PSF");
    check = psfdir->cd();
    if( ! check) throw(strcat("could not cd to ",psfdir->GetPath()));
    m_psf->fillParameterTables();

    TDirectory* fisheyedir=td->mkdir("FISHEYE");
    check = fisheyedir->cd();
    if( ! check) throw(strcat("could not cd to ",fisheyedir->GetPath()));
    m_fisheye->fillParameterTables();

    TDirectory* dispdir=td->mkdir("EDISP");
    check = dispdir->cd();
    if( ! check) throw(strcat("could not cd to ",dispdir->GetPath()));
    m_disp->fillParameterTables();

    TDirectory* aeffdir=td->mkdir("AEFF");
    check = aeffdir->cd();
    if( ! check) throw(strcat("could not cd to ",aeffdir->GetPath()));
    m_aeff->fillParameterTables();

    TDirectory* phidir=td->mkdir("PHIDEP");
    check = phidir->cd();
    if( ! check) throw(strcat("could not cd to ",phidir->GetPath()));
    m_phi_dep->fillParameterTables();

    delete file;
    current_time(out());
}

