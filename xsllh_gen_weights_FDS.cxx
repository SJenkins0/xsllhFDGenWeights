#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <unistd.h>
#include <list>

#include "TAxis.h"
#include "TChain.h"
#include "TClonesArray.h"
#include "TFile.h"
#include "TH1.h"
#include "TObjString.h"
#include "TTree.h"
#include "TVector3.h"

#include "toml/toml_helper.h"
#include "neutvect.h"

#include "ND__GRooTrackerVtx.h"
#include "ND__NRooTrackerVtx.h"

#include "T2KReWeight/Interface/T2KSyst.h"
#include "T2KReWeight/WeightEngines/NIWG/T2KNIWGUtils.h"
#include "T2KReWeight/WeightEngines/NEUT/T2KNEUTUtils.h"
#include "T2KReWeight/WeightEngines/T2KReWeightFactory.h"


const std::string RESET("\033[0m");
const std::string RED("\033[31;1m");
const std::string GREEN("\033[92m");
const std::string COMMENT_CHAR("#");

const std::string TAG = GREEN + "[xsReWeight]: " + RESET;
const std::string ERR = RED + "[ERROR]: " + RESET;


/****************************
 * This struct is the basic event container holding the variables needed for
 * both the weight generation and the later steps for generating splines.
 * Add the variables needed for your analysis that are contained in your
 * HighLand files.
 */
struct XSEvent
{
  //Basic set of variables
  //Recon tree
  int topology, reaction, sample, target, fgdtarget;
  int nutype, layer;
  float enu_true, enu_reco;
  float dirnu_true[3];
  float pmu_true, pmu_reco, pmu_seltrue, ppi_true, ppi_reco, ppi_seltrue;
  float cosmu_true, cosmu_reco, cosmu_seltrue, cospi_true, cospi_reco, cospi_seltrue;
  float dirmu_seltrue[3];
  float q2_true;
  float w_true;
  int selection;
  
  //Do NOT modify these
  float weight_nom;
  float* weight_syst;

  float weight_truth;
  
  //Do NOT modify these
  int run, subrun;
  int input_RooVtxIndex;
  int input_RooVtxEntry;
  int TruthVertexID;
  ND::NRooTrackerVtx* vtx;
  
  XSEvent()
  {
    weight_syst = nullptr;
    vtx = nullptr;
  }
};

/****************************
 * Updated this from Andrew's original implementation, to serve as a struct
 * for each indivdual fake data study. impl_name, nominal value and dial
 * values are now vectors, to take values for each dial used in the FDS.
 *
 * This struct holds information about the dial(s) being used for weight
 * generation. The 'user_name' is a informal name defined by the analyzer
 * in the TOML config, and is meant for better ease of use. The 'impl_name'
 * is the exact name as defined in T2KRW. The values array holds the dial
 * values when generating weights.
 */
struct FDS
{
  std::string user_name;
  std::vector<std::string> impl_name;
  int n_dials;
  int steps;
  std::vector<float> nominal;
  std::vector<std::vector<float> > values;
};

/****************************
 * The following functions are where most analysis modifications need to take
 * place. The output branches are for the output ROOT files containing the
 * weights and variables needed for splines. The input branches are for your
 * HighLand files. Both functions use the XSEvent as a method to pass information
 */
void setOutputBranches(TTree& t, XSEvent& evt, bool use_truth_tree);
void setInputBranches(TChain& c, XSEvent& evt, bool use_truth_tree, int nsel); //Add selection number for different category names between FGD1 and 2


/****************************
 * This function takes an XSEvent, the accum_level array, and a flag for if you are
 * using the truth tree (where selection cuts do not make sense). It returns a bool
 * corresponding to if the event passed (true) or failed (false) the cuts.
 */
template <typename A>
bool doesEventPassCut(XSEvent& xs_evt, A accum_level, int isel, bool truth_tree)
{
    if(truth_tree)
        return true;

    bool event_pass = false;
    const std::vector<int> cut_level = {10, 10, 10, 9, 8, 8}; //hard coded - set this!

    for(unsigned int s = 0; s < cut_level.size(); ++s)
      {
	if(accum_level[0][isel][s+1] > cut_level[s]) //s+1 to offset for 0th branch which is 1pi total - this will need altering for most analysis to remove the +1
	  {
	    event_pass = true;
	    xs_evt.sample = s; //Do the layer splitting in the tree converter rather than here
	    xs_evt.selection = isel;
	    break;
	  }
      }
    
    return event_pass;
}

/****************************
 * A couple utility functions to facilitate building the list of input ROOT files.
 */
std::string findFileExtension(const std::string& file);
std::vector<std::string> getInputRootFiles(const std::vector<std::string>& file_list);

/****************************
 * Welcome to the main function.
 */
int main(int argc, char** argv)
{
    std::cout << "--------------------------------------------------------\n"
              << TAG << "Welcome to the Super-xsLLhReWeight Interface.\n"
              << TAG << "Initializing the reweight machinery..." << std::endl;

    bool update_toml = false;
    bool use_truth_tree = false;
    int max_events = -9999;
    std::string config_toml;

    char option;
    while((option = getopt(argc, argv, "n:t:Th")) != -1)
      {
        switch(option)
	  {
	  case 'n':
	    max_events = std::stol(optarg);
	    break;
	  case 't':
	    config_toml = optarg;
	    break;
	  case 'T':
	    update_toml = true;
	    break;
	  case 'h':
	    std::cout << "USAGE: " << argv[0] << "\nOPTIONS\n"
		      << "-n : Max number of events to process\n"
		      << "-t : TOML config file\n"
		      << "-T : Update TOML config file with nominal values\n"
		      << "-h : Display this help message\n";
	  default:
	    return 0;
	  }
      }
    
    if(config_toml.empty())
      {
        std::cout << ERR << "Missing necessary command line arguments.\n" << std::endl;
        return 1;
      }
    
    std::cout << TAG << "Reading toml file..." << std::endl;
    auto opts = toml::parse(config_toml);

    auto neut_card = toml::find<std::string>(opts, "neutcard");
    std::cerr << neut_card << std::endl;
    if(neut_card.size())
        t2krew::T2KNEUTUtils::SetCardFile(neut_card);

    std::cout << TAG << "Getting T2KRW instance from the factory." << std::endl;
    auto T2KRW = t2krew::MakeT2KReWeightInstance(t2krew::Event::kNRooTracker);
    auto generate_list = toml::find<std::vector<std::string>>(opts, "generate");
    auto& fakedata_list = toml::find(opts, "FDS");

    std::vector<FDS> v_FDS;
    for(const auto& fakedata : generate_list)
    {
        std::cout << TAG << "Trying to find '" << fakedata << "' ..." << std::endl;
        toml::value d_toml;
        try
        {
            d_toml = toml::find(fakedata_list, fakedata);
        }
        catch(...)
        {
            std::cout << TAG << "Fake data study '" << fakedata << "' not found in toml file! Skipping."
                      << std::endl;
            continue;
        }

        std::cout << "Fake data study found. Reading parameters." << std::endl;
	
        //auto d_steps  = toml::find<double>(d_toml, "steps");
	auto d_steps = -1;
        auto d_values = toml::find_or(d_toml, "values", std::vector<std::vector<float> >{});
        auto d_impl   = toml::find<std::vector<std::string> >(d_toml, "impl_name");
	int d_dials   = d_impl.size(); //Get the number of dials we're using in the FDS


	std::cout << d_dials << std::endl;
	std::cout << d_values.size() << std::endl;
	std::cout << d_impl.size() << std::endl;

        T2KRW->Reset();
        T2KRW->Reconfigure();

	std::vector<float> d_nominal;

	for(int d=0; d<d_dials; d++){
	  t2krew::T2KSyst_t temp_rw_dial = T2KRW->DialFromString(d_impl[d]);
	  d_nominal.push_back(T2KRW->GetDial_From_Value(temp_rw_dial));
	}

        auto& table = toml::get<toml::table>(fakedata_list).at(fakedata);
        table["nominal"] = d_nominal;

        if(d_values.empty())
        {
	    std::cout << "Fake data cannot be generated. Exiting..." << std::endl;
	    exit(1);
        }
        else
        {
            std::cout << "Found dial values array." << std::endl;
	    d_steps = d_values.at(0).size(); //Currently making sure all dials have the same number of steps
        }

        FDS temp;
        temp.user_name = fakedata;
        temp.impl_name = d_impl;
        temp.nominal = d_nominal;
	temp.steps = d_steps;
        temp.values = d_values;
	temp.n_dials = d_dials;
        v_FDS.push_back(temp);

	std::cout << "Settings for " << temp.user_name << ":" << std::endl;
	std::cout << "Number of dials used: " << temp.n_dials << std::endl;
	for(int d=0; d<d_dials; d++){
	  std::cout << "Dial name: " << d_impl.at(d) << std::endl;
	  std::cout << "Steps    : " << d_steps << std::endl;
	  std::cout << "Nominal  : " << d_nominal.at(d) << std::endl;
	  std::cout << "Values   : ";
	  for(const auto& v : d_values.at(d))
            std::cout << v << " ";
	  std::cout << std::endl;
	}
        std::cout << std::endl;
    }

    if(update_toml)
    {
        std::cout << TAG << "Writing updated toml configuration file." << std::endl;
        std::ofstream toml_out;
        toml_out.open("output.toml");
        toml_out << std::setw(80) << opts;
        toml_out.close();
    }

    for(const auto& f : v_FDS)
      {      
	
	//Set up the output tree for the specific fake data study
	auto output_prefix = toml::find_or(opts, "output_prefix", std::string(""));
	const std::string output_name = output_prefix + f.user_name + ".root";
	TFile* file_output = TFile::Open(output_name.c_str(), "RECREATE");

	auto file_list = toml::find<std::vector<std::string>>(opts, "filelist");

	//Now we do a loop of two iterations
	//First is over the selected events, the second is for the truth
	for (unsigned int it=0; it<2; it++){
      
	  use_truth_tree = (bool)it; //This isn't strictly necessary..
	  std::string tree_name;
	  if(use_truth_tree){
	    std::cout << TAG << "Running over truth tree events." << std::endl;
	    tree_name = "truth";
	  }
	  else{
	    std::cout << TAG << "Running over default events, accumulation level cuts will be applied." << std::endl;
	    tree_name = "default";
	  }

	  if(max_events < 0)
	    max_events = toml::find_or(opts, "max_events", -9999);

	  TChain tree_input(tree_name.c_str());

	  auto root_list = getInputRootFiles(file_list);
	  if(root_list.empty())
	    {
	      std::cerr << ERR << "No input ROOT files. Exiting." << std::endl;
	      return 121;
	    }

	  for(const auto& file : root_list)
	    {
	      std::cout << TAG << "Adding " << file << " to TChain." << std::endl;
	      tree_input.Add(file.c_str());
	    }

	  const int total_events = tree_input.GetEntries();
	  const int nevents = max_events > 0 && max_events < total_events ? max_events : total_events;
	  std::cout << TAG << "Total events in chain: " << total_events << std::endl;
	  std::cout << TAG << "Processing " << nevents << " events." << std::endl;
	  if(nevents < 0)
	    {
	      std::cerr << ERR << "Number of events is negative \u203d" << std::endl;
	      std::cerr << ERR << "Exiting..." << std::endl;
	      return 121;
	    }
	  std::vector<XSEvent> vec_evts;
	  vec_evts.reserve(nevents);
    
	  const unsigned int nsamples = 6; //Hardcoded
	  const unsigned int nselections = 2; //Hardcoded
	  int accum_level[1][nselections][nsamples+1]; //Offset for total branch in the highland input tree - alter for most analyses

    
	  int NRooVtx{0}, vtx_run{0}, vtx_subrun{0};
	  std::list<ND::NRooTrackerVtx> v_vtxs;

	  TClonesArray* nRooVtxs = new TClonesArray("ND::NRooTrackerVtx");

	  TFile* current_file = nullptr;
	  TTree* tree_RooVtx = nullptr;
	  UInt_t current_hash = 0;


   
	  //Loop over number of selections, as we need to set different
	  //branch addresses dependent on FGD.
	  int nFail = 0;
	  for (unsigned int isel = 0; isel < nselections; isel++){
	    XSEvent xs_evt;
	    std::cout << TAG << "Setting branch addresses for selection " << isel << std::endl;
      
	    //Beware - this now resets all tree_input branch addresses before setting new ones 
	    setInputBranches(tree_input, xs_evt, use_truth_tree, isel); //Set the branch address for the ith sel
	    tree_input.SetBranchAddress("accum_level", accum_level);
      
	    std::cout << TAG << "Reading events from chain..." << std::endl;
	    for(int i = 0; i < nevents; ++i)
	      {
		tree_input.GetEntry(i);

		/************************    
		 * Check if the TChain has changed files by comparing the file hash        
		 * to the previous known file hash. If different, the file has changed 
		 * and set up the new RooVtx TTree.
		 */
		if(tree_input.GetFile()->Hash() != current_hash)
		  {
		    std::cout << TAG << "Changing file to " << tree_input.GetFile()->GetName() << std::endl;
		    current_file = tree_input.GetFile();
		    current_hash = tree_input.GetFile()->Hash();
		    tree_RooVtx = (TTree*)current_file->Get("NRooTrackerVtx");
		    tree_RooVtx->SetBranchAddress("Vtx", &nRooVtxs);
		    tree_RooVtx->SetBranchAddress("NVtx", &NRooVtx);
		    tree_RooVtx->SetBranchAddress("RunID", &vtx_run);
		    tree_RooVtx->SetBranchAddress("SubrunID", &vtx_subrun);
		  }

		tree_RooVtx->LoadTree(xs_evt.input_RooVtxEntry);
		tree_RooVtx->GetEntry(xs_evt.input_RooVtxEntry);
	  
		if(i % 1000 == 0){
		  std::cout << TAG << "Processing event " << i << std::endl;
		  std::cout << TAG << "Current file " << tree_input.GetFile()->GetName() << std::endl;
		  std::cout << TAG << "Current hash: " << tree_input.GetFile()->Hash() << std::endl;
		}
		
		//Lets try doing the event selection first, to see how this affects
		//number of events passing or failing cuts
		if(!doesEventPassCut(xs_evt, accum_level, isel, use_truth_tree)) continue;

#ifdef VTX_FAST
		ND::NRooTrackerVtx* vtx = (ND::NRooTrackerVtx*)nRooVtxs->At(xs_evt.input_RooVtxIndex);
#else
		ND::NRooTrackerVtx* vtx = nullptr;
		for(int v = 0; v < NRooVtx; ++v)
		  {
		    vtx = (ND::NRooTrackerVtx*)nRooVtxs->At(v);

		    if(vtx->TruthVertexID == xs_evt.TruthVertexID)
		      {

#ifdef DEBUG_MSG
			std::cout << "NVtx = " << NRooVtx << std::endl;
			std::cout << "Vtx -> " << vtx->TruthVertexID << " vs. " << xs_evt.TruthVertexID << std::endl;
			std::cout << "RooVtxIdx: " << xs_evt.input_RooVtxIndex << " vs. " << v << std::endl;
			std::cout << "RunID   : " << vtx_run << " vs. " << xs_evt.run << std::endl;
			std::cout << "SubRunID: " << vtx_subrun << " vs. " << xs_evt.subrun << std::endl;
			std::cout << "Nu E    : " << vtx->StdHepP4[0][3] << " vs. " << xs_evt.enu_true << std::endl;
			std::cout << "EvtCode : " << vtx->EvtCode->GetString().Atoi() << std::endl;
			std::cout << "HepPdg  : " << vtx->StdHepPdg[0] << ", " << vtx->StdHepPdg[1] << std::endl;
#endif
			break;
		      }
		    vtx = nullptr;
		  }//end loop over NRooTrackerVtxs
#endif

		if(vtx == nullptr){
		  nFail++;
		  continue;
		}


		//Calculate the true outgoing angle for the muon candidate
		TVector3 mu(xs_evt.dirmu_seltrue[0], xs_evt.dirmu_seltrue[1], xs_evt.dirmu_seltrue[2]);
		TVector3 nu(xs_evt.dirnu_true[0], xs_evt.dirnu_true[1], xs_evt.dirnu_true[2]);
		double mu_nu = mu.Angle(nu);
		xs_evt.cosmu_seltrue = TMath::Cos(mu_nu);

		    //Take all of these out to begin with - reweighted tree should be kinematically
		    //identical to the input. Will add them back in if there are issues

		    /*
		      if(xs_evt.pmu_reco < 0 || xs_evt.pmu_reco > 30000 || xs_evt.ppi_reco < 0)
		      continue;
		    */
		    /*
		    //If the true pion kinematics are invalid, fall back to the selected pion truth info
		    //Shouldn't need to do this for muon but can be implemented if necessary
		    if (xs_evt.ppi_true < 0 || xs_evt.cospi_true < -2){
		    xs_evt.ppi_true = xs_evt.ppi_seltrue;
		    xs_evt.cospi_true = xs_evt.cospi_seltrue;
		    }

		    //Small number of events that this doesnt fix - ignore these
		    if (xs_evt.ppi_true < 0) continue;
		    */
	      
		//Add the vertices
		v_vtxs.push_back(*vtx);
		vec_evts.push_back(xs_evt);
		vec_evts.back().vtx = &v_vtxs.back();
	      
		  
	      }
	  }

	  std::string tree_output_name;
	  if(use_truth_tree) tree_output_name = "trueEvents";
	  else tree_output_name = "selectedEvents";
	  std::cout << TAG << "Events passing cuts: " << vec_evts.size() << std::endl;
	  std::cout << TAG << "Events where the TruthVertexIDs couldn't be matched: " << nFail << std::endl;

	  TTree tree_output(tree_output_name.c_str(), tree_output_name.c_str());
	  tree_output.SetDirectory(file_output);
	
	  XSEvent tmp_evt;
	  setOutputBranches(tree_output, tmp_evt, use_truth_tree);

	  const int dialnum = f.n_dials;
	  std::cout << "Number of dial values = " << f.values.at(0).size() << std::endl;
	  std::cout << "n_dials here = " << f.n_dials << std::endl;
	  std::cout << "dial_num here = " << dialnum << std::endl;
	  float w[f.values.at(0).size()];
	  float x[dialnum][f.values.at(0).size()];
	  const std::string weight_str = "weight_syst[" + std::to_string(f.values.at(0).size()) + "]/F";
	  const std::string dialval_str = "dial_values[" + std::to_string(f.n_dials) + "][" + std::to_string(f.values.at(0).size()) +"]/F";
	  tree_output.Branch("weight_syst", w, weight_str.c_str());
	  tree_output.Branch("dial_values", x, dialval_str.c_str());

	  std::cout << TAG << "Generating weights for " << f.user_name << std::endl;
	  
	  std::vector<t2krew::T2KSyst_t> t2ksyst_dials;
	  for(int d_idx = 0; d_idx<f.n_dials; d_idx++)
	    t2ksyst_dials.push_back(T2KRW->DialFromString(f.impl_name[d_idx]));
	  
	  for(auto& e : vec_evts)
	    {
	  
	      auto rw_evt = t2krew::Event::Make(e.vtx);
	      //Loop over index in the dial value array - this is the same for each dial
	      //Each position will correspond to a single fake data set, whether it requires 1 or more dials set at once
	      for(unsigned int v = 0; v < f.values.at(0).size(); ++v)
		{
		  for(int d_idx=0; d_idx<f.n_dials; d_idx++){
		    
		    T2KRW->SetDial_To_Value(t2ksyst_dials.at(d_idx), f.values.at(d_idx).at(v));
		    x[d_idx][v] = f.values.at(d_idx).at(v);
		  }
		  T2KRW->Reconfigure();
		  
		  w[v] = T2KRW->CalcWeight(rw_evt);
		
		}
	      
	      tmp_evt = e;
	      file_output->cd();
	      tree_output.Fill();
	      
	    }
	  
	  T2KRW->Reset();
	  
	  tree_output.Write();
	} 
	file_output->Close();
	
	std::cout << TAG << "Finished weights for dial." << std::endl;
	std::cout << TAG << "File written: " << output_name << std::endl;
      }    
    
    /************************
     * Ideally the execution reaches this point.
     */
    std::cout << TAG << "Finished." << std::endl;
    return 0;
}

/****************************
 * Set the variable to branch mapping for the output weights ROOT file. Each
 * variable that goes in the ROOT file needs a corresponding variable in the
 * XSEvent structure. Below is a small set to get started.
 */
void setOutputBranches(TTree& t, XSEvent& xs_evt, bool use_truth_tree)
{
  std::cout << TAG << "Setting output branches." << std::endl;

  t.Branch("Enutrue", &xs_evt.enu_true, "enu_true/F");
  t.Branch("Enureco", &xs_evt.enu_reco, "enu_reco/F");
  t.Branch("sample", &xs_evt.sample, "sample/I");
  t.Branch("target", &xs_evt.target, "target/I");
  t.Branch("fgdtarget", &xs_evt.fgdtarget, "fgdtarget/I");
  t.Branch("reaction", &xs_evt.reaction, "reaction/I");
  t.Branch("topology", &xs_evt.topology, "topology/I");
  t.Branch("nutype", &xs_evt.nutype, "nutype/I");
  t.Branch("layer", &xs_evt.layer, "layer/I");
  t.Branch("Q2True", &xs_evt.q2_true, "Q2True/F");
  t.Branch("WTrue", &xs_evt.w_true, "WTrue/F");

  //Info to go into truth tree
  if(use_truth_tree){
    t.Branch("Pmutrue", &xs_evt.pmu_true, "pmu_true/F");
    t.Branch("CTHmutrue", &xs_evt.cosmu_true, "cosmu_true/F");
    t.Branch("Ppitrue", &xs_evt.ppi_true, "ppi_true/F");
    t.Branch("CTHpitrue", &xs_evt.cospi_true, "cospi_true/F");
  }
  //Reco information
  else{
    t.Branch("Pmureco", &xs_evt.pmu_reco, "pmu_reco/F");
    t.Branch("CTHmureco", &xs_evt.cosmu_reco, "cosmu_reco/F");
    t.Branch("Ppireco", &xs_evt.ppi_reco, "ppi_reco/F");
    t.Branch("CTHpireco", &xs_evt.cospi_reco, "cospi_reco/F");
    t.Branch("Pmutrue", &xs_evt.pmu_seltrue, "pmu_true/F");
    t.Branch("CTHmutrue", &xs_evt.cosmu_seltrue, "cosmu_true/F");
    t.Branch("Ppitrue", &xs_evt.ppi_seltrue, "ppi_true/F");
    t.Branch("CTHpitrue", &xs_evt.cospi_seltrue, "cospi_true/F");
    t.Branch("selection", &xs_evt.selection, "selection/I");
  }
      
  if(use_truth_tree) t.Branch("weight_nom", &xs_evt.weight_truth, "weight_nom/F");
  else t.Branch("weight_nom", &xs_evt.weight_nom, "weight_nom/F");
}

/****************************
 * Set the variable to branch mapping for the input HighLand files. Each
 * variable needs a corresponding variable in the XSEvent structure.
 */
void setInputBranches(TChain& c, XSEvent& xs_evt, bool use_truth_tree, int nsel)
{
  std::cout << TAG << "Setting input branches." << std::endl;
  c.ResetBranchAddresses();
  //Names of these categories depends on selection run
  if (nsel==0)//FGD1
    {
      c.SetBranchAddress("topology", &xs_evt.topology);
      c.SetBranchAddress("reducedreaction", &xs_evt.reaction);
      c.SetBranchAddress("moleculartarget", &xs_evt.target);
      c.SetBranchAddress("fgdtargetCC1Pi", &xs_evt.fgdtarget);
      c.SetBranchAddress("nutype", &xs_evt.nutype);
    }
  else if (nsel==1)//FGD2
    {
      c.SetBranchAddress("fgd2topology", &xs_evt.topology);
      c.SetBranchAddress("fgd2reducedreaction", &xs_evt.reaction);
      c.SetBranchAddress("fgd2moleculartarget", &xs_evt.target);
      c.SetBranchAddress("fgd2fgdtargetCC1Pi", &xs_evt.fgdtarget);
      c.SetBranchAddress("fgd2nutype", &xs_evt.nutype);
    }
  else std::cout << ERR << "Incorrect selection number provided in setInputBranches." << std::endl;
    
  c.SetBranchAddress("run", &xs_evt.run);
  c.SetBranchAddress("subrun", &xs_evt.subrun);
  c.SetBranchAddress("truelepton_mom", &xs_evt.pmu_true);
  c.SetBranchAddress("truelepton_costheta", &xs_evt.cosmu_true);
  c.SetBranchAddress("truepi_mom", &xs_evt.ppi_true);
  c.SetBranchAddress("truepi_costheta", &xs_evt.cospi_true);
  c.SetBranchAddress("nu_trueE", &xs_evt.enu_true);
  c.SetBranchAddress("RooVtxIndex", &xs_evt.input_RooVtxIndex);
  c.SetBranchAddress("RooVtxEntry", &xs_evt.input_RooVtxEntry);
  c.SetBranchAddress("TruthVertexID", &xs_evt.TruthVertexID);
  c.SetBranchAddress("true_Q2", &xs_evt.q2_true);
  c.SetBranchAddress("true_W", &xs_evt.w_true);
  
  if(use_truth_tree){
    c.SetBranchAddress("fgd2layertruth", &xs_evt.layer);
    c.SetBranchAddress("weight", &xs_evt.weight_truth);
  }
  else
    {
      c.SetBranchAddress("fgd2layer", &xs_evt.layer);
      c.SetBranchAddress("selmu_mom", &xs_evt.pmu_reco);
      c.SetBranchAddress("selmu_truemom", &xs_evt.pmu_seltrue);
      c.SetBranchAddress("selmu_costheta", &xs_evt.cosmu_reco);
      c.SetBranchAddress("selmu_truedir", &xs_evt.dirmu_seltrue);
      c.SetBranchAddress("selpi_mom2", &xs_evt.ppi_reco);
      c.SetBranchAddress("selpi_truemom", &xs_evt.ppi_seltrue);
      c.SetBranchAddress("selpi_costheta2", &xs_evt.cospi_reco);
      c.SetBranchAddress("selpi_truecostheta", &xs_evt.cospi_seltrue);
      c.SetBranchAddress("selmu_nuErecQE", &xs_evt.enu_reco);
      c.SetBranchAddress("nu_truedir", &xs_evt.dirnu_true);
      c.SetBranchAddress("weight_corr_total", &xs_evt.weight_nom);
    }

}


std::string findFileExtension(const std::string& file)
{
  std::size_t dot_pos = file.rfind('.');
  if(dot_pos != std::string::npos)
    {
      return file.substr(dot_pos);
    }
  else
    return "";
}

std::vector<std::string> getInputRootFiles(const std::vector<std::string>& file_list)
{
  std::vector<std::string> temp_list;
  for(const auto& file : file_list)
    {
      const auto file_ext = findFileExtension(file);
      if(file_ext == ".txt")
        {
	  std::ifstream fin(file, std::ios::in);
	  if(!fin.is_open())
            {
	      std::cerr << ERR << "Failed to open " << file << std::endl;
	      std::cerr << ERR << "Skipping file..." << std::endl;
            }
	  else
            {
	      std::cout << TAG << "Reading " << file << " for input ROOT files." << std::endl;
	      std::string name;
	      while(std::getline(fin, name))
                {
		  temp_list.emplace_back(name);
                }
            }
        }
      else if(file_ext == ".root")
        {
	  temp_list.emplace_back(file);
        }
      else
        {
	  std::cout << "Unsupported input file type: " << file << std::endl;
        }
    }

  return temp_list;
}

