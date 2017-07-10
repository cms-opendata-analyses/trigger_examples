// -*- C++ -*-
//
// Package:    TriggerInfoAnalyzer
// Class:      TriggerInfoAnalyzer
// 
/**\class TriggerInfoAnalyzer TriggerInfoAnalyzer.cc TriggerInfo/TriggerInfoAnalyzer/src/TriggerInfoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Edgar Carrera
//         Created:  Mon Jul  3 15:59:18 CEST 2017
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//Following the HLTEventAnalyzerAOD.h, 
//include the following headers:
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"

//Also include headers from  HLTEventAnalyzerAOD.cc
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"
#include <cassert>


//
// class declaration
//

class TriggerInfoAnalyzer : public edm::EDAnalyzer {
   public:
      explicit TriggerInfoAnalyzer(const edm::ParameterSet&);
      ~TriggerInfoAnalyzer();

      virtual void beginRun(edm::Run const&, edm::EventSetup const&);
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void analyzeTrigger(const edm::Event&, const edm::EventSetup&, const std::string& triggerName);
      //the follwing are not being used here
      virtual void beginJob() ;
      virtual void endJob() ;
      virtual void endRun(edm::Run const&, edm::EventSetup const&);
      virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
     
      // from HLTEventAnalyzerAOD.h
      /// module config parameters
      std::string   processName_;
      std::string   triggerName_;
      std::string   datasetName_;
      edm::InputTag triggerResultsTag_;
      edm::InputTag triggerEventTag_;


      // additional class data memebers
      // these are actually the containers where we will store
      // the trigger information
      edm::Handle<edm::TriggerResults>   triggerResultsHandle_;
      edm::Handle<trigger::TriggerEvent> triggerEventHandle_;
      HLTConfigProvider hltConfig_;



      // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
// Notice that here, using the parameter set tool, 
// you need to point the code to
// the right branch (in the EDM root files) where the trigger information
// is stored.  

// Also, at configuration time
// you will need to point to the appropiate triggers
// you want to look at. Alternatively (is also shown below), you
// could select the trigger names dynamically; for example getting them
// from the HLTConfigProvider.

// To start out, you need to define a processName, which is the name of
// the CMSSW computing process that originally wrote the products in the root
// file. Originally, this is always "HLT", by default.  
// In triggerName, you can
// use wildcards, which will be described later.
// As for the InputTags, these shall match the name of the ROOT branches
// where the information is stored.  This was essentially fixed and will
// most likely be the same always. 

//This should match your configuration python file
TriggerInfoAnalyzer::TriggerInfoAnalyzer(const edm::ParameterSet& ps):
processName_(ps.getParameter<std::string>("processName")),
triggerName_(ps.getParameter<std::string>("triggerName")),
datasetName_(ps.getParameter<std::string>("datasetName")),
triggerResultsTag_(ps.getParameter<edm::InputTag>("triggerResults")),
triggerEventTag_(ps.getParameter<edm::InputTag>("triggerEvent"))
{
   //now do what ever initialization is needed
  using namespace std;
  using namespace edm;
  
  //Print the configuration just to check
  cout << "Here is the information passed to the constructor:" <<endl;
  cout << "HLTEventAnalyzerAOD configuration: " << endl
       << "   ProcessName = " << processName_ << endl
       << "   TriggerName = " << triggerName_ << endl
       << "   DataSetName = " << datasetName_ << endl
       << "   TriggerResultsTag = " << triggerResultsTag_.encode() << endl
       << "   TriggerEventTag = " << triggerEventTag_.encode() << endl;

}


TriggerInfoAnalyzer::~TriggerInfoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//




// ------------ method called when starting to processes a run  ------------
void TriggerInfoAnalyzer::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
//--------------------------------------------------------------------------
{
  using namespace std;
  using namespace edm;


  //If the hltConfig can be initialized, then the below is an example of
  //how to extract the config information for the trigger from the 
  //so-called provenance.

  // The trigger configuration can change from 
  // run to run (during the run is the same),
  // so it needs to be called here.

  ///   "init" return value indicates whether intitialisation has succeeded
  ///   "changed" parameter indicates whether the config has actually changed

  bool changed(true);
  if (hltConfig_.init(iRun,iSetup,processName_,changed)) {
    if (changed) {
      // check if trigger name in (new) config
      if (triggerName_!="@") { // "@" means: analyze all triggers in config
	const unsigned int n(hltConfig_.size());
	const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName_));
	if (triggerIndex>=n) {
	  cout << "HLTEventAnalyzerAOD::analyze:"
	       << " TriggerName " << triggerName_ 
	       << " not available in (new) config!" << endl;
	  cout << "Available TriggerNames are: " << endl;
	  hltConfig_.dump("Triggers");
	}
      }

      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      //EXAMPLE: How to dump information from the Provenance
      // Uncomment the following lines if needed.
      //the hltConfig has many accessors
      //For details see the header of the class,HLTConfigProvider.h
      //To check the example, uncomment the lines below
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      //hltConfig_.dump("Streams");
      //hltConfig_.dump("Datasets"); 
      //hltConfig_.dump("PrescaleTable");
      //hltConfig_.dump("ProcessPSet");
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    }


  
  } else {
    cout << "HLTEventAnalyzerAOD::analyze:"
	 << " config extraction failure with process name "
	 << processName_ << endl;
  }
  

}//------------------- beginRun()





// ------------ method called for each event  ------------------------------
// As with any EDAnalyzer, this method is the heart of the analysis
void TriggerInfoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//--------------------------------------------------------------------------
{
   using namespace edm;
   using namespace std;

   // Get event products: 
   // In the following, the code is trying to access the information 
   // from the ROOT files and point the containers (that we created), 
   // namely triggerResultsHandle_ and triggerEVentHandle_, 
   // to the correct "address", given at configuration time 
   // and assigned to triggerResultsTag_ and triggerEventTag_
 
   // After that, a simple sanity check is done.
 
   iEvent.getByLabel(triggerResultsTag_,triggerResultsHandle_);
   if (!triggerResultsHandle_.isValid()) {
     cout << "HLTEventAnalyzerAOD::analyze: Error in getting TriggerResults product from Event!" << endl;
     return;
   }
   iEvent.getByLabel(triggerEventTag_,triggerEventHandle_);
   if (!triggerEventHandle_.isValid()) {
     cout << "HLTEventAnalyzerAOD::analyze: Error in getting TriggerEvent product from Event!" << endl;
     return;
   }
   // sanity check
   assert(triggerResultsHandle_->size()==hltConfig_.size());
   

   //The following two examples should be used separately or somehow
   //combine them according the analysis needs
   
   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   //EXAMPLE: analyze this event for triggers requested at config time
   //Uncomment the lines below
   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   //if (triggerName_=="@") {
   //  const unsigned int n(hltConfig_.size());
   //  for (unsigned int i=0; i!=n; ++i) {
   //    analyzeTrigger(iEvent,iSetup,hltConfig_.triggerName(i));
   //  }
   //} else {
   //  analyzeTrigger(iEvent,iSetup,triggerName_);
   //}
   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   //EXAMPLE: analyze this event for triggers that make up the dataset
   //selected at configuration time.
   //Notice that we can find out which triggers go into any stream
   //and/or dataset using the hltConfig_
   //Uncomment the lines below
   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   //Get all the trigger populating the datasetName_ dataset
   const vector<string> triggerNamesInDS = hltConfig_.datasetContent(datasetName_);
   for (unsigned i = 0; i < triggerNamesInDS.size(); i++) {
     analyzeTrigger(iEvent,iSetup,triggerNamesInDS[i]);
   }
   //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  return;

}//---------------------------analyze()




//---------------------------Actual trigger analysis-------------
void TriggerInfoAnalyzer::analyzeTrigger(const edm::Event& iEvent, const edm::EventSetup& iSetup, const std::string& triggerName) 
//-----------------------------------------------------------------
{

  using namespace std;
  using namespace edm;
  using namespace reco;
  using namespace trigger;

  cout<<"Currently analyzing trigger "<<triggerName<<endl;

  //Check the current configuration to see how many total triggers there are
  const unsigned int n(hltConfig_.size());
  //Get the trigger index for the current trigger
  const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  //check that the trigger in the event and in the configuration agree
  assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));
  // abort on invalid trigger name
  if (triggerIndex>=n) {
    cout << "HLTEventAnalyzerAOD::analyzeTrigger: path "
	 << triggerName << " - not found!" << endl;
    return;
  }

  
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // EXAMPLE: L1 and HLT prescale values via (L1) EventSetup
  // Current (default) prescale set index - to be taken from L1GtUtil via Event.
  // Try to get the L1 and HLT prescale values that were actually used 
  // for this event.
  // This example needs the conditions stored in the Global Tag,
  // which is some sort of snapshot of the by-then current detector
  // conditions.  They need to be extracted from the database
  // and for that the "Frontier Conditions" lines need to
  // be added in the python configuration file along with the 
  // name for the global tag.
  // This can make the job very slow at the very begining....
  //Uncomment the lines below
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  const std::pair<int,int> prescales(hltConfig_.prescaleValues(iEvent,iSetup,triggerName));
  cout << "HLTEventAnalyzerAOD::analyzeTrigger: path "
      << triggerName << " [" << triggerIndex << "] "
      << "prescales L1T,HLT: " << prescales.first << "," << prescales.second
      << endl;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // EXAMPLE: Find out if the trigger was active, accepted, or in error.
  // We could also find out whether the trigger was active (wasrun), 
  // if it accepted the event (accept) or if it gave an error (error).
  // Results from TriggerResults product
  //Uncomment the lines below
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  cout << " Trigger path status:"
       << " WasRun=" << triggerResultsHandle_->wasrun(triggerIndex)
       << " Accept=" << triggerResultsHandle_->accept(triggerIndex)
       << " Error =" << triggerResultsHandle_->error(triggerIndex)
       << endl;
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // EXAMPLE: Which module was last run within the trigger path
  // Additionally, one can find which module was last run within the trigger,
  // i.e., which module was responsible for stopping the trigger path.
  // One could find the list of modules in a given trigger from the
  // HLTConfigProvider as explained somewhere else in this code.
  // Get index (slot position) of module giving the decision of the path
  // as described in "DataFormats/Common/interface/HLTGlobalStatus.h"
  //Uncomment the lines below
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  const unsigned int m(hltConfig_.size(triggerIndex));
  const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));
  const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
  cout << " Last active module - label/type: "
       << moduleLabels[moduleIndex] << "/" << hltConfig_.moduleType(moduleLabels[moduleIndex])
       << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
       << endl;
  assert (moduleIndex<m);
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  // EXAMPLE: How to extract some additional information 
  // from TriggerEvent product (L3 information) 
  // Must look only for  modules actually run in this path for this event,
  // so you loop only up to moduleIndex obtained above.
  // Uncomment the lines below
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 //  for (unsigned int j=0; j<=moduleIndex; ++j) {
//     const string& moduleLabel(moduleLabels[j]);
//     const string  moduleType(hltConfig_.moduleType(moduleLabel));
//     // check whether the module is packed up in TriggerEvent product
//     // find index of filter in data-member vector from filter tag
//     // Look at DataFormats/HLTReco/interface/TriggerEvent.h
//     const unsigned int filterIndex(triggerEventHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
//     if (filterIndex<triggerEventHandle_->sizeFilters()) {
//       cout << " 'L3' filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl;
//       const Vids& VIDS (triggerEventHandle_->filterIds(filterIndex));
//       const Keys& KEYS(triggerEventHandle_->filterKeys(filterIndex));
//       const size_type nI(VIDS.size());
//       const size_type nK(KEYS.size());
//       assert(nI==nK);
//       const size_type n(max(nI,nK));
//       cout << "   " << n  << " accepted 'L3' objects found: " << endl;
//       const TriggerObjectCollection& TOC(triggerEventHandle_->getObjects());
//       for (size_type i=0; i!=n; ++i) {
// 	const TriggerObject& TO(TOC[KEYS[i]]);
// 	cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "
// 	     << TO.id() << " " << TO.pt() << " " << TO.eta() << " " 
// 	     << TO.phi() << " " << TO.mass()
// 	     << endl;
//       }
//     }
//   } 
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
  return;
}//--------------------------analyzeTrigger()



// ------------ method called once each job just before starting event loop  ------------
void 
TriggerInfoAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
TriggerInfoAnalyzer::endJob() 
{
}


// ------------ method called when ending the processing of a run  ------------
void TriggerInfoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}


// ------------ method called when starting to processes a luminosity block  ------------
void 
TriggerInfoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
TriggerInfoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
TriggerInfoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(TriggerInfoAnalyzer);
