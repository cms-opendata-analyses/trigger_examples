import FWCore.ParameterSet.Config as cms

process = cms.Process("TriggerInfo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#    'root://eospublic.cern.ch//eos/opendata/cms/Run2010B/Mu/AOD/Apr21ReReco-v1/0000/00459D48-EB70-E011-AF09-90E6BA19A252.root'
    'file:0048DD36-6570-E011-A43C-485B39800BA2.root' 
    )
)

#needed to get the actual prescale values used from the global tag
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'GR_R_42_V25::All'

#configure the analyzer
process.gettriggerinfo = cms.EDAnalyzer('TriggerInfoAnalyzer',
                              processName = cms.string("HLT"),
                              triggerName = cms.string("@"),         
                              datasetName = cms.string("Mu"),           
                              triggerResults = cms.InputTag("TriggerResults","","HLT"),
                              triggerEvent   = cms.InputTag("hltTriggerSummaryAOD","","HLT")                             
                              )


process.triggerinfo = cms.Path(process.gettriggerinfo)
process.schedule = cms.Schedule(process.triggerinfo)
