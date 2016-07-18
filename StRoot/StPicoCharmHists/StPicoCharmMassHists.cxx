/* **************************************************
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include <cmath>

#include "St_base/StMessMgr.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TFile.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoPrescales/StPicoPrescales.h"
#include "StPicoCharmContainers/StPicoKPiX.h"

#include "StPicoCharmMassHists.h"

//-----------------------------------------------------------------------
StPicoCharmMassHists::StPicoCharmMassHists(std::string fileBaseName, std::string prescalesDir)
{
   mPrescales = new StPicoPrescales(prescalesDir);

   mOutFile = new TFile(Form("%s.hists.root", fileBaseName.c_str()), "RECREATE");

   int nRuns = mPrescales->numberOfRuns();
   TH1::AddDirectory(0);
   TH1::SetDefaultSumw2();
   mh1TotalEventsInRun         = new TH1F("mh1TotalEventsInRun", "totalEventsInRun;runIndex;totalEventsInRun", nRuns + 1, 0, nRuns + 1);
   mh1TotalEventsInRunBeforeCut = new TH1F("mh1TotalEventsInRunBeforeCut", "totalEventsInRun;runIndex;totalEventsInRun", nRuns + 1, 0, nRuns + 1);

   //add centrality
   mh1Cent          = new TH1F("mh1Cent", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5);
   mh1CentWg        = new TH1F("mh1CentWg", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5);
   mh1gRefmultCor   = new TH1F("mh1gRefmultCor", "gRefmultCor;gRefmult;Counts", 700, 0, 700);
   mh1gRefmultCorWg = new TH1F("mh1gRefmultCorWg", "gRefmultCorWg;gRefmultCorWg;Counts", 700, 0, 700);
   mh2CentVz        = new TH2F("mh2CentVz", "CentralityVsVz;cent;Vz", 10, -1.5, 8.5, 200, -10, 10);
   mh2CentVzWg      = new TH2F("mh2CentVzWg", "CentralityVsVzWg;cent;Vz", 10, -1.5, 8.5, 200, -10, 10);

   mh3InvariantMassVsPtVsCent    = new TH3F("mh3InvariantMassVsPtVsCent", "invariantMassVsPtVsCent;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})", 150, 0, 15, 10, -1.5, 8.5, 100, 1.6, 2.6);
   mh3InvariantMassVsPtVsCentBkg = new TH3F("mh3InvariantMassVsPtVsCentBkg", "invariantMassVsPtVsCentBkg;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})", 150, 0, 15, 10, -1.5, 8.5, 100, 1.6, 2.6);
}

StPicoCharmMassHists::~StPicoCharmMassHists()
{
   delete mPrescales;
   // note that histograms are owned by mOutFile. They will be destructed
   // when the file is closed.
}

void StPicoCharmMassHists::addEvent(StPicoEvent const& picoEvent)
{
   int runIndex = mPrescales->runIndex(picoEvent.runId());
   mh1TotalEventsInRun->Fill(runIndex);
}

void StPicoCharmMassHists::addEventBeforeCut(StPicoEvent const& picoEvent)
{
   int runIndex = mPrescales->runIndex(picoEvent.runId());
   mh1TotalEventsInRunBeforeCut->Fill(runIndex);
}

void StPicoCharmMassHists::addCent(double const refmultCor, int const centrality, double const reweight, float const vz)
{
   mh1gRefmultCor->Fill(refmultCor);
   mh1gRefmultCorWg->Fill(refmultCor, reweight);
   mh1Cent->Fill(centrality);
   mh1CentWg->Fill(centrality, reweight);
   mh2CentVz->Fill(centrality, vz);
   mh2CentVzWg->Fill(centrality, vz, reweight);
}

void StPicoCharmMassHists::addKPiX(StLorentzVectorF const& fMom, bool fg, int centrality, const double reweight)
{
   if (fg) mh3InvariantMassVsPtVsCent->Fill(fMom.perp(), centrality, fMom.m(), reweight);
   else mh3InvariantMassVsPtVsCentBkg->Fill(fMom.perp(), centrality, fMom.m(), reweight);
}

void StPicoCharmMassHists::closeFile()
{
   mOutFile->cd();

   mh1TotalEventsInRun->Write();
   mh1TotalEventsInRunBeforeCut->Write();

   //centrality
   mh1Cent->Write();
   mh1CentWg->Write();
   mh1gRefmultCor->Write();
   mh1gRefmultCorWg->Write();
   mh2CentVz->Write();
   mh2CentVzWg->Write();

   mh3InvariantMassVsPtVsCent->Write();
   mh3InvariantMassVsPtVsCentBkg->Write();

   mOutFile->Close();
   mOutFile->Delete();
}
