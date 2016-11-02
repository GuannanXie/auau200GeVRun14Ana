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
#include "StarClassLibrary/StLorentzVectorF.hh"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoPrescales/StPicoPrescales.h"
#include "StPicoCharmContainers/StPicoKPiX.h"

#include "StPicoKPiXAnaMaker/StKPiXAnaCuts.h"
#include "StPicoCharmMassHists.h"

//-----------------------------------------------------------------------
StPicoCharmMassHists::StPicoCharmMassHists(std::string fileBaseName, std::string prescalesDir, bool const fillTopoDistHistograms): mFillTopoDistHistograms(fillTopoDistHistograms)
{
   TH1::AddDirectory(0);

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

   if(mFillTopoDistHistograms)
   {
     float const maxDca = 0.3;
     int const nDcaBins = 300;

     float const maxDca12 = 0.1;
     int const nDca12Bins = 200;

     mCosTheta[0] = new TH3F("h3FgCosTheta", ";p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 1000, 0.9, 1.0);
     mDecayL[0] = new TH3F("h3FgdecayL", ";p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
     mDca12[0] = new TH3F("h3FgdcaDaughters", ";p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDca12Bins, 0, maxDca12);
     mKaonDca2Vtx[0] = new TH3F("h3FgkaonDca", ";p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
     mPionDca2Vtx[0] = new TH3F("h3FgpionDca", ";p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
     mXaonDca2Vtx[0] = new TH3F("h3FgxaonDca", ";p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
     mPerpDca2Vtx[0] = new TH3F("h3FgPerpDca2Vtx", ";p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 100, 0, 0.05);

     mCosTheta[1] = new TH3F("h3BgCosTheta", ";p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 1000, 0.9, 1.0);
     mDecayL[1] = new TH3F("h3BgdecayL", ";p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
     mDca12[1] = new TH3F("h3BgdcaDaughters", ";p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDca12Bins, 0, maxDca12);
     mKaonDca2Vtx[1] = new TH3F("h3BgkaonDca", ";p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
     mPionDca2Vtx[1] = new TH3F("h3BgpionDca", ";p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
     mXaonDca2Vtx[1] = new TH3F("h3BgxaonDca", ";p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, nDcaBins,0,maxDca);
     mPerpDca2Vtx[1] = new TH3F("h3BgPerpDca2Vtx", ";p_{T} (GeV/c);centrality", 150, 0, 15, 9, 0, 9, 100, 0, 0.05);
   }
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

void StPicoCharmMassHists::fillTopoDistHistograms(StPicoKPiX const&  kpx, bool fg, int const centrality, kPiXAnaCuts::TopologicalCuts const& cuts)
{
   if(!mFillTopoDistHistograms) return;

   StLorentzVectorF fourMom = kpx.fourMom(cuts.xMassHypothesis);
   if (fourMom.m() <  cuts.minMass || fourMom.m() > cuts.maxMass) return;

   int iArr = fg? 0 : 1;
   float const pt = fourMom.perp();

   int ptIndex = -1;
   for (size_t i = 0; i<cuts.ptBinsEdge.size(); ++i)
   {
     if (pt >= cuts.ptBinsEdge[i] && pt < cuts.ptBinsEdge[i + 1])
     {
       ptIndex = i;
       break;
     }
   }
   if(ptIndex < 0) return;

   //CosTheta
   if (kpx.kaonDca() > cuts.kDca[ptIndex] &&
       kpx.pionDca() > cuts.pDca[ptIndex] &&
       kpx.xaonDca() > cuts.xDca[ptIndex] &&
       kpx.dcaDaughters() < cuts.dcaDaughters[ptIndex] &&
       kpx.decayLength() > cuts.decayLength[ptIndex] &&
       //std::cos(kpx.pointingAngle()) > cuts.cosTheta[ptIndex] &&
       kpx.perpDcaToVtx() < cuts.dcaV0ToPv[ptIndex])
     mCosTheta[iArr]->Fill(pt, centrality, std::cos(kpx.pointingAngle()));

   //Decay Length
   if (kpx.kaonDca() > cuts.kDca[ptIndex] &&
       kpx.pionDca() > cuts.pDca[ptIndex] &&
       kpx.xaonDca() > cuts.xDca[ptIndex] &&
       kpx.dcaDaughters() < cuts.dcaDaughters[ptIndex] &&
       // kpx.decayLength() > cuts.decayLength[ptIndex] &&
       std::cos(kpx.pointingAngle()) > cuts.cosTheta[ptIndex] &&
       kpx.perpDcaToVtx() < cuts.dcaV0ToPv[ptIndex])
     mDecayL[iArr]->Fill(pt, centrality, kpx.decayLength());

   //Dca12
   if (kpx.kaonDca() > cuts.kDca[ptIndex] &&
       kpx.pionDca() > cuts.pDca[ptIndex] &&
       kpx.xaonDca() > cuts.xDca[ptIndex] &&
       // kpx.dcaDaughters() < cuts.dcaDaughters[ptIndex] &&
       kpx.decayLength() > cuts.decayLength[ptIndex] &&
       std::cos(kpx.pointingAngle()) > cuts.cosTheta[ptIndex] &&
       kpx.perpDcaToVtx() < cuts.dcaV0ToPv[ptIndex])
     mDca12[iArr]->Fill(pt, centrality, kpx.dcaDaughters());

   //kaonDca
   if (//kpx.kaonDca() > cuts.kDca[ptIndex] &&
       kpx.pionDca() > cuts.pDca[ptIndex] &&
       kpx.xaonDca() > cuts.xDca[ptIndex] &&
       kpx.dcaDaughters() < cuts.dcaDaughters[ptIndex] &&
       kpx.decayLength() > cuts.decayLength[ptIndex] &&
       std::cos(kpx.pointingAngle()) > cuts.cosTheta[ptIndex] &&
       kpx.perpDcaToVtx() < cuts.dcaV0ToPv[ptIndex])
     mKaonDca2Vtx[iArr]->Fill(pt, centrality, kpx.kaonDca());

   //pionDca
   if (kpx.kaonDca() > cuts.kDca[ptIndex] &&
       // kpx.pionDca() > cuts.pDca[ptIndex] &&
       kpx.xaonDca() > cuts.xDca[ptIndex] &&
       kpx.dcaDaughters() < cuts.dcaDaughters[ptIndex] &&
       kpx.decayLength() > cuts.decayLength[ptIndex] &&
       std::cos(kpx.pointingAngle()) > cuts.cosTheta[ptIndex] &&
       kpx.perpDcaToVtx() < cuts.dcaV0ToPv[ptIndex])
     mPionDca2Vtx[iArr]->Fill(pt, centrality, kpx.pionDca());

   //xaonDca
   if (kpx.kaonDca() > cuts.kDca[ptIndex] &&
       kpx.pionDca() > cuts.pDca[ptIndex] &&
       // kpx.xaonDca() > cuts.xDca[ptIndex] &&
       kpx.dcaDaughters() < cuts.dcaDaughters[ptIndex] &&
       kpx.decayLength() > cuts.decayLength[ptIndex] &&
       std::cos(kpx.pointingAngle()) > cuts.cosTheta[ptIndex] &&
       kpx.perpDcaToVtx() < cuts.dcaV0ToPv[ptIndex])
     mXaonDca2Vtx[iArr]->Fill(pt, centrality, kpx.xaonDca());

   //perpDcaToVtx
   if (kpx.kaonDca() > cuts.kDca[ptIndex] &&
       kpx.pionDca() > cuts.pDca[ptIndex] &&
       kpx.xaonDca() > cuts.xDca[ptIndex] &&
       kpx.dcaDaughters() < cuts.dcaDaughters[ptIndex] &&
       kpx.decayLength() > cuts.decayLength[ptIndex] &&
       std::cos(kpx.pointingAngle()) > cuts.cosTheta[ptIndex]
       // kpx.perpDcaToVtx() < cuts.dcaV0ToPv[ptIndex]
      )
     mPerpDca2Vtx[iArr]->Fill(pt, centrality, kpx.perpDcaToVtx());
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

   if(mFillTopoDistHistograms)
   {
     for(int i=0; i<2; ++i)
     {
       mCosTheta[i]->Write();
       mDecayL[i]->Write();
       mDca12[i]->Write();
       mKaonDca2Vtx[i]->Write();
       mPionDca2Vtx[i]->Write();
       mXaonDca2Vtx[i]->Write();
       mPerpDca2Vtx[i]->Write();
     }
   }

   mOutFile->Close();
   mOutFile->Delete();
}
