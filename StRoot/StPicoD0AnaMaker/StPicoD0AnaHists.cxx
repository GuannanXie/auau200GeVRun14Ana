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
#include "TString.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoPrescales/StPicoPrescales.h"
#include "StPicoCharmContainers/StKaonPion.h"
#include "StAnaCuts.h"
#include "TNtuple.h"
#include "THnSparse.h"
#include "THn.h"

#include "StPicoD0AnaHists.h"

//-----------------------------------------------------------------------
StPicoD0AnaHists::StPicoD0AnaHists(TString fileBaseName, bool fillQaHists) : mFillQaHists(fillQaHists), mPrescales(NULL), mOutFile(NULL),
   mh2InvariantMassVsPt(NULL), mh2InvariantMassVsPtLike(NULL), mh2InvariantMassVsPtTof(NULL), mh2InvariantMassVsPtTofLike(NULL),
   mh1Cent(NULL), mh1CentWg(NULL), mh1gRefmultCor(NULL), mh1gRefmultCorWg(NULL), mh2CentVz(NULL), mh2CentVzWg(NULL), mh3InvariantMassVsPtVsCent(NULL), mh3InvariantMassVsPtVsCentLike(NULL), mh3InvariantMassVsPtVsCentTof(NULL), mh3InvariantMassVsPtVsCentTofLike(NULL)
   // mh2Tpc1PtCent(NULL),  mh2Tpc1PhiVz(NULL), mh2HFT1PtCent(NULL),  mh2HFT1PhiVz(NULL),  mh3DcaXyPtCent(NULL), mh3DcaZPtCent(NULL),
{
   mPrescales = new StPicoPrescales(anaCuts::prescalesFilesDirectoryName);

   mOutFile = new TFile(Form("%s.hists.root", fileBaseName.Data()), "RECREATE");

   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < anaCuts::nEtasDca; iEta++)
      {
         for (int iVz = 0; iVz < anaCuts::nVzsDca; iVz++)
         {
            for (int iCent = 0; iCent < anaCuts::nCentsDca; iCent++)
            {
               for (int iPhi = 0; iPhi < anaCuts::nPhisDca; iPhi++)
               {
                  mh3DcaXyZPtCentPartEtaVzPhi[iParticle][iEta][iVz][iCent][iPhi] = NULL;
               }
            }
         }
      }
   }

   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < anaCuts::nEtasRatio; iEta++)
      {
         for (int iVz = 0; iVz < anaCuts::nVzsRatio; iVz++)
         {
            for (int iPhi = 0; iPhi < anaCuts::nPhisRatio; iPhi++)
            {
               mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi] = NULL;
               mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi] = NULL;
            }
         }
      }
   }


   int nRuns = mPrescales->numberOfRuns();
   TH1::SetDefaultSumw2();
   mh1TotalEventsInRun         = new TH1F("mh1TotalEventsInRun", "totalEventsInRun;runIndex;totalEventsInRun", nRuns + 1, 0, nRuns + 1);
   mh1TotalEventsInRunBeforeCut = new TH1F("mh1TotalEventsInRunBeforeCut", "totalEventsInRun;runIndex;totalEventsInRun", nRuns + 1, 0, nRuns + 1);
   mh2InvariantMassVsPt        = new TH2F("mh2InvariantMassVsPt", "invariantMassVsPt;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})", 120, 0, 12, 50, 1.6, 2.1);
   mh2InvariantMassVsPtLike    = new TH2F("mh2InvariantMassVsPtLike", "invariantMassVsPtLike;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})", 120, 0, 12, 50, 1.6, 2.1);
   mh2InvariantMassVsPtTof     = new TH2F("mh2InvariantMassVsPtTof", "invariantMassVsPtTof;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})", 120, 0, 12, 50, 1.6, 2.1);
   mh2InvariantMassVsPtTofLike = new TH2F("mh2InvariantMassVsPtTofLike", "invariantMassVsPtTofLike;p_{T}(K#pi)(GeV/c);m_{K#pi}(GeV/c^{2})", 120, 0, 12, 50, 1.6, 2.1);
   //add centrality
   mh1Cent         = new TH1F("mh1Cent", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5);
   mh1CentWg         = new TH1F("mh1CentWg", "EventsVsCentrality;cent;Counts", 10, -1.5, 8.5);
   mh1gRefmultCor  = new TH1F("mh1gRefmultCor", "gRefmultCor;gRefmult;Counts", 700, 0, 700);
   mh1gRefmultCorWg  = new TH1F("mh1gRefmultCorWg", "gRefmultCorWg;gRefmultCorWg;Counts", 700, 0, 700);
   mh2CentVz         = new TH2F("mh2CentVz", "CentralityVsVz;cent;Vz", 10, -1.5, 8.5, 200, -10, 10);
   mh2CentVzWg         = new TH2F("mh2CentVzWg", "CentralityVsVzWg;cent;Vz", 10, -1.5, 8.5, 200, -10, 10);

   mh3InvariantMassVsPtVsCent        = new TH3F("mh3InvariantMassVsPtVsCent", "invariantMassVsPtVsCent;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})", 120, 0, 12, 10, -1.5, 8.5, 50, 1.6, 2.1);
   mh3InvariantMassVsPtVsCentLike    = new TH3F("mh3InvariantMassVsPtVsCentLike", "invariantMassVsPtVsCentLike;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})", 120, 0, 12, 10, -1.5, 8.5, 50, 1.6, 2.1);
   mh3InvariantMassVsPtVsCentTof     = new TH3F("mh3InvariantMassVsPtVsCentTof", "invariantMassVsPtVsCentTof;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})", 120, 0, 12, 10, -1.5, 8.5, 50, 1.6, 2.1);
   mh3InvariantMassVsPtVsCentTofLike = new TH3F("mh3InvariantMassVsPtVsCentTofLike", "invariantMassVsPtVsCentTofLike;p_{T}(K#pi)(GeV/c);Cent;m_{K#pi}(GeV/c^{2})", 120, 0, 12, 10, -1.5, 8.5, 50, 1.6, 2.1);

   /******************************************************************************************/
   /*             NOTE: All histograms below will not be defined if mFillQaHists is not true */
   /******************************************************************************************/

   if (!mFillQaHists) return;

   //Add some HFT ratio plots
   mh2Tpc1PtCent  = new TH2F("mh2Tpc1PtCent", "Tpc tacks;p_{T}(GeV/c);cent", 120, 0, 12, 10, -1.5, 8.5); //Dca 1.5cm
   mh2HFT1PtCent  = new TH2F("mh2HFT1PtCent", "HFT tacks;p_{T}(GeV/c);cent", 120, 0, 12, 10, -1.5, 8.5); //Dca 1.5cm
   mh2Tpc1PhiVz  = new TH2F("mh2Tpc1PhiVz", "Tpc tacks;#Phi;Vz", 100, -3.1415, 3.1415, 20, -10, 10); //Dca 1.5cm
   mh2HFT1PhiVz  = new TH2F("mh2HFT1PhiVz", "HFT tacks;#Phi;Vz", 100, -3.1415, 3.1415, 20, -10, 10); //Dca 1.5cm

   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < anaCuts::nEtasRatio; iEta++)
      {
         for (int iVz = 0; iVz < anaCuts::nVzsRatio; iVz++)
         {
            for (int iPhi = 0; iPhi < anaCuts::nPhisRatio; iPhi++)
            {
               mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]  = new TH2F(Form("mh2Tpc1PtCentPartEtaVzPhi_%d_%d_%d_%d", iParticle, iEta, iVz, iPhi), Form("mh2Tpc1PtCent_%s_Eta%2.1f_Vz%2.1f_Phi%2.1f;p_{T}(GeV/c);cent", anaCuts::ParticleName[iParticle], anaCuts::EtaEdgeRatio[iEta], anaCuts::VzEdgeRatio[iVz], anaCuts::PhiEdgeRatio[iPhi]), anaCuts::nPtsRatio, anaCuts::PtEdgeRatio, anaCuts::nCentsRatio, anaCuts::CentEdgeRatio); //Dca 1.cm
               mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]  = new TH2F(Form("mh2HFT1PtCentPartEtaVzPhi_%d_%d_%d_%d", iParticle, iEta, iVz, iPhi), Form("mh2HFT1PtCent_%s_Eta%2.1f_Vz%2.1f_Phi%2.1f;p_{T}(GeV/c);cent", anaCuts::ParticleName[iParticle], anaCuts::EtaEdgeRatio[iEta], anaCuts::VzEdgeRatio[iVz], anaCuts::PhiEdgeRatio[iPhi]), anaCuts::nPtsRatio, anaCuts::PtEdgeRatio, anaCuts::nCentsRatio, anaCuts::CentEdgeRatio); //Dca 1.cm
            }
         }
      }
   }


   // Add some Dca, resolution
   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < anaCuts::nEtasDca; iEta++)
      {
         for (int iVz = 0; iVz < anaCuts::nVzsDca; iVz++)
         {
            for (int iCent = 0; iCent < anaCuts::nCentsDca; iCent++)
            {
               for (int iPhi = 0; iPhi < anaCuts::nPhisDca; iPhi++)
               {
                  mh3DcaXyZPtCentPartEtaVzPhi[iParticle][iEta][iVz][iCent][iPhi]  = new TH3F(Form("mh3DcaXyZPtCentPartEtaVzPhi_%d_%d_%d_%d_%d", iParticle, iEta, iVz, iCent, iPhi), Form("mh3DcaXyZPt_%s_Eta%2.1f_Vz%2.1f_Cent%2.1f_Phi%2.1f;p_{T}(GeV/c);DcaXy(cm);DcaZ(cm)", anaCuts::ParticleName[iParticle], anaCuts::EtaEdgeDca[iEta], anaCuts::VzEdgeDca[iVz], anaCuts::CentEdgeDca[iCent], anaCuts::PhiEdgeDca[iPhi]), anaCuts::nPtsDca, anaCuts::PtEdgeDca, anaCuts::nDcasDca, anaCuts::DcaEdgeDca, anaCuts::nDcasDca, anaCuts::DcaEdgeDca); //Dca 1.cm
               }
            }
         }
      }
   }


   mh3DcaPtCent  = new TH3F("mh3DcaPtCent", "mh3DcaPtCent;p_{T}(GeV/c);cent;Dca(cm)", 120, 0, 12, 10, -1.5, 8.5, 1000, -1, 1); //Dca 1.cm
   mh3DcaXyPtCent  = new TH3F("mh3DcaXyPtCent", "mh3DcaXyPtCent;p_{T}(GeV/c);cent;DcaXy(cm)", 120, 0, 12, 10, -1.5, 8.5, 1000, -1, 1); //Dca 1.cm
   mh3DcaZPtCent  = new TH3F("mh3DcaZPtCent", "mh3DcaZPtCent;p_{T}(GeV/c);cent;DcaZ(cm)", 120, 0, 12, 10, -1.5, 8.5, 1000, -1, 1); //Dca 1.cm

//  nt = new TNtuple("nt","nt","runnumber:dca:vz:pt:eta:phi:centrality:grefmultCor:zdcCoincidance:tofMatchFlag:hftMatchFlag");
}
StPicoD0AnaHists::~StPicoD0AnaHists()
{
   delete mPrescales;
   // note that histograms are owned by mOutFile. They will be destructed
   // when the file is closed.
}
//-----------------------------------------------------------------------
void StPicoD0AnaHists::addEvent(StPicoEvent const* const picoEvent)
{
   int runIndex = mPrescales->runIndex(picoEvent->runId());
   mh1TotalEventsInRun->Fill(runIndex);
}
//-----------------------------------------------------------------------
void StPicoD0AnaHists::addEventBeforeCut(StPicoEvent const* const picoEvent)
{
   int runIndex = mPrescales->runIndex(picoEvent->runId());
   mh1TotalEventsInRunBeforeCut->Fill(runIndex);
}
//-----------------------------------------------------------------------
void StPicoD0AnaHists::addCent(const double refmultCor, int centrality, const double reweight, const float vz)
{
   mh1gRefmultCor->Fill(refmultCor);
   mh1gRefmultCorWg->Fill(refmultCor, reweight);
   mh1Cent->Fill(centrality);
   mh1CentWg->Fill(centrality, reweight);
   mh2CentVz->Fill(centrality, vz);
   mh2CentVzWg->Fill(centrality, vz, reweight);
}
//-----------------------------------------------------------------------
void StPicoD0AnaHists::addTpcDenom1(bool IsPion, bool IsKaon, float pt, int centrality, float Eta, float Phi, float Vz, float ZdcX)
{
   if (!mFillQaHists)
   {
      LOG_ERROR << " You are trying to fill QA histograms but StPicoD0AnaHists::mFillQaHists is false -- ignoring attemp! " << endm;
   }

   int EtaIndex = getEtaIndexRatio(Eta);
   int PhiIndex = getPhiIndexRatio(Phi);
   int VzIndex = getVzIndexRatio(Vz);
   if (IsPion)
   {
      mh2Tpc1PtCentPartEtaVzPhi[0][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   if (IsKaon)
   {
      mh2Tpc1PtCentPartEtaVzPhi[1][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   mh2Tpc1PtCent->Fill(pt, centrality);
   if (fabs(Eta) < 0.1 && pt > 3.0) mh2Tpc1PhiVz->Fill(Phi, Vz);

}
//-----------------------------------------------------------------------
void StPicoD0AnaHists::addHFTNumer1(bool IsPion, bool IsKaon, float pt, int centrality, float Eta, float Phi, float Vz, float ZdcX)
{
   if (!mFillQaHists)
   {
      LOG_ERROR << " You are trying to fill QA histograms but StPicoD0AnaHists::mFillQaHists is false -- ignoring attemp! " << endm;
   }

   int EtaIndex = getEtaIndexRatio(Eta);
   int PhiIndex = getPhiIndexRatio(Phi);
   int VzIndex = getVzIndexRatio(Vz);
   if (IsPion)
   {
      mh2HFT1PtCentPartEtaVzPhi[0][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   if (IsKaon)
   {
      mh2HFT1PtCentPartEtaVzPhi[1][EtaIndex][VzIndex][PhiIndex]->Fill(pt, centrality);
   }
   mh2HFT1PtCent->Fill(pt, centrality);
   if (fabs(Eta) < 0.1 && pt > 3.0) mh2HFT1PhiVz->Fill(Phi, Vz);
}
//-----------------------------------------------------------------------
void StPicoD0AnaHists::addKaonPion(StKaonPion const* const kp, bool unlike, bool tpc, bool tof, int centrality, const double reweight)
{
   if (unlike)
   {
      if (tpc) mh2InvariantMassVsPt->Fill(kp->pt(), kp->m(), reweight);
      if (tof) mh2InvariantMassVsPtTof->Fill(kp->pt(), kp->m(), reweight);
      if (tpc) mh3InvariantMassVsPtVsCent->Fill(kp->pt(), centrality, kp->m(), reweight);
      if (tof) mh3InvariantMassVsPtVsCentTof->Fill(kp->pt(), centrality, kp->m(), reweight);
   }
   else
   {
      if (tpc) mh2InvariantMassVsPtLike->Fill(kp->pt(), kp->m(), reweight);
      if (tof) mh2InvariantMassVsPtTofLike->Fill(kp->pt(), kp->m(), reweight);
      if (tpc) mh3InvariantMassVsPtVsCentLike->Fill(kp->pt(), centrality, kp->m(), reweight);
      if (tof) mh3InvariantMassVsPtVsCentTofLike->Fill(kp->pt(), centrality, kp->m(), reweight);
   }
}
//---------------------------------------------------------------------
void StPicoD0AnaHists::addDcaPtCent(float dca, float dcaXy, float dcaZ, bool IsPion, bool IsKaon, float pt,  int centrality, float Eta, float Phi, float Vz, float ZdcX)
{
   if (!mFillQaHists)
   {
      LOG_ERROR << " You are trying to fill QA histograms but StPicoD0AnaHists::mFillQaHists is false -- ignoring attemp! " << endm;
   }

   int EtaIndex = getEtaIndexDca(Eta);
   int PhiIndex = getPhiIndexDca(Phi);
   int VzIndex = getVzIndexDca(Vz);

   if (centrality < 7) return; // Only analysis 0-10% dor DCA, memory Issue
   centrality = -1;

   if (IsPion)
   {
      mh3DcaXyZPtCentPartEtaVzPhi[0][EtaIndex][VzIndex][centrality + 1][PhiIndex]->Fill(pt, dcaXy, dcaZ);
   }
   if (IsKaon)
   {
      mh3DcaXyZPtCentPartEtaVzPhi[1][EtaIndex][VzIndex][centrality + 1][PhiIndex]->Fill(pt, dcaXy, dcaZ);
   }
   mh3DcaPtCent->Fill(pt, centrality, dca);
   mh3DcaXyPtCent->Fill(pt, centrality, dcaXy);
   mh3DcaZPtCent->Fill(pt, centrality, dcaZ);
}
//---------------------------------------------------------------------
int StPicoD0AnaHists::getEtaIndexDca(float Eta)
{
   for (int i = 0; i < anaCuts::nEtasDca; i++)
   {
      if ((Eta >= anaCuts::EtaEdgeDca[i]) && (Eta < anaCuts::EtaEdgeDca[i + 1]))
         return i;
   }
   return anaCuts::nEtasDca - 1;
}
//---------------------------------------------------------------------
int StPicoD0AnaHists::getPhiIndexDca(float Phi)
{
   for (int i = 0; i < anaCuts::nPhisDca; i++)
   {
      if ((Phi >= anaCuts::PhiEdgeDca[i]) && (Phi < anaCuts::PhiEdgeDca[i + 1]))
         return i;
   }
   return anaCuts::nPhisDca - 1;
}
//---------------------------------------------------------------------
int StPicoD0AnaHists::getVzIndexDca(float Vz)
{
   for (int i = 0; i < anaCuts::nVzsDca; i++)
   {
      if ((Vz >= anaCuts::VzEdgeDca[i]) && (Vz < anaCuts::VzEdgeDca[i + 1]))
         return i;
   }
   return anaCuts::nVzsDca - 1;
}
//---------------------------------------------------------------------
// int StPicoD0AnaHists::getZdcxIndex(float ZdcX)
// {
//    for (int i = 0; i < anaCuts::nZdcxs; i++)
//    {
//       if ((ZdcX >= anaCuts::ZdcxEdge[i]) && (ZdcX < anaCuts::ZdcxEdge[i + 1]))
//          return i;
//    }
//    return anaCuts::nZdcxs - 1;
// }
//---------------------------------------------------------------------
int StPicoD0AnaHists::getEtaIndexRatio(float Eta)
{
   for (int i = 0; i < anaCuts::nEtasRatio; i++)
   {
      if ((Eta >= anaCuts::EtaEdgeRatio[i]) && (Eta < anaCuts::EtaEdgeRatio[i + 1]))
         return i;
   }
   return anaCuts::nEtasRatio - 1;
}
//---------------------------------------------------------------------
int StPicoD0AnaHists::getPhiIndexRatio(float Phi)
{
   for (int i = 0; i < anaCuts::nPhisRatio; i++)
   {
      if ((Phi >= anaCuts::PhiEdgeRatio[i]) && (Phi < anaCuts::PhiEdgeRatio[i + 1]))
         return i;
   }
   return anaCuts::nPhisRatio - 1;
}
//---------------------------------------------------------------------
int StPicoD0AnaHists::getVzIndexRatio(float Vz)
{
   for (int i = 0; i < anaCuts::nVzsRatio; i++)
   {
      if ((Vz >= anaCuts::VzEdgeRatio[i]) && (Vz < anaCuts::VzEdgeRatio[i + 1]))
         return i;
   }
   return anaCuts::nVzsRatio - 1;
}
//---------------------------------------------------------------------
void StPicoD0AnaHists::addQaNtuple(int runnumber, float dca, float vz, float pt, float eta, float phi, int centrality, const double refmultCor, float zdcx, int tofMatchFlag, int hftMatchFlag)
{
//  nt->Fill(runnumber, dca, vz, pt, eta, phi, centrality, refmultCor, zdcx, tofMatchFlag, hftMatchFlag);
}
//---------------------------------------------------------------------
void StPicoD0AnaHists::closeFile()
{
   mOutFile->cd();

   mh1TotalEventsInRun->Write();
   mh1TotalEventsInRunBeforeCut->Write();
   mh2InvariantMassVsPt->Write();
   mh2InvariantMassVsPtLike->Write();
   mh2InvariantMassVsPtTof->Write();
   mh2InvariantMassVsPtTofLike->Write();
   //centrality
   mh1Cent->Write();
   mh1CentWg->Write();
   mh1gRefmultCor->Write();
   mh1gRefmultCorWg->Write();
   mh2CentVz->Write();
   mh2CentVzWg->Write();
   mh3InvariantMassVsPtVsCent->Write();
   mh3InvariantMassVsPtVsCentLike->Write();
   mh3InvariantMassVsPtVsCentTof->Write();
   mh3InvariantMassVsPtVsCentTofLike->Write();

   if (!mFillQaHists)
   {
      mOutFile->Close();
      mOutFile->Delete();
      return;
   }

   //HFT ratio QA
   mh2Tpc1PtCent->Write();
   mh2Tpc1PhiVz->Write();
   mh2HFT1PhiVz->Write();
   mh2HFT1PtCent->Write();

   //HFT DCA Ratio
   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < anaCuts::nEtasDca; iEta++)
      {
         for (int iVz = 0; iVz < anaCuts::nVzsDca; iVz++)
         {
            for (int iCent = 0; iCent < anaCuts::nCentsDca; iCent++)
            {
               for (int iPhi = 0; iPhi < anaCuts::nPhisDca; iPhi++)
               {
                  mh3DcaXyZPtCentPartEtaVzPhi[iParticle][iEta][iVz][iCent][iPhi]->Write();
               }
            }
         }
      }
   }

   for (int iParticle = 0; iParticle < anaCuts::nParticles; iParticle++)
   {
      for (int iEta = 0; iEta < anaCuts::nEtasRatio; iEta++)
      {
         for (int iVz = 0; iVz < anaCuts::nVzsRatio; iVz++)
         {
            for (int iPhi = 0; iPhi < anaCuts::nPhisRatio; iPhi++)
            {
               mh2Tpc1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]->Write();
               mh2HFT1PtCentPartEtaVzPhi[iParticle][iEta][iVz][iPhi]->Write();
            }
         }
      }
   }

   mh3DcaPtCent->Write();
   mh3DcaXyPtCent->Write();
   mh3DcaZPtCent->Write();

   // nt->Write();

   mOutFile->Close();
   mOutFile->Delete();
}
