/* **************************************************
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "TFile.h"
#include "TMath.h"
#include "TClonesArray.h"
#include "TTree.h"
#include "TNtuple.h"

#include "StEvent/StDcaGeometry.h"
#include "StPhysicalHelixD.hh"
#include "phys_constants.h"
#include "StRoot/StPicoDstMaker/StPicoBTofPidTraits.h"
#include "StBTofUtil/tofPathLength.hh"
#include "StPicoDstMaker/StPicoDstMaker.h"
#include "StPicoDstMaker/StPicoDst.h"
#include "StPicoDstMaker/StPicoEvent.h"
#include "StPicoDstMaker/StPicoTrack.h"
#include "StPicoCharmContainers/StPicoD0Event.h"
#include "StPicoCharmContainers/StKaonPion.h"
#include "StPicoCharmContainers/StD0SoftPion.h"
#include "StPicoDstarAnaMaker.h"
#include "StPicoDstarAnaHists.h"
#include "StAnaCuts.h"
#include "StRoot/StRefMultCorr/StRefMultCorr.h"

#include "StMemStat.h"

ClassImp(StPicoDstarAnaMaker)

StPicoDstarAnaMaker::StPicoDstarAnaMaker(char const * name, TString const inputFilesList,
                                   TString const outFileBaseName, StPicoDstMaker* picoDstMaker, StRefMultCorr* grefmultCorrUtil):
   StMaker(name), mPicoDstMaker(picoDstMaker), mPicoD0Event(NULL), mGRefMultCorrUtil(grefmultCorrUtil),
   mInputFilesList(inputFilesList), mOutFileBaseName(outFileBaseName),
   mChain(NULL), mEventCounter(0), mFillQaHists(false), mFillBackgroundTrees(false), mHists(NULL)
{}

Int_t StPicoDstarAnaMaker::Init()
{
   mPicoD0Event = new StPicoD0Event();

   mChain = new TChain("T");
   std::ifstream listOfFiles(mInputFilesList.Data());
   if (listOfFiles.is_open())
   {
      std::string file;
      while (getline(listOfFiles, file))
      {
         LOG_INFO << "StPicoDstarAnaMaker - Adding :" << file << endm;
         mChain->Add(file.c_str());
      }
   }
   else
   {
      LOG_ERROR << "StPicoDstarAnaMaker - Could not open list of files. ABORT!" << endm;
      return kStErr;
   }

   mChain->GetBranch("dEvent")->SetAutoDelete(kFALSE);
   mChain->SetBranchAddress("dEvent", &mPicoD0Event);

   mOutFileBaseName = mOutFileBaseName.ReplaceAll(".root", "");
   // -------------- USER VARIABLES -------------------------
   mHists = new StPicoDstarAnaHists(mOutFileBaseName, mFillQaHists, mFillBackgroundTrees);

   return kStOK;
}
//-----------------------------------------------------------------------------
StPicoDstarAnaMaker::~StPicoDstarAnaMaker()
{
}
//-----------------------------------------------------------------------------
Int_t StPicoDstarAnaMaker::Finish()
{
   mHists->closeFile();
   return kStOK;
}
//-----------------------------------------------------------------------------
Int_t StPicoDstarAnaMaker::Make()
{
   StMemStat mem;
   readNextEvent();

   if (!mPicoDstMaker)
   {
      LOG_WARN << " StPicoDstarAnaMaker - No PicoDstMaker! Skip! " << endm;
      return kStWarn;
   }

   StPicoDst const* picoDst = mPicoDstMaker->picoDst();

   if (!picoDst)
   {
      LOG_WARN << "StPicoDstarAnaMaker - No PicoDst! Skip! " << endm;
      return kStWarn;
   }

   if (mPicoD0Event->runId() != picoDst->event()->runId() ||
         mPicoD0Event->eventId() != picoDst->event()->eventId())
   {
      LOG_ERROR << " StPicoDstarAnaMaker - !!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!" << "\n";
      LOG_ERROR << " StPicoDstarAnaMaker - SOMETHING TERRIBLE JUST HAPPENED. StPicoEvent and StPicoD0Event are not in sync." << endm;
      exit(1);
   }

   // -------------- USER ANALYSIS -------------------------

   mGRefMultCorrUtil->init(picoDst->event()->runId());

   if (!mGRefMultCorrUtil)
   {
      LOG_WARN << " No mGRefMultCorrUtil! Skip! " << endl;
      return kStWarn;
   }

   if (mGRefMultCorrUtil->isBadRun(picoDst->event()->runId()))
   {
      //cout<<"This is a bad run from mGRefMultCorrUtil! Skip! " << endl;
      return kStOK;
   }

   if (isGoodTrigger(picoDst->event()))
   {
      mHists->addEventBeforeCut(picoDst->event());

      StThreeVectorF pVtx = picoDst->event()->primaryVertex();
      if (isGoodEvent(picoDst->event(), pVtx))
      {
         TClonesArray const* aKaonPion = mPicoD0Event->kaonPionArray();
         mHists->addEvent(picoDst->event());

         mGRefMultCorrUtil->initEvent(picoDst->event()->grefMult(), pVtx.z(), picoDst->event()->ZDCx()) ;

         int centrality  = mGRefMultCorrUtil->getCentralityBin9();
         const double reweight = mGRefMultCorrUtil->getWeight();
         const double refmultCor = mGRefMultCorrUtil->getRefMultCorr();
         mHists->addCent(refmultCor, centrality, reweight, pVtx.z());

         //Basiclly add some QA plots
         UInt_t nTracks = picoDst->numberOfTracks();

         if (mFillQaHists)
         {
            for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
            {
               StPicoTrack const* trk = picoDst->track(iTrack);
               if (!trk) continue;
               StPhysicalHelixD helix = trk->helix();
               float dca = float(helix.geometricSignedDistance(pVtx));
               StThreeVectorF momentum = trk->gMom(pVtx, picoDst->event()->bField());

               if (!isGoodQaTrack(trk, momentum, dca)) continue;

               StThreeVectorF dcaPoint = helix.at(helix.pathLength(pVtx.x(), pVtx.y()));
               float dcaZ = dcaPoint.z() - pVtx.z();
               double dcaXy = helix.geometricSignedDistance(pVtx.x(), pVtx.y());

               bool tpcPion = isTpcPion(trk);
               bool tpcKaon = isTpcKaon(trk);
               bool tpcProton = isTpcProton(trk);
               float piBeta = getTofBeta(trk, pVtx);
               float kBeta = piBeta;
               float pBeta = piBeta;
               bool piTofAvailable = !isnan(piBeta) && piBeta > 0;
               bool kTofAvailable = !isnan(kBeta) && kBeta > 0;
               bool pTofAvailable = !isnan(pBeta) && pBeta > 0;
               bool tofPion = isTofPion(trk, piBeta, pVtx);
               bool tofKaon = isTofKaon(trk, kBeta, pVtx);
               bool tofProton = isTofProton(trk, pBeta, pVtx);

               bool goodPion = (piTofAvailable && tofPion && tpcPion) || (!piTofAvailable && tpcPion);//Always require TPC
               bool goodKaon = (kTofAvailable && tofKaon && tpcKaon) || (!kTofAvailable && tpcKaon);
               bool goodProton = (pTofAvailable && tofProton && tpcProton) || (!pTofAvailable && tpcProton);
               // bool goodKaon = (momentum.perp() <= 1.6 && kTofAvailable && tofKaon && tpcKaon) || (momentum.perp() > 1.6 && tpcKaon);//strict Kaon pid

               // cout<<"mem.Used() = " << mem.Used()<<endl;
               // cout<<"mem.Free() = " << mem.Free()<<endl;
               // cout<<"mem.ProgSize() = " << mem.ProgSize()<<endl;

               // bool goodPion = piTofAvailable && tofPion && tpcPion;//strict Pion pid
               // bool goodKaon = kTofAvailable && tofKaon && tpcKaon;//strict Kaon pid
               // bool goodProton = pTofAvailable && tofProton && tpcProton;//strict Proton pid

               if (trk  && fabs(dca) < 1.5 && trk->isHFTTrack() && (goodPion || goodKaon || goodProton))
               {
                  mHists->addDcaPtCent(dca, dcaXy, dcaZ, goodPion, goodKaon, goodProton, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), pVtx.z(), picoDst->event()->ZDCx() / 1000.); //add Dca distribution
               }
               if (trk  && fabs(dca) < 1.5 && (goodPion || goodKaon || goodProton))
               {
                  mHists->addTpcDenom1(goodPion, goodKaon, goodProton, momentum.perp(), centrality, momentum.pseudoRapidity(), momentum.phi(), pVtx.z(), picoDst->event()->ZDCx() / 1000.); //Dca cut on 1.5cm, add Tpc Denominator
               }
               if (trk && fabs(dca) < 1.5 && trk->isHFTTrack() && (goodPion || goodKaon || goodProton) && fabs(dcaXy) < 1. && fabs(dcaZ) < 1.)
               {
                  mHists->addHFTNumer1(goodPion, goodKaon, goodProton, momentum.perp(), centrality,  momentum.pseudoRapidity(), momentum.phi(), pVtx.z(), picoDst->event()->ZDCx() / 1000.); //Dca cut on 1.5cm, add HFT Numerator
               }
            } // .. end tracks loop
         }

         for (int idx = 0; idx < aKaonPion->GetEntries(); ++idx)
         {
            StKaonPion const* kp = (StKaonPion*)aKaonPion->UncheckedAt(idx);

            bool goodPair = isGoodPair(kp);
            if (!goodPair && !mFillBackgroundTrees) continue;

            StPicoTrack const* kaon = picoDst->track(kp->kaonIdx());
            StPicoTrack const* pion = picoDst->track(kp->pionIdx());

            if (!isGoodTrack(kaon, pVtx) || !isGoodTrack(pion, pVtx)) continue;

            // PID
            if (!isTpcPion(pion) || !isTpcKaon(kaon)) continue;
            float piBeta = getTofBeta(pion, pVtx);
            float kBeta = getTofBeta(kaon, pVtx);
            bool piTofAvailable = !isnan(piBeta) && piBeta > 0;
            bool kTofAvailable = !isnan(kBeta) && kBeta > 0;
            bool tofPion = piTofAvailable ? isTofPion(pion, piBeta, pVtx) : true;//this is bybrid pid, not always require tof
            bool tofKaon = kTofAvailable ? isTofKaon(kaon, kBeta, pVtx) : true;//this is bybrid pid, not always require tof
            bool tof = tofPion && tofKaon;
            if(!tof) continue;

            bool unlikeD0 = kaon->charge() * pion->charge() < 0 ? true : false;

            if (!goodPair) continue;
            if (!unlikeD0) continue;
            if ( (kp->m()<1.82) || (kp->m()>1.92) ) continue;

            for (unsigned short iTrack = 0; iTrack < nTracks; ++iTrack)
            {
               StPicoTrack const* sptrk = picoDst->track(iTrack);
               if (!sptrk) continue;
               StPhysicalHelixD sphelix = sptrk->helix();
               float spdca = float(sphelix.geometricSignedDistance(pVtx));
               StThreeVectorF spmomentum = sptrk->gMom(pVtx, picoDst->event()->bField());//global momentum
               // StThreeVectorF spmomentum = sptrk->pMom();//primary momentum

               if (!isGoodSoftPionTrack(sptrk, spmomentum, spdca)) continue;
               if (fabs(spdca)>1.0) continue;
               
               if (!isTpcPion(sptrk)) continue;
               float spBeta = getTofBeta(sptrk, pVtx);
               bool spTofAvailable = !isnan(spBeta) && spBeta > 0;
               bool tofSoftPion = spTofAvailable ? isTofPion(sptrk, spBeta, pVtx) : true;//this is bybrid pid, not always require tof
               if(!tofSoftPion) continue;

               if((kp->kaonIdx() == iTrack) || (kp->pionIdx() == iTrack)) continue;

               StD0SoftPion const D0SoftPion(kp, sptrk, kp->kaonIdx(), kp->pionIdx(), iTrack, pVtx, picoDst->event()->bField());

               bool unlikeDstar = pion->charge() * sptrk->charge() > 0 ? true : false;

               if( (kp->pt()/spmomentum.perp()) < 10 || (kp->pt()/spmomentum.perp()) > 20) continue;

               mHists->addD0SoftPion(D0SoftPion, kp, unlikeDstar, true, true, centrality, reweight);
            }

            // if (goodPair) mHists->addKaonPion(kp, unlikeD0, true, tof, centrality, reweight);

            // if (mFillBackgroundTrees && tof)
            // {
            //    if (centrality < 0) continue;
            //    int const ptBin = getD0PtIndex(kp);
            //
            //    if (unlike)
            //    {
            //       if (kp->m() > anaCuts::likeSignMassRange.first && kp->m() < anaCuts::likeSignMassRange.second) mHists->addBackground(kp, kaon, pion, ptBin, false);
            //    }
            //    else if (isSideBand(kp->m()))
            //    {
            //       mHists->addBackground(kp, kaon, pion, ptBin, true);
            //    }
            // }

         } // end of kaonPion loop
      } // end of isGoodEvent
   } // end of isGoodTrigger

   return kStOK;
}
//-----------------------------------------------------------------------------
int StPicoDstarAnaMaker::getD0PtIndex(StKaonPion const* const kp) const
{
   for (int i = 0; i < anaCuts::nPtBins; i++)
   {
      if ((kp->pt() >= anaCuts::PtBinsEdge[i]) && (kp->pt() < anaCuts::PtBinsEdge[i + 1]))
         return i;
   }
   return anaCuts::nPtBins - 1;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isGoodEvent(StPicoEvent const* const picoEvent, StThreeVectorF const& pVtx) const
{
   return fabs(pVtx.z()) < anaCuts::vz &&
          fabs(pVtx.z() - picoEvent->vzVpd()) < anaCuts::vzVpdVz &&
          !(fabs(pVtx.x()) < anaCuts::Verror && fabs(pVtx.y()) < anaCuts::Verror && fabs(pVtx.z()) < anaCuts::Verror) &&
          sqrt(TMath::Power(pVtx.x(), 2) + TMath::Power(pVtx.y(), 2)) <=  anaCuts::Vrcut;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isGoodTrigger(StPicoEvent const* const picoEvent) const
{
   for (auto trg : anaCuts::triggers)
   {
      if (picoEvent->isTrigger(trg)) return true;
   }

   return false;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isGoodSoftPionTrack(StPicoTrack const* const trk, StThreeVectorF const& momentum, double const dca) const
{
   // return trk->gPt() > anaCuts::minSoftPt&& 
   return momentum.perp() > anaCuts::minSoftPt&& 
     trk->nHitsFit() >= anaCuts::nHitsFit&& 
     fabs(momentum.pseudoRapidity()) <= anaCuts::Eta;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isGoodQaTrack(StPicoTrack const* const trk, StThreeVectorF const& momentum, double const dca) const
{
   return trk->gPt() > anaCuts::qaGPt && trk->nHitsFit() >= anaCuts::qaNHitsFit && fabs(momentum.pseudoRapidity()) <= anaCuts::Eta;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isGoodTrack(StPicoTrack const* const trk, StThreeVectorF const& vtx) const
{
   StThreeVectorF mom = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField());


   return mom.perp() > anaCuts::minPt &&
          trk->nHitsFit() >= anaCuts::nHitsFit &&
          fabs(mom.pseudoRapidity()) <= anaCuts::Eta;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isTpcPion(StPicoTrack const* const trk) const
{
   return fabs(trk->nSigmaPion()) < anaCuts::nSigmaPion;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isTpcKaon(StPicoTrack const* const trk) const
{
   return fabs(trk->nSigmaKaon()) < anaCuts::nSigmaKaon;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isTpcProton(StPicoTrack const* const trk) const
{
   return fabs(trk->nSigmaProton()) < anaCuts::nSigmaProton;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isGoodPair(StKaonPion const* const kp) const
{
   int tmpIndex = getD0PtIndex(kp);
   return cos(kp->pointingAngle()) > anaCuts::cosTheta[tmpIndex] &&
          // kp->cosThetaStar() < 0.8 && 
          kp->pionDca() > anaCuts::pDca[tmpIndex] && kp->kaonDca() > anaCuts::kDca[tmpIndex] &&
          kp->dcaDaughters() < anaCuts::dcaDaughters[tmpIndex] &&
          kp->decayLength() > anaCuts::decayLength[tmpIndex] &&
          fabs(kp->lorentzVector().rapidity()) < anaCuts::RapidityCut &&
          ((kp->decayLength()) * sin(kp->pointingAngle())) < anaCuts::dcaV0ToPv[tmpIndex];
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isSideBand(float const m) const
{
   if (m > anaCuts::sideBandMassRange0.first && m < anaCuts::sideBandMassRange0.second) return true;
   if (m > anaCuts::sideBandMassRange1.first && m < anaCuts::sideBandMassRange1.second) return true;

   return false;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isTofKaon(StPicoTrack const* const trk, float beta, StThreeVectorF const& vtx) const
{
   bool tofKaon = false;

   if (beta > 0)
   {
      double ptot = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField()).mag();
      float beta_k = ptot / sqrt(ptot * ptot + M_KAON_PLUS * M_KAON_PLUS);
      tofKaon = fabs(1 / beta - 1 / beta_k) < anaCuts::kTofBetaDiff ? true : false;
   }

   return tofKaon;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isTofPion(StPicoTrack const* const trk, float beta, StThreeVectorF const& vtx) const
{
   bool tofPion = false;

   if (beta > 0)
   {
      double ptot = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField()).mag();
      float beta_pi = ptot / sqrt(ptot * ptot + M_PION_PLUS * M_PION_PLUS);
      tofPion = fabs(1 / beta - 1 / beta_pi) < anaCuts::piTofBetaDiff ? true : false;
   }

   return tofPion;
}
//-----------------------------------------------------------------------------
bool StPicoDstarAnaMaker::isTofProton(StPicoTrack const* const trk, float beta, StThreeVectorF const& vtx) const
{
   bool tofProton = false;

   if (beta > 0)
   {
      double ptot = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField()).mag();
      float beta_p = ptot / sqrt(ptot * ptot + M_PROTON * M_PROTON);
      tofProton = fabs(1 / beta - 1 / beta_p) < anaCuts::pTofBetaDiff ? true : false;
   }

   return tofProton;
}
//-----------------------------------------------------------------------------
float StPicoDstarAnaMaker::getTofBeta(StPicoTrack const* const trk, StThreeVectorF const& vtx) const
{
   int index2tof = trk->bTofPidTraitsIndex();

   float beta = std::numeric_limits<float>::quiet_NaN();

   if (index2tof >= 0)
   {
      StPicoBTofPidTraits const* const tofPid = mPicoDstMaker->picoDst()->btofPidTraits(index2tof);

      if (tofPid)
      {
         beta = tofPid->btofBeta();

         if (beta < 1e-4)
         {
            StThreeVectorF const btofHitPos = tofPid->btofHitPos();

            StPhysicalHelixD helix = trk->helix();
            float L = tofPathLength(&vtx, &btofHitPos, helix.curvature());
            float tof = tofPid->btof();
            if (tof > 0) beta = L / (tof * (C_C_LIGHT / 1.e9));
            else beta = std::numeric_limits<float>::quiet_NaN();
         }
      }
   }

   return beta;
}
//-----------------------------------------------------------------------------
int StPicoDstarAnaMaker::trkHalf(StPicoTrack const* const trk, StThreeVectorF const& vtx) const
{
   StThreeVectorF mom = trk->gMom(vtx, mPicoDstMaker->picoDst()->event()->bField());
   if (mom.phi() > anaCuts:: rightHalfLowEdge && mom.phi() < anaCuts:: rightHalfHighEdge) return +1; //right side
   else return -1;//lest side

}
