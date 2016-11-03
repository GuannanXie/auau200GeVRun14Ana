#ifndef StPicoD0AnaMaker_h
#define StPicoD0AnaMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoD0Event
 *  simultaneously and do analysis.
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include "TChain.h"
#include "StMaker.h"
#include "StThreeVectorF.hh"

class TString;
class TFile;
class TNtuple;
class StPicoEvent;
class StPicoD0Event;
class StKaonPion;
class StPicoTrack;
class StPicoDstMaker;
class StPicoD0AnaHists;
class StRefMultCorr;

class StPicoD0AnaMaker : public StMaker
{
public:
   StPicoD0AnaMaker(char const * name, TString const inputFilesList,
                    TString const outBaseName, StPicoDstMaker* picoDstMaker, StRefMultCorr* grefmultCorrUtil);
   virtual ~StPicoD0AnaMaker();

   virtual Int_t Init();
   virtual Int_t Make();
   virtual Int_t Finish();

    int getEntries() const;
    void fillQaHistograms(bool b = true);
    void fillBackgroundTrees(bool b = true);

private:
   StPicoD0AnaMaker() {}
   void readNextEvent();

    int  getD0PtIndex(StKaonPion const* ) const;
    bool isGoodTrigger(StPicoEvent const*) const;
    bool isGoodEvent(StPicoEvent const*, StThreeVectorF const& vtx) const;
    bool isGoodQaTrack(StPicoTrack const* ,StThreeVectorF const& momentum ,double dca) const;
    bool isGoodTrack(StPicoTrack const*, StThreeVectorF const&) const;
    bool isTpcPion(StPicoTrack const*) const;
    bool isTpcKaon(StPicoTrack const*) const;
    bool isTpcProton(StPicoTrack const*) const;
    bool isTofPion(StPicoTrack const*, float beta, StThreeVectorF const& vtx) const;
    bool isTofKaon(StPicoTrack const*, float beta, StThreeVectorF const& vtx) const;
    bool isTofProton(StPicoTrack const*, float beta, StThreeVectorF const& vtx) const;
    bool isGoodPair(StKaonPion const*) const;
    bool isSideBand(float m) const;
    float getTofBeta(StPicoTrack const*,StThreeVectorF const& vtx) const;
    int trkHalf(StPicoTrack const*, StThreeVectorF const& vtx) const;

   StPicoDstMaker* mPicoDstMaker;
   StPicoD0Event* mPicoD0Event;
   StRefMultCorr* mGRefMultCorrUtil;

    TString mInputFilesList;
    TString mOutFileBaseName;
    TChain* mChain;
    int mEventCounter;
    bool mFillQaHists;
    bool mFillBackgroundTrees;

   // -------------- USER variables -------------------------
   // add your member variables here.
   // Remember that ntuples size can be really big, use histograms where appropriate
   StPicoD0AnaHists* mHists;

   ClassDef(StPicoD0AnaMaker, 1)
};

inline void StPicoD0AnaMaker::fillQaHistograms(bool b) { mFillQaHists = b;}
inline void StPicoD0AnaMaker::fillBackgroundTrees(bool b) { mFillBackgroundTrees = b;}

inline int StPicoD0AnaMaker::getEntries() const
{
   return mChain ? mChain->GetEntries() : 0;
}

inline void StPicoD0AnaMaker::readNextEvent()
{
   mChain->GetEntry(mEventCounter++);
}
#endif
