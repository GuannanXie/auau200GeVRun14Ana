#ifndef StPicoKPiXAnaMaker_h
#define StPicoKPiXAnaMaker_h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoKPiEvent
 *  simultaneously and do analysis.
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include <string>
#include "TChain.h"
#include <vector>
#include "StMaker.h"
#include "StThreeVectorF.hh"

class TString;
class TFile;
class TNtuple;
class StPicoDstMaker;
class StPicoEvent;
class StPicoTrack;
class StPicoKPiXEvent;
class StPicoKPiX;
class StRefMultCorr;
class StPicoCharmMassHists;

namespace kPiXAnaCuts
{
  class TopologicalCuts;
}

class StPicoKPiXAnaMaker : public StMaker
{
  public:
    StPicoKPiXAnaMaker(char const * name, TString const inputFilesList, std::string const outBaseName,
                       StPicoDstMaker* picoDstMaker, StRefMultCorr* grefmultCorrUtil);
    virtual ~StPicoKPiXAnaMaker();

    virtual Int_t Init();
    virtual Int_t Make();
    virtual Int_t Finish();

    int getEntries() const;
    void fillDpmHists(bool b=true);
    void fillDsHists(bool b=true);
    void fillLcHists(bool b=true);

  private:
    StPicoKPiXAnaMaker() {}
    void readNextEvent();

    bool isGoodTrigger(StPicoEvent const*) const;
    bool isGoodEvent(StPicoEvent const*, StThreeVectorF const& vtx) const;
    bool isGoodTrack(StPicoTrack const*, StThreeVectorF const&) const;
    bool isTpcPion(StPicoTrack const*) const;
    bool isTpcKaon(StPicoTrack const*) const;
    bool isTpcProton(StPicoTrack const*) const;
    bool isTofPion(StPicoTrack const*, float beta, StThreeVectorF const& vtx) const;
    bool isTofKaon(StPicoTrack const*, float beta, StThreeVectorF const& vtx) const;
    bool isTofProton(StPicoTrack const*, float beta, StThreeVectorF const& vtx) const;
    float getTofBeta(StPicoTrack const*,StThreeVectorF const& vtx) const;
    int getIndex(float const value, std::vector<float> const& edges) const;
    bool isGoodKPiX(StPicoKPiX const* kpx, kPiXAnaCuts::TopologicalCuts const& cuts) const;

    int mEventCounter;
    StPicoDstMaker*  mPicoDstMaker;
    StPicoKPiXEvent* mPicoKPiXEvent;
    StRefMultCorr*   mGRefMultCorrUtil;
    TChain*          mChain;

    std::string mInputFilesList;
    std::string mOutFileBaseName;

    bool mFillDpmHists;
    bool mFillDsHists;
    bool mFillLcHists;

    StPicoCharmMassHists* mDpmHists;
    StPicoCharmMassHists* mDsHists;
    StPicoCharmMassHists* mLcHists;

    ClassDef(StPicoKPiXAnaMaker, 0)
};

inline int StPicoKPiXAnaMaker::getEntries() const
{
  return mChain? mChain->GetEntries() : 0;
}

inline void StPicoKPiXAnaMaker::fillDpmHists(bool b) { mFillDpmHists = b; }
inline void StPicoKPiXAnaMaker::fillDsHists(bool b)  { mFillDsHists = b; }
inline void StPicoKPiXAnaMaker::fillLcHists(bool b)  { mFillLcHists = b; }
inline void StPicoKPiXAnaMaker::readNextEvent()
{
  mChain->GetEntry(mEventCounter++);
}
#endif
