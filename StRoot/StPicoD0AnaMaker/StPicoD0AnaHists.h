#ifndef StPicoD0AnaHists__h
#define StPicoD0AnaHists__h

/* **************************************************
 *  A Maker to read a StPicoEvent and StPicoD0Event
 *  simultaneously and do analysis.
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include "StAnaCuts.h"
#include "THnSparse.h"
#include "THn.h"

class TH1F;
class TH2F;
class TH3F;
class TFile;
class TString;
class StPicoPrescales;
class StPicoEvent;
class StKaonPion;
class TNtuple;


class StPicoD0AnaHists
{
public:
   StPicoD0AnaHists(TString fileBaseName, bool fillQaHists = true);
   virtual ~StPicoD0AnaHists();
   void addEvent(StPicoEvent const *);
   void addEventBeforeCut(StPicoEvent const *);
   void addCent(const double refmultCor, int centrality, const double reweight, const float vz);
   void addKaonPion(StKaonPion const*, bool unlike, bool tpc, bool tof, int centrality, const double reweight);
   void addTpcDenom1(bool IsPion, bool IsKaon, float pt, int centrality, float Eta, float Phi, float Vz, float ZdcX);
   void addHFTNumer1(bool IsPion, bool IsKaon, float pt, int centrality, float Eta, float Phi, float Vz, float ZdcX);
   void addQaNtuple(int, float, float, float, float, float, int, const double, float, int, int);
   void addDcaPtCent(float dca, float dcaXy, float  dcaZ, bool IsPion, bool IsKaon, float pt,  int centrality, float Eta, float Phi, float Vz, float ZdcX);
   int getEtaIndexDca(float Eta) ;
   int getPhiIndexDca(float Phi) ;
   int getVzIndexDca(float Vz) ;
   // int getZdcxIndex(float ZdcX) ;
   int getEtaIndexRatio(float Eta) ;
   int getPhiIndexRatio(float Phi) ;
   int getVzIndexRatio(float Vz) ;
   void closeFile();

private:
   StPicoD0AnaHists() {}

   bool mFillQaHists;
   StPicoPrescales* mPrescales;
   TFile* mOutFile;
   TH1F* mh1TotalEventsInRun;
   TH1F* mh1TotalEventsInRunBeforeCut;
   TH2F* mh2InvariantMassVsPt;
   TH2F* mh2InvariantMassVsPtLike;
   TH2F* mh2InvariantMassVsPtTof;
   TH2F* mh2InvariantMassVsPtTofLike;
   //centrality
   TH1F* mh1Cent;
   TH1F* mh1CentWg;
   TH1F* mh1gRefmultCor;
   TH1F* mh1gRefmultCorWg;
   TH2F* mh2CentVz;
   TH2F* mh2CentVzWg;
   TH3F* mh3InvariantMassVsPtVsCent;
   TH3F* mh3InvariantMassVsPtVsCentLike;
   TH3F* mh3InvariantMassVsPtVsCentTof;
   TH3F* mh3InvariantMassVsPtVsCentTofLike;
   //d0
   TH3F* mh3d0InvariantMassVsPtVsCent;
   TH3F* mh3d0InvariantMassVsPtVsCentTof;
   //d0 bar
   TH3F* mh3d0barInvariantMassVsPtVsCent;
   TH3F* mh3d0barInvariantMassVsPtVsCentTof;
   //LikeSign Pos+Neg
   TH3F* mh3InvariantMassVsPtVsCentLikePos;
   TH3F* mh3InvariantMassVsPtVsCentTofLikePos;
   //LikeSign Pos+Neg
   TH3F* mh3InvariantMassVsPtVsCentLikeNeg;
   TH3F* mh3InvariantMassVsPtVsCentTofLikeNeg;
   //HFT Left
   TH3F* mh3InvariantMassVsPtVsCentLeft;
   TH3F* mh3InvariantMassVsPtVsCentLikeLeft;
   TH3F* mh3InvariantMassVsPtVsCentTofLeft;
   TH3F* mh3InvariantMassVsPtVsCentTofLikeLeft;
   //HFT Right
   TH3F* mh3InvariantMassVsPtVsCentRight;
   TH3F* mh3InvariantMassVsPtVsCentLikeRight;
   TH3F* mh3InvariantMassVsPtVsCentTofRight;
   TH3F* mh3InvariantMassVsPtVsCentTofLikeRight;
   //HFT Mixed
   TH3F* mh3InvariantMassVsPtVsCentMixed;
   TH3F* mh3InvariantMassVsPtVsCentLikeMixed;
   TH3F* mh3InvariantMassVsPtVsCentTofMixed;
   TH3F* mh3InvariantMassVsPtVsCentTofLikeMixed;
   //HFT ratio QA
   TH2F* mh2Tpc1PtCent;
   TH2F* mh2Tpc1PhiVz;
   TH2F* mh2HFT1PtCent;
   TH2F* mh2HFT1PhiVz;
   TH2F* mh2Tpc1PtCentPartEtaVzPhi[anaCuts::nParticles][anaCuts::nEtasRatio][anaCuts::nVzsRatio][anaCuts::nPhisRatio];
   TH2F* mh2HFT1PtCentPartEtaVzPhi[anaCuts::nParticles][anaCuts::nEtasRatio][anaCuts::nVzsRatio][anaCuts::nPhisRatio];

   //HFT Dca
   TH3F* mh3DcaXyZPtCentPartEtaVzPhi[anaCuts::nParticles][anaCuts::nEtasDca][anaCuts::nVzsDca][anaCuts::nCentsDca];

   TH3F* mh3DcaPtCent;
   TH3F* mh3DcaXyPtCent;
   TH3F* mh3DcaZPtCent;

   TNtuple* nt;
};
#endif
