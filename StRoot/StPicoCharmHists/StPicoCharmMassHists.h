#ifndef StPicoCharmMassHists__h
#define StPicoCharmMassHists__h

/* **************************************************
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include <string>
#include "StarClassLibrary/StLorentzVectorF.hh"

class TH1F;
class TH2F;
class TH3F;
class TFile;
class StPicoPrescales;
class StPicoEvent;
class StPicoCharmMassHists;
class TNtuple;


/* Don't #include this class in another class header file (particularly in a ROOT class)
 * it won't compile with ROOT 5.xx. Forward declare it instead. */

class StPicoCharmMassHists
{
  public:
   StPicoCharmMassHists() = delete;
   StPicoCharmMassHists(std::string fileBaseName, std::string prescalesDir);
   virtual ~StPicoCharmMassHists();

   void addEvent(StPicoEvent const&);
   void addEventBeforeCut(StPicoEvent const&);
   void addCent(double refmultCor,int centrality, double reweight, float vz);
   void addKPiX(StLorentzVectorF const&, bool fg, int centrality, double reweight);
   void closeFile();

  private:

   StPicoPrescales* mPrescales = nullptr;
   TFile* mOutFile = nullptr;
   TH1F*  mh1TotalEventsInRun = nullptr;
   TH1F*  mh1TotalEventsInRunBeforeCut = nullptr;

   //centrality
   TH1F* mh1Cent = nullptr;
   TH1F* mh1CentWg = nullptr;
   TH1F* mh1gRefmultCor = nullptr;
   TH1F* mh1gRefmultCorWg = nullptr;
   TH2F* mh2CentVz = nullptr;
   TH2F* mh2CentVzWg = nullptr;

   TH3F* mh3InvariantMassVsPtVsCent = nullptr;
   TH3F* mh3InvariantMassVsPtVsCentBkg = nullptr;
};
#endif
