#include <string>

using std::string;

void runPicoKPiXAnaMaker(string kPixList, string outFileName, string badRunListFileName = "picoList_bad_MB.list")
{
  //Check STAR Library. Please set SL_version to the original star library used in the production from http://www.star.bnl.gov/devcgi/dbProdOptionRetrv.pl
  string SL_version = "SL16d";
  string env_SL = getenv("STAR");
  if (env_SL.find(SL_version) == string::npos)
  {
    cout << "Environment Star Library does not match the requested library in runPicoD0EventMaker.C. Exiting..." << endl;
    exit(1);
  }

  gROOT->LoadMacro("$STAR/StRoot/StMuDSTMaker/COMMON/macros/loadSharedLibraries.C");
  loadSharedLibraries();

  gSystem->Load("StPicoDstMaker");
  gSystem->Load("StPicoPrescales");
  gSystem->Load("StPicoCharmContainers");
  gSystem->Load("StPicoCharmHists");
  gSystem->Load("StPicoKPiXAnaMaker");
  gSystem->Load("StRefMultCorr");

  chain = new StChain();

  // create list of picoDst files
  TString command = "sed 's/hft/picodsts/g' " + kPixList + " >correspondingPico.list";
  gSystem->Exec(command.Data());
  command = "sed -i 's/picoKPiX/picoDst/g' correspondingPico.list";
  gSystem->Exec(command.Data());
  command = "sed -i 's/Pico16a/physics2/g' correspondingPico.list";
  gSystem->Exec(command.Data());
  command = "sed -i 's/KPiX//g' correspondingPico.list";
  gSystem->Exec(command.Data());

  StPicoDstMaker*     picoDstMaker = new StPicoDstMaker(0, "correspondingPico.list", "picoDstMaker");
  StRefMultCorr*      grefmultCorrUtil  = CentralityMaker::instance()->getgRefMultCorr();
  StPicoKPiXAnaMaker* picoKPiXAnaMaker = new StPicoKPiXAnaMaker("picoD0AnaMaker", kPixList, outFileName, picoDstMaker, grefmultCorrUtil);
  picoKPiXAnaMaker->fillDpmHists(true);
  picoKPiXAnaMaker->fillDsHists(true);
  picoKPiXAnaMaker->fillLcHists(false);

  grefmultCorrUtil->setVzForWeight(6, -6.0, 6.0);
  grefmultCorrUtil->readScaleForWeight("StRoot/StRefMultCorr/macros/weight_grefmult_vpd30_vpd5_Run14.txt");

  for(Int_t i=0;i<6;i++)
  {
    cout << i << " " << grefmultCorrUtil->get(i, 0) << endl;
  }

  // -------------- USER variables -------------------------

  chain->Init();

  int nEntries = picoKPiXAnaMaker->getEntries();
  cout<<"Processing "<<nEntries<<" events..."<<endl;
  for (int iEvent = 0; iEvent < nEntries; ++iEvent)
  {
    chain->Clear();
    if(iEvent && iEvent%2000 == 0) cout<<"... finished processing "<<iEvent<<" events."<<endl;

    int iret = chain->Make();
    if (iret)
    {
      cout << "Bad return code!" << iret << endl;
      break;
    }
  }
  cout<<"Finished processing "<<nEntries<<" events."<<endl;

  chain->Finish();
  delete chain;

  // delete list of picos
  command = "rm -f correspondingPico.list";
  gSystem->Exec(command.Data());
}
