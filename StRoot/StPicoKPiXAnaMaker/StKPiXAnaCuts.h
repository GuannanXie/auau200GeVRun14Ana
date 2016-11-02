#ifndef StAnaCuts_H
#define StAnaCuts_H

/* **************************************************
 *
 *  Authors: Guannan Xie <guannanxie@lbl.gov>
 *           Mustafa Mustafa <mmustafa@lbl.gov>
 *
 * **************************************************
 */

#include "Rtypes.h"
#include <string>
#include <array>
#include <vector>

#include "phys_constants.h"

namespace kPiXAnaCuts
{
  // path to lists of triggers prescales
  // lists are obtained from http://www.star.bnl.gov/protected/common/common2014/trigger2014/plots_au200gev/
  std::string const prescalesFilesDirectoryName = "./run14AuAu200GeVPrescales";

  //trigger
  std::array<unsigned int, 5> const triggers = { 450050,    // vpdmb-5-p-nobsmd-hlt
                                                 450060,    // vpdmb-5-p-nobsmd-hlt
                                                 450005,    // vpdmb-5-p-nobsmd
                                                 450015,    // vpdmb-5-p-nobsmd
                                                 450025};    // vpdmb-5-p-nobsmd

  // event
  float const vz = 6.0;// cm.
  float const vzVpdVz = 3.0; // 3 cm.
  float const Verror = 1.0e-5; //
  float const Vrcut = 2.0; //

  // all tracks
  float const minPt = 0.5;
  int   const nHitsFit = 20;
  float const Eta = 1.0;

  //pions
  float const nSigmaPion = 3.0;
  float const pTofBetaDiff = 0.03;

  //kaons
  float const nSigmaKaon = 2.0;
  float const kTofBetaDiff = 0.03;

  //protons
  float const nSigmaProton = 2.0;
  float const protonTofBetaDiff = 0.03;

  struct TopologicalCuts
  {
    double xMassHypothesis;
    float  minMass;
    float  maxMass;
    float  rapidityCut;
    int    nPtBins;
    std::vector<float> ptBinsEdge;
    std::vector<float> dcaV0ToPv;
    std::vector<float> decayLength;
    std::vector<float> cosTheta;
    std::vector<float> dcaDaughters;
    std::vector<float> kDca;
    std::vector<float> pDca;
    std::vector<float> xDca;

    TopologicalCuts(double xMassHypothesis,
                    float  minMass,
                    float  maxMass,
                    float  rapidityCut,
                    int    nPtBins,
                    std::vector<float> const& ptBinsEdge,
                    std::vector<float> const& dcaV0ToPv,
                    std::vector<float> const& decayLength,
                    std::vector<float> const& cosTheta,
                    std::vector<float> const& dcaDaughters,
                    std::vector<float> const& kDca,
                    std::vector<float> const& pDca,
                    std::vector<float> const& xDca):
                    xMassHypothesis(xMassHypothesis),
                    minMass(minMass),
                    maxMass(maxMass),
                    rapidityCut(rapidityCut),
                    nPtBins(nPtBins),
                    ptBinsEdge(ptBinsEdge),
                    dcaV0ToPv(dcaV0ToPv),
                    decayLength(decayLength),
                    cosTheta(cosTheta),
                    dcaDaughters(dcaDaughters),
                    kDca(kDca),
                    pDca(pDca),
                    xDca(xDca)
    {}
  };

  TopologicalCuts const DpmCuts(M_PION_PLUS, 1.843, 1.898,
                          1.0,                                      // rapidity
                          5,                                        // nPtBins
                          {0.    , 1.,     2.,     3.,     5., 15.},// ptEdges
                          {0.0061, 0.0049, 0.0038, 0.0038, 0.0040}, // dcaV0ToPv
                          {0.0300, 0.0300, 0.0300, 0.0247, 0.0259}, // decayLength
                          {0.9500, 0.9500, 0.9500, 0.9500, 0.9500}, // cosTheta
                          {0.0084, 0.0066, 0.0057, 0.0050, 0.0060}, // dcaDaughters
                          {0.0103, 0.0091, 0.0095, 0.0079, 0.0058}, // kDca
                          {0.0110, 0.0111, 0.0086, 0.0081, 0.0062}, // p1Dca
                          {0.0110, 0.0111, 0.0086, 0.0081, 0.0062});// p2Dca

  TopologicalCuts const DsCuts(M_KAON_MINUS, 1.940, 1.995,
                          1.0,                                      // rapidity
                          5,                                        // nPtBins
                          {0.    , 1.,     2.,     3.,     5., 15.},// ptEdges
                          {0.0061, 0.0049, 0.0038, 0.0038, 0.0040}, // dcaV0ToPv
                          {0.0150, 0.0150, 0.0150, 0.0150, 0.0150}, // decayLength
                          {0.9500, 0.9500, 0.9500, 0.9500, 0.9500}, // cosTheta
                          {0.0084, 0.0066, 0.0057, 0.0050, 0.0060}, // dcaDaughters
                          {0.0103, 0.0091, 0.0095, 0.0079, 0.0058}, // k1Dca
                          {0.0110, 0.0111, 0.0086, 0.0081, 0.0062}, // pDca
                          {0.0103, 0.0091, 0.0095, 0.0079, 0.0058});// k2Dca

  TopologicalCuts const LcCuts(M_PROTON, 2.258, 2.308,
                          1.0,                                      // rapidity
                          5,                                        // nPtBins
                          {0.    , 1.,     2.,     3.,     5., 15.},// ptEdges
                          {0.0100, 0.0100, 0.0100, 0.0100, 0.0100}, // dcaV0ToPv
                          {0.0080, 0.0236, 0.0236, 0.0250, 0.0250}, // decayLength
                          {0.9900, 0.9915, 0.9915, 0.9950, 0.9950}, // cosTheta
                          {0.0060, 0.0086, 0.0086, 0.0050, 0.0050}, // dcaDaughters
                          {0.0070, 0.0114, 0.0114, 0.0075, 0.0075}, // k1Dca
                          {0.0070, 0.0076, 0.0076, 0.0080, 0.0080}, // piDca
                          {0.0070, 0.0105, 0.0105, 0.0065, 0.0065});// pDhca
}
#endif
