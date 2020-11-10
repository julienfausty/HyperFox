#include <catch2/catch.hpp>
#include "HDGnGammaModel.h"
#include "HDGConvection.h"
#include "HDGDiffusion.h"
#include "Reaction.h"
#include "Source.h"
#include "Euler.h"

using namespace hfox;
using namespace nGamma;

TEST_CASE("Testing the HDGnGammaModel", "[unit][model][HDGBurgersModel]"){

  double tol = 1e-12;
  int maxOrd = 7;

  for(int j = 0; j < maxOrd; j++){
    ReferenceElement refEl(2, j+1, "simplex");
    SECTION("Testing setup (order=" + std::to_string(j+1) + ")"){
      CHECK_NOTHROW(HDGnGammaModel(&refEl));
      HDGnGammaModel mod(&refEl);
      CHECK_THROWS(mod.compute());
      CHECK_THROWS(mod.allocate(3));
      CHECK_NOTHROW(mod.allocate(2));
      CHECK_THROWS(mod.compute());
      CHECK_NOTHROW(mod.setElementNodes(refEl.getNodes()));
      std::map<std::string, std::vector<double> > fm;
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["Solution"] = std::vector<double>(refEl.getNumNodes()*2, 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["Tau"] = std::vector<double>(refEl.getNumFaces() * (refEl.getFaceElement()->getNumNodes()) * std::pow(2, 2), 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["Trace"] = std::vector<double>(refEl.getNumFaces() * (refEl.getFaceElement()->getNumNodes()) * 2, 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["b"] = std::vector<double>(refEl.getNumNodes()*3, 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["D"] = std::vector<double>(refEl.getNumNodes()*4, 1.0);
      CHECK_THROWS(mod.setFieldMap(&fm));
      fm["G"] = std::vector<double>(refEl.getNumNodes()*4, 1.0);
      CHECK_NOTHROW(mod.setFieldMap(&fm));
    };

    //setup field map
    std::map<std::string, std::vector<double> > fm;
    std::vector<double> tau(refEl.getNumFaces() * (refEl.getFaceElement()->getNumNodes()) * std::pow(2, 2), 0.0);
    std::vector<double> trace(refEl.getNumFaces() * (refEl.getFaceElement()->getNumNodes()) * 2, 0.0);
    std::vector<double> sol(refEl.getNumNodes() * 2, 2.0);
    std::vector<double> diff(refEl.getNumNodes()*4, 0);
    std::vector<double> b(refEl.getNumNodes() * 3, 1.0);
    for(int k = 0; k < refEl.getNumNodes(); k++){
      EMap<EMatrix>(diff.data() + 4*k, 2, 2) = EMatrix::Identity(2, 2)*3.0;
    }
    for(int k = 0; k < refEl.getNumFaces()*refEl.getFaceElement()->getNumNodes(); k++){
      EMap<EMatrix>(tau.data() + 4*k, 2, 2) = EMatrix::Identity(2, 2);
    }
    int nFaces = refEl.getNumFaces();
    int nNodesFc = refEl.getFaceElement()->getNumNodes();
    for(int iFace = 0; iFace < nFaces; iFace++){
      for(int iN = 0; iN < nNodesFc; iN++){
        for(int d = 0; d < 2; d++){
          trace[(iFace*nNodesFc + iN)*2 + d] = 2.0;
        }
      }
    }
    fm["Tau"] = tau;
    fm["Solution"] = sol;
    fm["Trace"] = trace;
    fm["D"] = diff;
    fm["G"] = diff;
    fm["b"] = b;

    SECTION("Testing stationary assembly (order = " + std::to_string(j+1) + ")"){
      HDGnGammaModel mod(&refEl);
      mod.allocate(2);
      mod.setFieldMap(&fm);
    };
  }

};
