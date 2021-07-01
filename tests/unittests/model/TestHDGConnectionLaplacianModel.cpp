#include <catch2/catch.hpp>
#include "HDGConnectionLaplacianModel.h"
#include "HDGLaplaceModel.h"

using namespace hfox;

TEST_CASE("Testing the HDGConnectionLaplacianModel", "[unit][model][HDGEmbeddedModel][HDGConnectionLaplacianModel]"){

  double tol = 1e-12;
  int maxDim = 2;
  int maxOrd = 3;

  for(int iD = 2; iD < maxDim+1; iD++){
    for(int iO = 1; iO < maxOrd+1; iO++){

      ReferenceElement refEl(iD, iO, "simplex");
      int nNodesEl = refEl.getNumNodes();
      int nFaces = refEl.getNumFaces();
      const ReferenceElement * fEl = refEl.getFaceElement();
      int nNodesFc = fEl->getNumNodes();
      std::vector< std::vector<double> > nodes(nNodesEl);
      SECTION("Testing construction and error handling (refDim = " + std::to_string(iD) + ", ord = " + std::to_string(iO) + ")"){
        CHECK_NOTHROW(HDGConnectionLaplacianModel(&refEl));
        for(int eD = 0; eD < iD; eD++){
          HDGConnectionLaplacianModel mod(&refEl);
          CHECK_THROWS(mod.setEmbeddingDimension(eD));
        }
        for(int eD = iD; eD < maxDim + 1; eD++){
          HDGConnectionLaplacianModel mod(&refEl);
          CHECK_NOTHROW(mod.setEmbeddingDimension(eD));
          CHECK_THROWS(mod.compute());
          CHECK_NOTHROW(mod.allocate(1));
          for(int iN = 0; iN < nNodesEl; iN++){
            nodes[iN].resize(eD, 0.0);
            std::fill(nodes[iN].begin(), nodes[iN].end(), 0.0);
            std::copy(refEl.getNodes()->at(iN).begin(), refEl.getNodes()->at(iN).end(), nodes[iN].begin());
          }
          CHECK_NOTHROW(mod.setElementNodes(&nodes));
          std::map<std::string, std::vector<double> > fm;
          CHECK_THROWS(mod.setFieldMap(&fm));
          fm["Tau"] = std::vector<double>(nFaces * nNodesFc, 1.0);
          CHECK_THROWS(mod.setFieldMap(&fm));
          fm["Jacobian"] = std::vector<double>(nNodesEl * eD * iD, 1.0);
          CHECK_THROWS(mod.setFieldMap(&fm));
          fm["Metric"] = std::vector<double>(nNodesEl * iD * iD, 1.0);
          CHECK_NOTHROW(mod.setFieldMap(&fm));
          fm["DiffusionTensor"] = std::vector<double>(nNodesEl * eD * eD, 1.0);
          CHECK_NOTHROW(mod.setFieldMap(&fm));
        }
      };

      SECTION("Testing reduction to Laplace model for embedding dimension = reference dimension (refDim = " + std::to_string(iD) + ", ord = " + std::to_string(iO) + ")"){
        HDGConnectionLaplacianModel embeddedMod(&refEl);
        embeddedMod.setEmbeddingDimension(iD);
        HDGLaplaceModel lapMod(&refEl);
        for(int iN = 0; iN < nNodesEl; iN++){
          nodes[iN].resize(iD, 0.0);
          std::fill(nodes[iN].begin(), nodes[iN].end(), 0.0);
          std::copy(refEl.getNodes()->at(iN).begin(), refEl.getNodes()->at(iN).end(), nodes[iN].begin());
        }
        embeddedMod.allocate(1);
        lapMod.allocate(1);
        embeddedMod.setElementNodes(&nodes);
        lapMod.setElementNodes(&nodes);
        std::map<std::string, std::vector<double> > fm;
        fm["Tau"] = std::vector<double>(nNodesFc*nFaces, 1.0);
        fm["Jacobian"] = std::vector<double>(nNodesEl*iD*iD, 0.0);
        fm["Metric"] = std::vector<double>(nNodesEl*iD*iD, 0.0);
        fm["DiffusionTensor"] = std::vector<double>(nNodesEl*iD*iD, 0.0);
        for(int iN = 0; iN < nNodesEl; iN++){
          EMap<EMatrix>(fm["Jacobian"].data() + iN*iD*iD, iD, iD) = EMatrix::Identity(iD, iD);
          EMap<EMatrix>(fm["Metric"].data() + iN*iD*iD, iD, iD) = EMatrix::Identity(iD, iD);
          EMap<EMatrix>(fm["DiffusionTensor"].data() + iN*iD*iD, iD, iD) = EMatrix::Identity(iD, iD);
        }
        embeddedMod.setFieldMap(&fm);
        lapMod.setFieldMap(&fm);
        CHECK_NOTHROW(embeddedMod.compute());
        CHECK_NOTHROW(lapMod.compute());
        EMatrix diffMat = *(embeddedMod.getLocalMatrix()) - *(lapMod.getLocalMatrix());
        std::cout << "Difference" << std::endl;
        std::cout << diffMat << std::endl;
        CHECK((diffMat*(diffMat.transpose())).trace() < tol);
      };
    }
  }
};
