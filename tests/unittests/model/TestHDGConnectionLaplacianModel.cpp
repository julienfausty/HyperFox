#include <catch2/catch.hpp>
#include "HDGConnectionLaplacianModel.h"
#include "HDGLaplaceModel.h"

using namespace hfox;

TEST_CASE("Testing the HDGConnectionLaplacianModel", "[unit][model][HDGEmbeddedModel][HDGConnectionLaplacianModel]"){

  double tol = 1e-12;
  int maxDim = 3;
  int maxOrd = 5;

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
        CHECK((diffMat*(diffMat.transpose())).trace() < tol);
      };
    }
  }

  maxOrd = 10;

  for(int iO = 1; iO < maxOrd+1; iO++){
    SECTION("Testing embedded x surface (ord=" + std::to_string(iO) + ")"){
      int nVal = 1;
      int refDim = 2;
      int eDim = 3;
      ReferenceElement refEl(refDim, iO, "simplex");
      HDGConnectionLaplacianModel mod(&refEl);
      mod.setEmbeddingDimension(eDim);
      mod.allocate(1);
      int nNodesEl = refEl.getNumNodes();
      int nFaces = refEl.getNumFaces();
      int nNodesFc = refEl.getFaceElement()->getNumNodes();
      std::vector< std::vector<double> > nodes(nNodesEl, std::vector<double>(eDim, 0.0));
      for(int iN = 0; iN < nNodesEl; iN++){
        nodes[iN][0] = refEl.getNodes()->at(iN)[0];
        nodes[iN][1] = refEl.getNodes()->at(iN)[1];
        nodes[iN][2] = std::pow(refEl.getNodes()->at(iN)[0], nVal);
      }
      mod.setElementNodes(&nodes);
      std::map<std::string, std::vector<double> > fm;
      fm["Tau"] = std::vector<double>(nNodesFc*nFaces, 1.0);
      fm["Jacobian"] = std::vector<double>(nNodesEl*refDim*eDim, 0.0);
      fm["Metric"] = std::vector<double>(nNodesEl*refDim*refDim, 0.0);
      fm["DiffusionTensor"] = std::vector<double>(nNodesEl*eDim*eDim, 0.0);
      for(int iN = 0; iN < nNodesEl; iN++){
        EMap<EMatrix> jac(fm["Jacobian"].data() + iN*refDim*eDim, refDim, eDim);
        jac.block(0, 0, refDim, refDim) = EMatrix::Identity(refDim, refDim);
        jac(0, 2) = nVal*std::pow(nodes[iN][0], nVal-1);
        EMap<EMatrix> metric(fm["Metric"].data() + iN*refDim*refDim, refDim, refDim);
        metric = EMatrix::Identity(refDim, refDim);
        metric(0, 0) += std::pow(nVal * std::pow(nodes[iN][0], nVal-1), 2);
        EMap<EMatrix>(fm["DiffusionTensor"].data() + iN*eDim*eDim, eDim, eDim) = EMatrix::Identity(eDim, eDim);
      }
      mod.setFieldMap(&fm);
      CHECK_NOTHROW(mod.compute());
      std::vector<double> J(nNodesEl*eDim, 0.0);
      std::vector<double> psi(nNodesEl, 0.0);
      for(int iN = 0; iN < nNodesEl; iN++){
        J[iN*eDim] = nodes[iN][0];
        psi[iN] = (std::pow(1.0 + std::pow(nVal*std::pow(nodes[iN][0], nVal-1), 2), 2));
      }
      EMatrix Suq = mod.getLocalMatrix()->block(0, nNodesEl, nNodesEl, nNodesEl*eDim);
      double computedInt = (EMap<EVector>(psi.data(), nNodesEl).transpose() * Suq * EMap<EVector>(J.data(), nNodesEl*eDim))(0,0);
      double anaInt = 0.0;
      for(int ip = 0; ip < refEl.getNumIPs(); ip++){
        double x = refEl.getIPCoords()->at(ip)[0];
        anaInt += std::pow(nVal, 2)*(nVal-1)*std::pow(x, 3*(nVal-1))*std::sqrt(1.0 + std::pow(nVal*std::pow(x, nVal-1), 2)) * (refEl.getIPWeights()->at(ip));
        anaInt -= (1 + std::pow(nVal * std::pow(x, nVal-1), 2))*std::sqrt(1.0 + std::pow(nVal*std::pow(x, nVal-1), 2)) * (refEl.getIPWeights()->at(ip));
      }
      CHECK(std::abs(computedInt - anaInt) < tol);
      std::fill(J.begin(), J.end(), 0.0);
      std::fill(psi.begin(), psi.end(), 0.0);
      std::vector<double> lambda(nNodesFc*nFaces, 0.0);
      for(int iN = 0; iN < nNodesEl; iN++){
        double x = nodes[iN][0];
        J[iN*eDim] = 1.0/(1.0 + std::pow(nVal, 2)*std::pow(x, 2*nVal-2));
        J[iN*eDim + 2] = nVal*std::pow(x, nVal-1)*J[iN*eDim];
        psi[iN] = x;
      }
      for(int iF = 0; iF < nFaces; iF++){
        for(int iN = 0; iN < nNodesFc; iN++){
          double x = nodes[refEl.getFaceNodes()->at(iF)[iN]][0];
          lambda[iF*nNodesFc + iN] = x;
        }
      }
      EMatrix Squ = mod.getLocalMatrix()->block(nNodesEl, 0, nNodesEl*eDim, nNodesEl);
      EMatrix Sqq = mod.getLocalMatrix()->block(nNodesEl, nNodesEl, nNodesEl*eDim, nNodesEl*eDim);
      EMatrix Sql = mod.getLocalMatrix()->block(nNodesEl, nNodesEl*(eDim+1), nNodesEl*eDim, nNodesFc*nFaces);
      EVector Diff = Squ*EMap<EVector>(psi.data(), nNodesEl) + Sqq*EMap<EVector>(J.data(), nNodesEl*eDim) + Sql*EMap<EVector>(lambda.data(), nFaces*nNodesFc);
      CHECK(Diff.norm() < tol);
    };
  }

  for(int iO = 4; iO < maxOrd+1; iO++){
    SECTION("Testing embedded x^2 surface (ord=" + std::to_string(iO) + ")"){
      int nVal = 2;
      int refDim = 2;
      int eDim = 3;
      ReferenceElement refEl(refDim, iO, "simplex");
      HDGConnectionLaplacianModel mod(&refEl);
      mod.setEmbeddingDimension(eDim);
      mod.allocate(1);
      int nNodesEl = refEl.getNumNodes();
      int nFaces = refEl.getNumFaces();
      int nNodesFc = refEl.getFaceElement()->getNumNodes();
      std::vector< std::vector<double> > nodes(nNodesEl, std::vector<double>(eDim, 0.0));
      for(int iN = 0; iN < nNodesEl; iN++){
        nodes[iN][0] = refEl.getNodes()->at(iN)[0];
        nodes[iN][1] = refEl.getNodes()->at(iN)[1];
        nodes[iN][2] = std::pow(refEl.getNodes()->at(iN)[0], nVal);
      }
      mod.setElementNodes(&nodes);
      std::map<std::string, std::vector<double> > fm;
      fm["Tau"] = std::vector<double>(nNodesFc*nFaces, 1.0);
      fm["Jacobian"] = std::vector<double>(nNodesEl*refDim*eDim, 0.0);
      fm["Metric"] = std::vector<double>(nNodesEl*refDim*refDim, 0.0);
      fm["DiffusionTensor"] = std::vector<double>(nNodesEl*eDim*eDim, 0.0);
      for(int iN = 0; iN < nNodesEl; iN++){
        EMap<EMatrix> jac(fm["Jacobian"].data() + iN*refDim*eDim, refDim, eDim);
        jac.block(0, 0, refDim, refDim) = EMatrix::Identity(refDim, refDim);
        jac(0, 2) = nVal*std::pow(nodes[iN][0], nVal-1);
        EMap<EMatrix> metric(fm["Metric"].data() + iN*refDim*refDim, refDim, refDim);
        metric = EMatrix::Identity(refDim, refDim);
        metric(0, 0) += std::pow(nVal * std::pow(nodes[iN][0], nVal-1), 2);
        EMap<EMatrix>(fm["DiffusionTensor"].data() + iN*eDim*eDim, eDim, eDim) = EMatrix::Identity(eDim, eDim);
      }
      mod.setFieldMap(&fm);
      CHECK_NOTHROW(mod.compute());
      std::vector<double> J(nNodesEl*eDim, 0.0);
      std::vector<double> psi(nNodesEl, 0.0);
      for(int iN = 0; iN < nNodesEl; iN++){
        J[iN*eDim] = 1.0;
        psi[iN] = std::pow(1.0 + std::pow(nVal*std::pow(nodes[iN][0], nVal-1), 2), 2);
      }
      EMatrix Suq = mod.getLocalMatrix()->block(0, nNodesEl, nNodesEl, nNodesEl*eDim);
      double computedInt = (EMap<EVector>(psi.data(), nNodesEl).transpose() * Suq * EMap<EVector>(J.data(), nNodesEl*eDim))(0,0);
      double anaInt = 0.0;
      for(int ip = 0; ip < refEl.getNumIPs(); ip++){
        double x = refEl.getIPCoords()->at(ip)[0];
        anaInt += std::pow(nVal, 2)*(nVal-1.0)*std::pow(x, 2*nVal-3)*std::sqrt(1.0 + std::pow(nVal*std::pow(x, nVal-1), 2)) * (refEl.getIPWeights()->at(ip));
      }
      CHECK(std::abs(computedInt - anaInt) < 1e-2);
      std::fill(J.begin(), J.end(), 0.0);
      std::fill(psi.begin(), psi.end(), 0.0);
      std::vector<double> lambda(nNodesFc*nFaces, 0.0);
      for(int iN = 0; iN < nNodesEl; iN++){
        double x = nodes[iN][0];
        J[iN*eDim] = 1.0;
        J[iN*eDim + 2] = nVal*std::pow(x, nVal-1);
        psi[iN] = (x + (std::pow(nVal, 2)/(2*nVal-1))*std::pow(x, 2*nVal-1));
      }
      for(int iF = 0; iF < nFaces; iF++){
        for(int iN = 0; iN < nNodesFc; iN++){
          double x = nodes[refEl.getFaceNodes()->at(iF)[iN]][0];
          lambda[iF*nNodesFc + iN] = (x + (std::pow(nVal, 2)/(2*nVal-1))*std::pow(x, 2*nVal-1));
        }
      }
      EMatrix Squ = mod.getLocalMatrix()->block(nNodesEl, 0, nNodesEl*eDim, nNodesEl);
      EMatrix Sqq = mod.getLocalMatrix()->block(nNodesEl, nNodesEl, nNodesEl*eDim, nNodesEl*eDim);
      EMatrix Sql = mod.getLocalMatrix()->block(nNodesEl, nNodesEl*(eDim+1), nNodesEl*eDim, nNodesFc*nFaces);
      EVector Diff = Squ*EMap<EVector>(psi.data(), nNodesEl) + Sqq*EMap<EVector>(J.data(), nNodesEl*eDim) + Sql*EMap<EVector>(lambda.data(), nFaces*nNodesFc);
      CHECK(Diff.norm() < 2e-2);
    };
  }
};
