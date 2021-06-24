#include <catch2/catch.hpp>
#include <string>
#include <cmath>
#include "IP2NodeSolver.h"

using namespace hfox;

void monomials(int dim, int ord, std::vector<double> & nodes, std::vector<double> * vals){
  int nNodes = nodes.size()/dim;
  vals->resize(nNodes, 0.0);
  std::fill(vals->begin(), vals->end(), 0.0);
  for(int iN = 0; iN < nNodes; iN++){
    vals->at(iN) = std::pow(nodes[iN*dim], ord);
  }
};

TEST_CASE("Testing the IP2NodeSolver", "[unit][geometry][IP2NodeSolver]"){

  SECTION("Testing construction and error handling"){
    Mesh myMesh(1, 1, "simplex");
    CHECK_NOTHROW(IP2NodeSolver(&myMesh));
    IP2NodeSolver ip2n(&myMesh);
    CHECK_THROWS(ip2n.solve());
    Field faceField(&myMesh, Face, 1, 1);
    CHECK_THROWS(ip2n.setField(&faceField));
    Field cellField(&myMesh, Cell, 1, 1);
    CHECK_NOTHROW(ip2n.setField(&cellField));
    CHECK_THROWS(ip2n.solve());
  };

  int maxDim = 3;
  int maxOrder = 5;
  double tol = 1e-6;
  SECTION("Testing monomial reconstruction on nodes"){
    for(int dim = 2; dim < maxDim + 1; dim++){
      for(int ord = 1; ord < maxOrder + 1; ord++){
        Mesh myMesh(dim, ord, "simplex");
        const ReferenceElement * refEl = myMesh.getReferenceElement();
        int nNodesEl = refEl->getNumNodes();
        std::vector<double> nodeCoords(nNodesEl*dim, 0.0);
        for(int iN = 0; iN < nNodesEl; iN++){
          std::copy(refEl->getNodes()->at(iN).begin(), refEl->getNodes()->at(iN).end(), nodeCoords.begin() + iN*dim);
        }
        std::vector<int> connectivity(nNodesEl, 0);
        std::iota(connectivity.begin(), connectivity.end(), 0);
        myMesh.setMesh(dim, nodeCoords, connectivity);
        Field res(&myMesh, Cell, nNodesEl, 1);
        IP2NodeSolver ip2n(&myMesh);
        ip2n.setField(&res);
        ip2n.setFunction([dim, ord](std::vector<double> & xs, std::vector<double> * vals){monomials(dim,ord,xs,vals);});
        ip2n.solve();
        std::vector<double> realVals(nNodesEl, 0.0);
        monomials(dim, ord, nodeCoords, &realVals);
        for(int iN = 0; iN < nNodesEl; iN++){
          CHECK(std::abs(realVals[iN] - res.getValues()->at(iN)) < tol);
        }
      }
    }
  }

  SECTION("Testing monomial surface properties on nodes"){
    int dim = 2;
    for(int ord = 1; ord < maxOrder + 1; ord++){
      Mesh myMesh(dim, ord, "simplex");
      const ReferenceElement * refEl = myMesh.getReferenceElement();
      int nNodesEl = refEl->getNumNodes();
      std::vector<double> nodeCoords(nNodesEl*(dim+1), 0.0);
      for(int iN = 0; iN < nNodesEl; iN++){
        std::copy(refEl->getNodes()->at(iN).begin(), refEl->getNodes()->at(iN).end(), nodeCoords.begin() + iN*(dim+1));
      }
      std::vector<double> monoms(nNodesEl, 0.0);
      monomials(dim + 1, ord, nodeCoords, &monoms);
      for(int iN = 0; iN < nNodesEl; iN++){
        nodeCoords[iN*(dim+1) + 2] = monoms[iN];
      }
      std::vector<int> connectivity(nNodesEl, 0);
      std::iota(connectivity.begin(), connectivity.end(), 0);
      myMesh.setMesh(dim + 1, nodeCoords, connectivity);
      Field res(&myMesh, Cell, nNodesEl, dim*(dim+1));
      IP2NodeSolver ip2n(&myMesh);
      ip2n.setField(&res);
      ip2n.setFunction(IP2NodeSolver::jacobianVals);
      ip2n.solve();
      std::vector<EMatrix> jacs(nNodesEl, EMatrix::Zero(dim, dim+1));
      monomials(dim + 1, ord-1, nodeCoords, &monoms);
      for(int iN = 0; iN < nNodesEl; iN++){
        jacs[iN].block(0, 0, dim, dim) = EMatrix::Identity(dim, dim);
        jacs[iN](0, 2) = ord*monoms[iN];
      }
      for(int iN = 0; iN < nNodesEl; iN++){
        CHECK((EMap<EMatrix>(res.getValues()->data() + iN*dim*(dim+1), dim, dim+1) - jacs[iN]).sum() < tol);
      }
    }
  }

  SECTION("Testing monomial metric properties on nodes"){
    int dim = 2;
    for(int ord = 1; ord < maxOrder + 1; ord++){
      int refOrd = 2*(ord-1);
      if(refOrd == 0){
        refOrd = 1;
      }
      Mesh myMesh(dim, refOrd, "simplex");
      const ReferenceElement * refEl = myMesh.getReferenceElement();
      int nNodesEl = refEl->getNumNodes();
      std::vector<double> nodeCoords(nNodesEl*(dim+1), 0.0);
      for(int iN = 0; iN < nNodesEl; iN++){
        std::copy(refEl->getNodes()->at(iN).begin(), refEl->getNodes()->at(iN).end(), nodeCoords.begin() + iN*(dim+1));
      }
      std::vector<double> monoms(nNodesEl, 0.0);
      monomials(dim + 1, ord, nodeCoords, &monoms);
      for(int iN = 0; iN < nNodesEl; iN++){
        nodeCoords[iN*(dim+1) + 2] = monoms[iN];
      }
      std::vector<int> connectivity(nNodesEl, 0);
      std::iota(connectivity.begin(), connectivity.end(), 0);
      myMesh.setMesh(dim + 1, nodeCoords, connectivity);
      Field res(&myMesh, Cell, nNodesEl, dim*dim);
      IP2NodeSolver ip2n(&myMesh);
      ip2n.setField(&res);
      ip2n.setFunction(IP2NodeSolver::metricVals);
      ip2n.solve();
      std::vector<EMatrix> gs(nNodesEl, EMatrix::Zero(dim, dim));
      monomials(dim + 1, 2*(ord-1), nodeCoords, &monoms);
      for(int iN = 0; iN < nNodesEl; iN++){
        gs[iN] = EMatrix::Identity(dim, dim);
        gs[iN](0,0) += std::pow(ord, 2)*monoms[iN];
      }
      for(int iN = 0; iN < nNodesEl; iN++){
        CHECK((EMap<EMatrix>(res.getValues()->data() + iN*dim*dim, dim, dim) - gs[iN]).sum() < tol);
      }
    }
  }

};
