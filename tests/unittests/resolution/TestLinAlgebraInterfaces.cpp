#include <catch2/catch.hpp>
#include <string>
#include <map>
#include <tuple>
#include "LinAlgebraInterface.h"
#include "PetscInterface.h"

using namespace hfox;

TEST_CASE("Stock linear systems", "[unit][resolution][LinAlgebraInterface]"){
  //setup
  std::map<std::string, LinAlgebraInterface*> interfaces;
  //interfaces["Petsc"] = new PetscInterface();
  
  std::map<std::string, LinAlgebraInterface*>::iterator itIFace;
  for(itIFace = interfaces.begin(); itIFace != interfaces.end(); itIFace++){
    std::string name = itIFace->first;
    LinAlgebraInterface * iFace = itIFace->second;
    SECTION(name + ": checking initialize, configure and allocate ordering"){
      int i = 0, j = 0, ndofs = 1;
      double val = 1;
      std::vector<double> doubleVec = {1};
      std::vector<int> intVec = {0};
      std::vector< std::tuple<int, int> > tupleVec = {std::tuple<int, int>(0, 0)};
      //constructed interface
      CHECK_THROWS(iFace->solve(&doubleVec));
      CHECK_THROWS(iFace->setValsRHS(intVec, doubleVec));
      CHECK_THROWS(iFace->setValRHS(i, val));
      CHECK_THROWS(iFace->setValsMatrix(tupleVec, doubleVec));
      CHECK_THROWS(iFace->setValMatrix(i, j, val));
      CHECK_THROWS(iFace->allocate(ndofs));
      CHECK_THROWS(iFace->configure());

      CHECK_NOTHROW(iFace->initialize());
      //initialized interface
      CHECK_THROWS(iFace->solve(&doubleVec));
      CHECK_THROWS(iFace->setValsRHS(intVec, doubleVec));
      CHECK_THROWS(iFace->setValRHS(i, val));
      CHECK_THROWS(iFace->setValsMatrix(tupleVec, doubleVec));
      CHECK_THROWS(iFace->setValMatrix(i, j, val));
      CHECK_THROWS(iFace->allocate(ndofs));

      CHECK_NOTHROW(iFace->configure());
      //configured interface
      CHECK_THROWS(iFace->solve(&doubleVec));
      CHECK_THROWS(iFace->setValsRHS(intVec, doubleVec));
      CHECK_THROWS(iFace->setValRHS(i, val));
      CHECK_THROWS(iFace->setValsMatrix(tupleVec, doubleVec));
      CHECK_THROWS(iFace->setValMatrix(i, j, val));

      CHECK_NOTHROW(iFace->allocate(ndofs));
      //allocated interface
      CHECK_NOTHROW(iFace->setValsRHS(intVec, doubleVec));
      CHECK_NOTHROW(iFace->setValRHS(i, val));
      CHECK_NOTHROW(iFace->setValsMatrix(tupleVec, doubleVec));
      CHECK_NOTHROW(iFace->setValMatrix(i, j, val));
      CHECK_NOTHROW(iFace->solve(&doubleVec));

      //empty system
      CHECK_NOTHROW(iFace->clearSystem());
    }

    SECTION(name + ": Identity solving"){
      //insert bigger and bigger identity matrix problems here
    }
  }

  for(itIFace = interfaces.begin(); itIFace != interfaces.end(); itIFace++){
    delete (itIFace->second);
  }
};
