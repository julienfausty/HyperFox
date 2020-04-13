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
  interfaces["Petsc"] = new PetscInterface();
  
  std::map<std::string, LinAlgebraInterface*>::iterator itIFace;
  for(itIFace = interfaces.begin(); itIFace != interfaces.end(); itIFace++){
    std::string name = itIFace->first;
    LinAlgebraInterface * iFace = itIFace->second;
    SECTION(name + ": checking initialize, configure and allocate ordering"){
      int i = 0, j = 0, ndofs = 1;
      double val = 1;
      std::vector<double> doubleVec = {1};
      std::vector<int> intVec = {0};
      //constructed interface
      CHECK_THROWS(iFace->solve(&doubleVec));
      CHECK_THROWS(iFace->assemble());
      CHECK_THROWS(iFace->setValsRHS(intVec, doubleVec));
      CHECK_THROWS(iFace->setValRHS(i, val));
      CHECK_THROWS(iFace->setValsMatrix(intVec, intVec, doubleVec));
      CHECK_THROWS(iFace->setValMatrix(i, j, val));
      CHECK_THROWS(iFace->allocate(ndofs));
      CHECK_THROWS(iFace->configure());

      CHECK_NOTHROW(iFace->initialize());
      //initialized interface
      CHECK_THROWS(iFace->solve(&doubleVec));
      CHECK_THROWS(iFace->assemble());
      CHECK_THROWS(iFace->setValsRHS(intVec, doubleVec));
      CHECK_THROWS(iFace->setValRHS(i, val));
      CHECK_THROWS(iFace->setValsMatrix(intVec, intVec, doubleVec));
      CHECK_THROWS(iFace->setValMatrix(i, j, val));
      CHECK_THROWS(iFace->allocate(ndofs));

      CHECK_NOTHROW(iFace->configure());
      //configured interface
      CHECK_THROWS(iFace->solve(&doubleVec));
      CHECK_THROWS(iFace->assemble());
      CHECK_THROWS(iFace->setValsRHS(intVec, doubleVec));
      CHECK_THROWS(iFace->setValRHS(i, val));
      CHECK_THROWS(iFace->setValsMatrix(intVec, intVec, doubleVec));
      CHECK_THROWS(iFace->setValMatrix(i, j, val));

      CHECK_NOTHROW(iFace->allocate(ndofs));
      //allocated interface
      CHECK_NOTHROW(iFace->setValsRHS(intVec, doubleVec));
      CHECK_NOTHROW(iFace->setValRHS(i, val));
      CHECK_NOTHROW(iFace->setValsMatrix(intVec, intVec, doubleVec));
      CHECK_NOTHROW(iFace->setValMatrix(i, j, val));
      CHECK_NOTHROW(iFace->assemble());
      CHECK_NOTHROW(iFace->solve(&doubleVec));

      ////empty system
      //CHECK_NOTHROW(iFace->clearSystem());
    }

    //SECTION(name + ": Identity solving"){
      //int maxDofs = 1e2;
      //std::vector<int> is;
      //std::vector<int> js;
      //std::vector<double> vals, rhs, sol;
      //for(int k = 0; k < maxDofs; k++){
        //iFace->clearSystem();
        //iFace->allocate(k+1);
        //vals.push_back(1.0);
        //is.push_back(k);
        //js.push_back(k);
        //rhs.push_back(k);
        //for(int i = 0; i < is.size(); i++){
          //iFace->setValMatrix(is[i], js[i], vals[i]);
          //iFace->setValRHS(is[i], rhs[i]);
        //}
        //iFace->assemble();
        //CHECK_NOTHROW(iFace->solve(&sol));
        //for(int j = 0; j < k+1; j++){   
          //CHECK(sol[j] == Approx(rhs[j]).margin(1e-12));
        //}
      //}
    //};

    //SECTION(name + ": Triangular matrix"){
      //int maxDofs = 1e2;
      //std::vector<int> is;
      //std::vector<int> js, ls;
      //std::vector<double> vals, rhs, sol, anaSol;
      //for(int k = 0; k < maxDofs; k++){
        //iFace->clearSystem();
        //iFace->allocate(k+1);
        //vals.push_back(1.0);
        //is.push_back(k);
        //js.push_back(k);
        //for(int j = 0; j < k; j++){
          //vals.push_back(1.0);
          //is.push_back(k);
          //js.push_back(j);
        //}
        //ls.push_back(k);
        //rhs.push_back(1.0);
        //anaSol.resize(k+1, 0.0);
        //anaSol[0] = 1.0;
        //for(int i = 0; i < is.size(); i++){
          //iFace->setValMatrix(is[i], js[i], vals[i]);
          //iFace->setValRHS(is[i], rhs[i]);
        //}
        //iFace->assemble();
        //CHECK_NOTHROW(iFace->solve(&sol));
        //for(int j = 0; j < k+1; j++){   
          //CHECK(sol[j] == Approx(anaSol[j]).margin(1e-12));
        //}
      //}
    //};

    //SECTION(name + ": Hinge matrix"){
      //int maxDofs = 1e2;
      //std::vector<int> is, ls;
      //std::vector<int> js;
      //std::vector<double> vals, rhs, sol, anaSol;
      //for(int k = 0; k < maxDofs; k++){
        //iFace->clearSystem();
        //iFace->allocate(k+1);
        //vals.push_back(1.0);
        //is.push_back(k);
        //js.push_back(k);
        //if(k != 0){
          //vals.push_back(1.0);
          //is.push_back(k);
          //js.push_back(k-1);
        //}
        //ls.push_back(k);
        //rhs.push_back(1.0);
        //anaSol.resize(k+1, 0.0);
        //if((k % 2) == 0){
          //anaSol.push_back(1.0);
        //} else {
          //anaSol.push_back(0.0);
        //}
        //for(int i = 0; i < is.size(); i++){
          //iFace->setValMatrix(is[i], js[i], vals[i]);
          //iFace->setValRHS(is[i], rhs[i]);
        //}
        //iFace->assemble();
        //CHECK_NOTHROW(iFace->solve(&sol));
        //for(int j = 0; j < k+1; j++){   
          //CHECK(sol[j] == Approx(anaSol[j]).margin(1e-12));
        //}
      //}
    //};
  }



  for(itIFace = interfaces.begin(); itIFace != interfaces.end(); itIFace++){
    delete (itIFace->second);
  }
};
