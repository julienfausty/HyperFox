#include <catch2/catch.hpp>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include "FieldIntegrator.h"
#include "Mesh.h"
#include "HDF5Io.h"
#include "TestUtils.h"

using namespace hfox;

double integralMonomial(int power){
  return 1.0/(power + 1);
};

TEST_CASE("Unit testing the Integrator class (using FieldIntegrator)", "[unit][field][Integrator][FieldIntegrator]"){

  int dimMax = 3;
  int orderMax = 5;
  double tol = 1e-2;

  std::string meshDirPath = TestUtils::getRessourcePath() + "/meshes/regression/";

  SECTION("Testing construcion of the integrators"){
    for(int iDim = 2; iDim < dimMax + 1; iDim++){
      for(int iOrd = 1; iOrd < orderMax + 1; iOrd++){
        Mesh myMesh(iDim, iOrd, "simplex");
        std::string meshname = "regression_dim-" + std::to_string(iDim) + "_h-2e-1_ord-" + std::to_string(iOrd) + ".h5";
        HDF5Io hdfio(&myMesh);
        hdfio.load(meshDirPath + meshname);
        const ReferenceElement * refEl = myMesh.getReferenceElement();
        Field myCellField(&myMesh, Cell, refEl->getNumNodes(), 1);
        Field myNodeField(&myMesh, Node, 1, 1);
        Field myFaceField(&myMesh, Face, refEl->getFaceElement()->getNumNodes(), 1);
        std::vector<double> dbuff(1, 0.0);
        CHECK_THROWS(FieldIntegrator(NULL, 1));
        CHECK_THROWS(FieldIntegrator(&myMesh, -1));
        CHECK_THROWS(FieldIntegrator(&myMesh, 0));
        CHECK_NOTHROW(FieldIntegrator(&myMesh, 1));
        FieldIntegrator myFI(&myMesh, 1);
        CHECK_THROWS(myFI.setField(NULL));
        Field myNoneField(&myMesh, None, 1, 1);
        CHECK_THROWS(myFI.setField(&myNoneField));
        CHECK_NOTHROW(myFI.setField(&myCellField));
        CHECK_NOTHROW(myFI.integrate(&dbuff));
        CHECK_NOTHROW(myFI.setField(&myFaceField));
        CHECK_NOTHROW(myFI.integrate(&dbuff));
        CHECK_NOTHROW(myFI.setField(&myNodeField));
        CHECK_NOTHROW(myFI.integrate(&dbuff));
      }
    }
  };

  SECTION("Testing constant value field integration"){
    for(int iDim = 2; iDim < dimMax + 1; iDim++){
      for(int iOrd = 1; iOrd < orderMax + 1; iOrd++){
        Mesh myMesh(iDim, iOrd, "simplex");
        std::string meshname = "regression_dim-" + std::to_string(iDim) + "_h-2e-1_ord-" + std::to_string(iOrd) + ".h5";
        HDF5Io hdfio(&myMesh);
        hdfio.load(meshDirPath + meshname);
        const ReferenceElement * refEl = myMesh.getReferenceElement();
        Field myCellField(&myMesh, Cell, refEl->getNumNodes(), 1);
        Field myNodeField(&myMesh, Node, 1, 1);
        std::vector<double> dbuff(1, 0.0);
        FieldIntegrator myFI(&myMesh, 1);
        std::fill(myCellField.getValues()->begin(), myCellField.getValues()->end(), 2.0);
        std::fill(myNodeField.getValues()->begin(), myNodeField.getValues()->end(), 5.75);
        //integrate cell field
        myFI.setField(&myCellField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - 2.0) < tol);
        myFI.setField(&myNodeField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - 5.75) < tol);
      }
    }
  };

  SECTION("Testing monomial value field integration"){
    for(int iDim = 2; iDim < dimMax + 1; iDim++){
      for(int iOrd = 1; iOrd < orderMax + 1; iOrd++){
        Mesh myMesh(iDim, iOrd, "simplex");
        std::string meshname = "regression_dim-" + std::to_string(iDim) + "_h-2e-1_ord-" + std::to_string(iOrd) + ".h5";
        HDF5Io hdfio(&myMesh);
        hdfio.load(meshDirPath + meshname);
        const ReferenceElement * refEl = myMesh.getReferenceElement();
        Field myCellField(&myMesh, Cell, refEl->getNumNodes(), 1);
        Field myNodeField(&myMesh, Node, 1, 1);
        std::vector<double> dbuff(1, 0.0);
        FieldIntegrator myFI(&myMesh, 1);
        std::vector<int> cell;
        std::vector<double> coords;
        for(int iEl = 0; iEl < myMesh.getNumberCells(); iEl++){
          myMesh.getCell(iEl, &cell);
          for(int iN = 0; iN < refEl->getNumNodes(); iN++){
            myMesh.getPoint(cell[iN], &coords);
            myCellField.getValues()->at(iEl * (refEl->getNumNodes()) + iN) = std::pow(coords[0], 2*iOrd);
            myNodeField.getValues()->at(cell[iN]) = std::pow(coords[0], 2*iOrd);
          }
        }
        //integrate cell field
        myFI.setField(&myCellField);
        myFI.integrate(&dbuff);
        double intVal = integralMonomial(2*iOrd);
        CHECK(std::abs(dbuff[0] - intVal) < tol);
        myFI.setField(&myNodeField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - intVal) < tol);
      }
    }
  };
};
