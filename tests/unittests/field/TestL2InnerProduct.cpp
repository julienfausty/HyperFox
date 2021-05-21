#include <catch2/catch.hpp>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include "L2InnerProduct.h"
#include "Mesh.h"
#include "HDF5Io.h"
#include "TestUtils.h"

using namespace hfox;

namespace l2innerprod{

double integralMonomial(int power){
  return 1.0/(power + 1);
};

void monomial(std::vector<double> & coords, std::vector<double> * val, int power){
  val->resize(1, 0.0);
  val->at(0) = std::pow(coords[0], power);
};

}

TEST_CASE("Unit testing the L2InnerProduct class", "[unit][field][Integrator][L2InnerProduct]"){

  int dimMax = 3;
  int orderMax = 5;
  double tol = 1e-6;

  std::string meshDirPath = TestUtils::getRessourcePath() + "/meshes/regression/";

  SECTION("Testing construcion of the L2InnerProduct"){
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
        CHECK_THROWS(L2InnerProduct(NULL, 1));
        CHECK_THROWS(L2InnerProduct(&myMesh, -1));
        CHECK_THROWS(L2InnerProduct(&myMesh, 0));
        CHECK_NOTHROW(L2InnerProduct(&myMesh, 1));
        L2InnerProduct myFI(&myMesh, 1);
        CHECK_THROWS(myFI.setProduct(NULL, NULL));
        Field myNoneField(&myMesh, None, 1, 1);
        CHECK_THROWS(myFI.setProduct(&myNoneField, &myCellField));
        CHECK_THROWS(myFI.setProduct(&myNoneField, [](std::vector<double> & coord, std::vector<double> * val){l2innerprod::monomial(coord, val, 1);}));
        CHECK_THROWS(myFI.setProduct(&myCellField, NULL));
        CHECK_THROWS(myFI.setProduct(&myCellField, &myFaceField));
        //FieldField
        CHECK_NOTHROW(myFI.setProduct(&myCellField, &myNodeField));
        CHECK_NOTHROW(myFI.integrate(&dbuff));
        CHECK_NOTHROW(myFI.setProduct(&myFaceField, &myFaceField));
        CHECK_NOTHROW(myFI.integrate(&dbuff));
        CHECK_NOTHROW(myFI.setProduct(&myNodeField, &myCellField));
        CHECK_NOTHROW(myFI.integrate(&dbuff));
        //FieldFunction
        CHECK_NOTHROW(myFI.setProduct(&myNodeField, [](std::vector<double> & coord, std::vector<double> * val){l2innerprod::monomial(coord, val, 1);}));
        CHECK_NOTHROW(myFI.integrate(&dbuff));
        CHECK_NOTHROW(myFI.setProduct([](std::vector<double> & coord, std::vector<double> * val){l2innerprod::monomial(coord, val, 1);}, &myCellField));
        CHECK_NOTHROW(myFI.integrate(&dbuff));
        //FunctionFunction
        CHECK_NOTHROW(myFI.setProduct([](std::vector<double> & coord, std::vector<double> * val){l2innerprod::monomial(coord, val, 1);}, [](std::vector<double> & coord, std::vector<double> * val){l2innerprod::monomial(coord, val, 1);}));
        CHECK_NOTHROW(myFI.integrate(&dbuff));
      }
    }
  };

  SECTION("Testing L2InnerProduct for constant members"){
    for(int iDim = 2; iDim < dimMax + 1; iDim++){
      for(int iOrd = 1; iOrd < orderMax + 1; iOrd++){
        Mesh myMesh(iDim, iOrd, "simplex");
        std::string meshname = "regression_dim-" + std::to_string(iDim) + "_h-2e-1_ord-" + std::to_string(iOrd) + ".h5";
        HDF5Io hdfio(&myMesh);
        hdfio.load(meshDirPath + meshname);
        const ReferenceElement * refEl = myMesh.getReferenceElement();
        Field myCellField(&myMesh, Cell, refEl->getNumNodes(), 2);
        Field myNodeField(&myMesh, Node, 1, 2);
        std::fill(myCellField.getValues()->begin(), myCellField.getValues()->end(), 2.0);
        std::fill(myNodeField.getValues()->begin(), myNodeField.getValues()->end(), 3.0);
        std::vector<double> dbuff(1, 0.0);
        L2InnerProduct myFI(&myMesh, 1);
        //FieldField
        myFI.setProduct(&myCellField, &myCellField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - 8.0) < tol);
        myFI.setProduct(&myNodeField, &myNodeField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - 18.0) < tol);
        myFI.setProduct(&myNodeField, &myCellField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - 12.0) < tol);
        myFI.setProduct(&myCellField, &myNodeField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - 12.0) < tol);
        //FieldFunction
        myFI.setProduct(&myCellField, [](std::vector<double> & coord, std::vector<double> * val){val->resize(2);std::fill(val->begin(), val->end(), 4.0);});
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - 16.0) < tol);
        myFI.setProduct([](std::vector<double> & coord, std::vector<double> * val){val->resize(2);std::fill(val->begin(), val->end(), 4.0);}, &myNodeField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - 24.0) < tol);
        //FunctionFunction
        myFI.setProduct([](std::vector<double> & coord, std::vector<double> * val){val->resize(2);std::fill(val->begin(), val->end(), 4.0);}, [](std::vector<double> & coord, std::vector<double> * val){val->resize(2);std::fill(val->begin(), val->end(), 4.0);});
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - 32.0) < tol);
      }
    }
  };

  SECTION("Testing L2InnerProduct for monomial scalars"){
    for(int iDim = 2; iDim < dimMax + 1; iDim++){
      for(int iOrd = 1; iOrd < orderMax + 1; iOrd++){
        Mesh myMesh(iDim, iOrd, "simplex");
        std::string meshname = "regression_dim-" + std::to_string(iDim) + "_h-2e-1_ord-" + std::to_string(iOrd) + ".h5";
        HDF5Io hdfio(&myMesh);
        hdfio.load(meshDirPath + meshname);
        const ReferenceElement * refEl = myMesh.getReferenceElement();
        Field myCellField(&myMesh, Cell, refEl->getNumNodes(), 1);
        Field myNodeField(&myMesh, Node, 1, 1);
        std::vector<int> cell;
        std::vector<double> coords;
        for(int iEl = 0; iEl < myMesh.getNumberCells(); iEl++){
          myMesh.getCell(iEl, &cell);
          for(int iN = 0; iN < refEl->getNumNodes(); iN++){
            myMesh.getPoint(cell[iN], &coords);
            myCellField.getValues()->at(iEl * (refEl->getNumNodes()) + iN) = std::pow(coords[0], iOrd);
            myNodeField.getValues()->at(cell[iN]) = std::pow(coords[0], iOrd);
          }
        }
        std::vector<double> dbuff(1, 0.0);
        L2InnerProduct myFI(&myMesh, 1);
        double res = l2innerprod::integralMonomial(2*iOrd);
        //FieldField
        myFI.setProduct(&myCellField, &myCellField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - res) < tol);
        myFI.setProduct(&myNodeField, &myNodeField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - res) < tol);
        myFI.setProduct(&myNodeField, &myCellField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - res) < tol);
        myFI.setProduct(&myCellField, &myNodeField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - res) < tol);
        //FieldFunction
        myFI.setProduct(&myCellField, [iOrd](std::vector<double> & coord, std::vector<double> * val){l2innerprod::monomial(coord, val, iOrd);});
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - res) < tol);
        myFI.setProduct([iOrd](std::vector<double> & coord, std::vector<double> * val){l2innerprod::monomial(coord, val, iOrd);}, &myNodeField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - res) < tol);
        //FunctionFunction
        myFI.setProduct([iOrd](std::vector<double> & coord, std::vector<double> * val){l2innerprod::monomial(coord, val, iOrd);}, [iOrd](std::vector<double> & coord, std::vector<double> * val){l2innerprod::monomial(coord, val, iOrd);});
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - res) < tol);
      }
    }
  };

  SECTION("Testing L2InnerProduct for monomial scalars"){
    for(int iDim = 2; iDim < dimMax + 1; iDim++){
      for(int iOrd = 1; iOrd < orderMax + 1; iOrd++){
        Mesh myMesh(iDim, iOrd, "simplex");
        std::string meshname = "regression_dim-" + std::to_string(iDim) + "_h-2e-1_ord-" + std::to_string(iOrd) + ".h5";
        HDF5Io hdfio(&myMesh);
        hdfio.load(meshDirPath + meshname);
        const ReferenceElement * refEl = myMesh.getReferenceElement();
        Field myCellField(&myMesh, Cell, refEl->getNumNodes(), 2);
        Field myNodeField(&myMesh, Node, 1, 2);
        std::vector<int> cell;
        std::vector<double> coords;
        for(int iEl = 0; iEl < myMesh.getNumberCells(); iEl++){
          myMesh.getCell(iEl, &cell);
          for(int iN = 0; iN < refEl->getNumNodes(); iN++){
            myMesh.getPoint(cell[iN], &coords);
            myCellField.getValues()->at((iEl * (refEl->getNumNodes()) + iN)*2) = std::pow(coords[0], iOrd);
            myCellField.getValues()->at((iEl * (refEl->getNumNodes()) + iN)*2 + 1) = std::pow(coords[0], iOrd);
            myNodeField.getValues()->at(cell[iN]*2) = std::pow(coords[0], iOrd);
            myNodeField.getValues()->at(cell[iN]*2 + 1) = std::pow(coords[0], iOrd);
          }
        }
        std::vector<double> dbuff(1, 0.0);
        L2InnerProduct myFI(&myMesh, 1);
        double res = 2.0*l2innerprod::integralMonomial(2*iOrd);
        //FieldField
        myFI.setProduct(&myCellField, &myCellField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - res) < tol);
        myFI.setProduct(&myNodeField, &myNodeField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - res) < tol);
        myFI.setProduct(&myNodeField, &myCellField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - res) < tol);
        myFI.setProduct(&myCellField, &myNodeField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - res) < tol);
        //FieldFunction
        myFI.setProduct(&myCellField, [iOrd](std::vector<double> & coord, std::vector<double> * val){l2innerprod::monomial(coord, val, iOrd);val->resize(2);val->at(1) = val->at(0);});
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - res) < tol);
        myFI.setProduct([iOrd](std::vector<double> & coord, std::vector<double> * val){l2innerprod::monomial(coord, val, iOrd);val->resize(2);val->at(1) = val->at(0);}, &myNodeField);
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - res) < tol);
        //FunctionFunction
        myFI.setProduct([iOrd](std::vector<double> & coord, std::vector<double> * val){l2innerprod::monomial(coord, val, iOrd);val->resize(2);val->at(1) = val->at(0);}, [iOrd](std::vector<double> & coord, std::vector<double> * val){l2innerprod::monomial(coord, val, iOrd);val->resize(2);val->at(1) = val->at(0);});
        myFI.integrate(&dbuff);
        CHECK(std::abs(dbuff[0] - res) < tol);
      }
    }
  };

};
