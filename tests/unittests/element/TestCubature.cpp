#include <catch2/catch.hpp>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>
#include "Cubature.h"

double x10(double x, double y, double z){
  return (std::pow(x, 10.0));
}

double testPoly1(double x, double y, double z){
  return (std::pow(x, 4.0) + std::pow(y, 3.0) + std::pow(z, 2.0) + (x*y*z));
}

double testPoly2(double x, double y, double z){
  return (std::pow(x*y, 3.0) + std::pow(y*z, 3.0) + std::pow(z*x, 3.0));
}

using namespace hfox;

TEST_CASE("Unittesting the Cubature class.", "[unit][cubature][element]"){
  elementGeometry s = simplex;
  elementGeometry o = orthotope;
  int dimMax = 3;
  int orderMax = 10;
  // dimension loop
  for(int i = 0; i < (dimMax+1); i++){
    // order loop
    for(int j = 0; j < (orderMax+1); j++){
      SECTION("Testing Construction (dim = "+std::to_string(i)+", order= "+std::to_string(j)+")"){
        CHECK_NOTHROW(Cubature(i, j, s));
        CHECK_NOTHROW(Cubature(i, j, o));
        CHECK_THROWS(Cubature(100, j, o));
        CHECK_THROWS(Cubature(i, 100, o));
      };
    }
  }

  std::vector<Cubature *> cubSimplexes((dimMax + 1)*(orderMax + 1));
  std::vector<Cubature *> cubOrthotopes((dimMax + 1)*(orderMax + 1));
  std::vector<int> numIPsSimplex((dimMax + 1)*(orderMax +1));
  std::vector<int> numIPsSimplexd2 = {1, 1, 3, 6, 6, 7, 12, 15, 16, 19, 25};
  std::vector<int> numIPsSimplexd3 = {1, 1, 4, 8, 14, 14, 24, 35, 46, 59, 81};
  std::vector<int> numIPsOrthotope((dimMax + 1)*(orderMax + 1));
  std::vector<int> numIPsOrthotoped2 = {1, 1, 4, 4, 8, 8, 12, 12, 20, 20, 28};
  std::vector<int> numIPsOrthotoped3 = {1, 1, 6, 6, 14, 14, 34, 34, 58, 58, 90};
  // dimension loop
  for(int i = 0; i < (dimMax+1); i++){
    // order loop
    for(int j = 0; j < (orderMax+1); j++){
      int index = (orderMax+1)*i + j;
      cubSimplexes[index] = new Cubature(i, j, s);
      cubOrthotopes[index] = new Cubature(i, j, o);
      if(i != 0){
        if(i == 1){
          numIPsSimplex[index] = std::ceil((j+1.0)/2.0);
          numIPsOrthotope[index] = std::ceil((j+1.0)/2.0);
        } else if(i == 2){
          numIPsSimplex[index] = numIPsSimplexd2[j];
          numIPsOrthotope[index] = numIPsOrthotoped2[j];
        } else if(i == 3){
          numIPsSimplex[index] = numIPsSimplexd3[j];
          numIPsOrthotope[index] = numIPsOrthotoped3[j];
        }
      } else{
        numIPsSimplex[index] = 0;
        numIPsOrthotope[index] = 0;
      }
    }
  }
  // dimension loop
  for(int i = 0; i < (dimMax +1); i++){
    // order loop
    for(int j = 0; j < (orderMax+1); j++){
      SECTION("Testing easy gets (dim = "+std::to_string(i)+", order= "+std::to_string(j)+")"){
        int index = (orderMax+1)*i + j;
        CHECK(cubSimplexes[index]->getDimension() == i);
        CHECK(cubOrthotopes[index]->getDimension() == i);
        CHECK(cubSimplexes[index]->getOrder() == j);
        CHECK(cubOrthotopes[index]->getOrder() == j);
        CHECK(cubSimplexes[index]->getGeometry() == s);
        CHECK(cubOrthotopes[index]->getGeometry() == o);
        CHECK(cubSimplexes[index]->getNumIPs() == numIPsSimplex[index]);
        CHECK(cubOrthotopes[index]->getNumIPs() == numIPsOrthotope[index]);
      };
    }
  }

  std::vector<double> simplexVolumes = {0.0, 2.0, 2.0, 4.0/3.0};
  std::vector<double> orthotopeVolumes = {0.0, 2.0, 4.0, 8.0};
  const std::vector< std::vector<double> > * IPs;
  const std::vector<double> * weights;
  double integral = 0.0;
  // dim loop
  for(int i = 1; i < (dimMax+1); i++){
    // order loop
    for(int j = 0; j < (orderMax+1); j++){
      SECTION("Testing volume integral (dim = "+std::to_string(i)+", order= "+std::to_string(j)+")"){
        int index = (orderMax + 1)*i + j;
        IPs = cubOrthotopes[index]->getIPCoords();
        weights = cubOrthotopes[index]->getIPWeights();
        for(int i = 0; i < cubOrthotopes[index]->getNumIPs(); i++){
          integral += (*weights)[i]*1.0;
        }
        CHECK(integral == Approx(orthotopeVolumes[i]));
        integral = 0.0;
        IPs = cubSimplexes[index]->getIPCoords();
        weights = cubSimplexes[index]->getIPWeights();
        for(int i = 0; i < cubSimplexes[index]->getNumIPs(); i++){
          integral += (*weights)[i]*1.0;
        }
        CHECK(integral == Approx(simplexVolumes[i]));
        integral = 0.0;
      };
    }
  }

  // x^10 in 1d
  double integralx10Orthotope1d = 2.0/11.0;
  // x^10 in 2d
  double integralx10Orthotope2d = 4.0/11.0;
  // x^10 in 3d
  double integralx10Orthotope3d = 8.0/11.0;
  //x^10 in simplexes
  double integralx10Simplex2d = 2.0/11.0;
  double integralx10Simplex3d = 24.0/143.0;;
  //testPoly1
  double integraltest1Simplex1d = 2.0/5.0;
  double integraltest1Orthotope2d = 4.0/5.0;
  double integraltest1Simplex2d = 0.0;
  double integraltest1Orthotope3d = 64.0/15.0;
  double integraltest1Simplex3d = 136.0/315.0;
  //testPoly2
  double integraltest2Simplex1d = 0.0;
  double integraltest2Orthotope2d = 0.0;
  double integraltest2Simplex2d = 0.0;
  double integraltest2Orthotope3d = 0.0;
  double integraltest2Simplex3d = 4.0/15.0;

  SECTION("Testing integration capabilities."){
    int cubIndex;
    const std::vector< std::vector<double> > * IPs;
    const std::vector<double> * weights;
    double integral = 0.0;
    //x10
    //1d
    cubIndex = 11 + 10;
    IPs = cubOrthotopes[cubIndex]->getIPCoords();
    weights = cubOrthotopes[cubIndex]->getIPWeights();
    for(int i = 0; i < cubOrthotopes[cubIndex]->getNumIPs(); i++){
      integral += (*weights)[i]*x10(((*IPs)[i])[0], 0.0, 0.0);
    }
    CHECK(integral == Approx(integralx10Orthotope1d));
    integral = 0.0;
    //2d
    cubIndex = 2*11 + 10;
    IPs = cubOrthotopes[cubIndex]->getIPCoords();
    weights = cubOrthotopes[cubIndex]->getIPWeights();
    for(int i = 0; i < cubOrthotopes[cubIndex]->getNumIPs(); i++){
      integral += (*weights)[i]*x10(((*IPs)[i])[0], ((*IPs)[i])[1], 0.0);
    }
    CHECK(integral == Approx(integralx10Orthotope2d));
    integral = 0.0;

    IPs = cubSimplexes[cubIndex]->getIPCoords();
    weights = cubSimplexes[cubIndex]->getIPWeights();
    for(int i = 0; i < cubSimplexes[cubIndex]->getNumIPs(); i++){
      integral += (*weights)[i]*x10(((*IPs)[i])[0], ((*IPs)[i])[1], 0.0);
    }
    CHECK(integral == Approx(integralx10Simplex2d));
    integral = 0.0;
    //3d
    cubIndex = 3*11 + 10;
    IPs = cubOrthotopes[cubIndex]->getIPCoords();
    weights = cubOrthotopes[cubIndex]->getIPWeights();
    for(int i = 0; i < cubOrthotopes[cubIndex]->getNumIPs(); i++){
      integral += (*weights)[i]*x10(((*IPs)[i])[0], ((*IPs)[i])[1], ((*IPs)[i])[2]);
    }
    CHECK(integral == Approx(integralx10Orthotope3d));
    integral = 0.0;

    IPs = cubSimplexes[cubIndex]->getIPCoords();
    weights = cubSimplexes[cubIndex]->getIPWeights();
    for(int i = 0; i < cubSimplexes[cubIndex]->getNumIPs(); i++){
      integral += (*weights)[i]*x10(((*IPs)[i])[0], ((*IPs)[i])[1], ((*IPs)[i])[2]);
    }
    CHECK(integral == Approx(integralx10Simplex3d));
    integral = 0.0;
    //testpoly1
    //1d
    cubIndex = 11 + 4;
    IPs = cubSimplexes[cubIndex]->getIPCoords();
    weights = cubSimplexes[cubIndex]->getIPWeights();
    for(int i = 0; i < cubSimplexes[cubIndex]->getNumIPs(); i++){
      integral += (*weights)[i]*testPoly1(((*IPs)[i])[0], 0.0, 0.0);
    }
    CHECK(integral == Approx(integraltest1Simplex1d));
    integral = 0.0;
    //2d
    cubIndex = 2*11 + 4;
    IPs = cubOrthotopes[cubIndex]->getIPCoords();
    weights = cubOrthotopes[cubIndex]->getIPWeights();
    for(int i = 0; i < cubOrthotopes[cubIndex]->getNumIPs(); i++){
      integral += (*weights)[i]*testPoly1(((*IPs)[i])[0], ((*IPs)[i])[1], 0.0);
    }
    CHECK(integral == Approx(integraltest1Orthotope2d));
    integral = 0.0;

    IPs = cubSimplexes[cubIndex]->getIPCoords();
    weights = cubSimplexes[cubIndex]->getIPWeights();
    for(int i = 0; i < cubSimplexes[cubIndex]->getNumIPs(); i++){
      integral += (*weights)[i]*testPoly1(((*IPs)[i])[0], ((*IPs)[i])[1], 0.0);
    }
    CHECK(integral == Approx(integraltest1Simplex2d).margin(1e-12));
    integral = 0.0;
    //3d
    cubIndex = 3*11 + 4;
    IPs = cubOrthotopes[cubIndex]->getIPCoords();
    weights = cubOrthotopes[cubIndex]->getIPWeights();
    for(int i = 0; i < cubOrthotopes[cubIndex]->getNumIPs(); i++){
      integral += (*weights)[i]*testPoly1(((*IPs)[i])[0], ((*IPs)[i])[1], ((*IPs)[i])[2]);
    }
    CHECK(integral == Approx(integraltest1Orthotope3d));
    integral = 0.0;

    IPs = cubSimplexes[cubIndex]->getIPCoords();
    weights = cubSimplexes[cubIndex]->getIPWeights();
    for(int i = 0; i < cubSimplexes[cubIndex]->getNumIPs(); i++){
      integral += (*weights)[i]*testPoly1(((*IPs)[i])[0], ((*IPs)[i])[1], ((*IPs)[i])[2]);
    }
    CHECK(integral == Approx(integraltest1Simplex3d));
    integral = 0.0;
   //testpoly2
    //1d
    cubIndex = 11 + 6;
    IPs = cubSimplexes[cubIndex]->getIPCoords();
    weights = cubSimplexes[cubIndex]->getIPWeights();
    for(int i = 0; i < cubSimplexes[cubIndex]->getNumIPs(); i++){
      integral += (*weights)[i]*testPoly2(((*IPs)[i])[0], 0.0, 0.0);
    }
    CHECK(integral == Approx(integraltest2Simplex1d));
    integral = 0.0;
    //2d
    cubIndex = 2*11 + 6;
    IPs = cubOrthotopes[cubIndex]->getIPCoords();
    weights = cubOrthotopes[cubIndex]->getIPWeights();
    for(int i = 0; i < cubOrthotopes[cubIndex]->getNumIPs(); i++){
      integral += (*weights)[i]*testPoly2(((*IPs)[i])[0], ((*IPs)[i])[1], 0.0);
    }
    CHECK(integral == Approx(integraltest2Orthotope2d));
    integral = 0.0;

    IPs = cubSimplexes[cubIndex]->getIPCoords();
    weights = cubSimplexes[cubIndex]->getIPWeights();
    for(int i = 0; i < cubSimplexes[cubIndex]->getNumIPs(); i++){
      integral += (*weights)[i]*testPoly2(((*IPs)[i])[0], ((*IPs)[i])[1], 0.0);
    }
    CHECK(integral == Approx(integraltest2Simplex2d).margin(1e-12));
    integral = 0.0;
    //3d
    cubIndex = 3*11 + 6;
    IPs = cubOrthotopes[cubIndex]->getIPCoords();
    weights = cubOrthotopes[cubIndex]->getIPWeights();
    for(int i = 0; i < cubOrthotopes[cubIndex]->getNumIPs(); i++){
      integral += (*weights)[i]*testPoly2(((*IPs)[i])[0], ((*IPs)[i])[1], ((*IPs)[i])[2]);
    }
    CHECK(integral == Approx(integraltest2Orthotope3d).margin(1e-12));
    integral = 0.0;

    IPs = cubSimplexes[cubIndex]->getIPCoords();
    weights = cubSimplexes[cubIndex]->getIPWeights();
    for(int i = 0; i < cubSimplexes[cubIndex]->getNumIPs(); i++){
      integral += (*weights)[i]*testPoly2(((*IPs)[i])[0], ((*IPs)[i])[1], ((*IPs)[i])[2]);
    }
    CHECK(integral == Approx(integraltest2Simplex3d));
    integral = 0.0;
  };
  // dim loop
  int index;
  for(int i = 0; i < (dimMax+1); i++){
    // order loop
    for(int j = 0; j < (orderMax+1); j++){
      index = (orderMax+1)*i + j;
      delete cubSimplexes[index];
      delete cubOrthotopes[index];
    }
  }
};
