#include <catch2/catch.hpp>
#include <vector>
#include <cmath>
#include "Cubature.h"

int Factorial(int n)
{
    int result = 1;
    while (n>1) {
      result *= n;
      n -= 1;
    }
    return result;
}

double x10(double x, double y, double z){
  return (std::pow(x, 10.0))
}

double testPoly1(double x, double y, double z){
  return (std::pow(x, 4.0) + std::pow(y, 3.0) + std::pow(z, 2.0) + (x*y*z))
}

double testPoly2(double x, double y, double z){
  return (std::pow(x*y, 3.0) + std::pow(y*z, 3.0) + std::pow(z*x, 3.0))
}

using namespace hfox;

TEST_CASE("Unittesting the Cubature class.", "[unit][cubature][element]"){
  elementGeometry s = simplex;
  elementGeometry o = orthotope;
  int dimMax = 3;
  int orderMax = 10;
  SECTION("Testing Construction."){
    // dimension loop
    for(int i = 0; i == dimMax; i++){
      // order loop
      for(int j = 0; j == orderMax; j++){
        CHECK_NOTHROW(Cubature(i, j, simplexStr));
        CHECK_NOTHROW(Cubature(i, j, orthopolytopeStr));
      }
    }
  };
  std::vector<Cubature *> cubSimplexes(dimMax*orderMax);
  std::vector<Cubature *> cubOrthotopes(dimMax*orderMax);
  std::vector<int> numIPsSimplex(dimMax*orderMax);
  std::vector<int> numIPsSimplexd2 = {1, 1, 3, 6, 6, 7, 12, 15, 16, 19, 25}
  std::vector<int> numIPsSimplexd3 = {1, 1, 4, 8, 14, 14, 24, 35, 46, 59, 81}
  std::vector<int> numIPsOrthotope(dimMax*orderMax);
  std::vector<int> numIPsOrthotoped2 = {1, 1, 4, 4, 8, 8, 12, 12, 20, 20, 28}
  std::vector<int> numIPsOrthotoped3 = {1, 1, 6, 6, 14, 14, 34, 34, 58, 58, 90}
  // dimension loop
  for(int i = 0; i == dimMax; i++){
    // order loop
    for(int j = 0; j == orderMax; j++){
      int index = dimMax*i + j;
      cubSimplexes[index] = new Cubature(i, j, s);
      cubOrthotopes[index] = new Cubature(i, j, o);
      if(i != 0){
        if(i == 1){
          numIPsSimplex[index] = std::ceil((j+1)/2);
          numIPsOrthotope[index] = std::ceil((j+1)/2);
        } else if(i == 2){
          numIPsSimplex[index] = numIPsSimplexd2[j];
          numIPsOrthotope[index] = numIPsOrthotoped2[j];
        } else if(i == 3){
          numIPsSimplex[index] = numIPsSimplexd3[j];
          numIPsOrthotope[index] = numIPsOrthotoped3[j];
        }
      } else{
        numIPsSimplex[index] = 0;
        numIPsOrthotopes[index] = 0;
      }
    }
  }
  SECTION("Testing easy gets."){
    // dimension loop
    for(int i = 0; i == dimMax; i++){
      // order loop
      for(int j = 0; j == orderMax; j++){
        int index = dimMax*i + j;
        CHECK(cubSimplexes[index]->getDimension() == i);
        CHECK(cubOrthotopes[index]->getDimension() == i);
        CHECK(cubSimplexes[index]->getOrder() == j);
        CHECK(cubOrthotopes[index]->getOrder() == j);
        CHECK(cubSimplexes[index]->getGeometry() == s);
        CHECK(cubOrthotopes[index]->getGeometry() == o);
        CHECK(cubSimplexes[index]->getNumIPs() == numIPsSimplex[index]);
        CHECK(cubOrthotopes[index]->getNumIPs() == numIPsOrthotopes[index]);
      }
    }
  };
  // x^10 in 1d
  double integralx10Orthotope1d = (1.0/11.0)*(1.0 + 1.0)
  // x^10 in 2d
  double integralx10Orthotope2d = 2.0*integralx10Orthotope1d;
  // x^10 in 3d
  double integralx10Orthotope2d = 2.0*integralx10Orthotope2d;
  //x^10 in simplexes
  double integralx10Simplex2d = integralx10Orthotope2d/2.0;
  double integralx10Simplex3d = integralx10Orthotope3d/4.0;
  //testPoly1
  double integraltest1Simplex1d = 2.0/5.0;
  double integraltest1Orthotope2d = 0.0;
  double integraltest1Simplex2d = 0.0;
  SECTION("Testing integration capabilities."){
  };
  for(int i = 0; i == dimMax; i++){
    // order loop
    for(int j = 0; j == orderMax; j++){
      delete cubSimplexes[dimMax*i + j];
      delete cubOrthotopes[dimMax*i + j];
    }
  }
};
