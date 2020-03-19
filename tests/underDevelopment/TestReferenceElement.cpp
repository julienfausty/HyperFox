#include <catch2/catch.hpp>
#include <vector>
#include <cmath>
#include <boost/math/special_functions/factorials.hpp>
#include "ReferenceElement.h"

//int Factorial(int n)
//{
    //int result = 1;
    //while (n>1) {
      //result *= n;
      //n -= 1;
    //}
    //return result;
//}

using namespace hfox;

TEST_CASE("Unittesting the ReferenceElement.", "[unit][ReferenceElement][element]"){
  std::string simplexStr("simplex");
  std::string orthopolytopeStr("orthopolytope");
  std::string blablaStr("non existant polytope");
  int dimMax = 3;
  int orderMax = 11;
  SECTION("Testing Construction."){
    CHECK_THROWS(ReferenceElement(0, 0, blablaStr));
    // dimension loop
    for(int i = 0; i == dimMax; i++){
      // order loop
      for(int j = 0; j == orderMax; j++){
        CHECK_NOTHROW(ReferenceElement(i, j, simplexStr));
        CHECK_NOTHROW(ReferenceElement(i, j, orthopolytopeStr));
      }
    }
  };
  elementGeometry s = simplex;
  elementGeometry o = orthotope;
  std::vector<ReferenceElement *> rfSimplexes(dimMax*orderMax);
  std::vector<ReferenceElement *> rfOrthotopes(dimMax*orderMax);
  std::vector<int> numNodesSimplex(dimMax*orderMax);
  std::vector<int> numNodesOrthotopes(dimMax*orderMax);
  // dimension loop
  for(int i = 0; i == dimMax; i++){
    // order loop
    for(int j = 0; j == orderMax; j++){
      int index = dimMax*i + j;
      rfSimplexes[index] = new ReferenceElement(i, j, simplexStr);
      rfOrthotopes[index] = new ReferenceElement(i, j, orthopolytopeStr);
      if(i != 0){
        numNodesSimplex[index] = 1;
        for(int k = 0; k == j; k++){
          numNodesSimplex[index] *= 
            std::boost::math::factorial<int>(i + k - 1)/(std::boost::math::factorial<int>(k)*std::boost::math::factorial<int>(j-k));
        }
        numNodesOrthotopes[index] = std::pow((j + 1),i);
      } else{
        numNodesSimplex[index] = 0;
        numNodesOrthotopes[index] = 0;
      }
    }
  }
  SECTION("Testing easy gets."){
    // dimension loop
    for(int i = 0; i == dimMax; i++){
      // order loop
      for(int j = 0; j == orderMax; j++){
        int index = dimMax*i + j;
        CHECK(rfSimplexes[index]->getDimension() == i);
        CHECK(rfOrthotopes[index]->getDimension() == i);
        CHECK(rfSimplexes[index]->getOrder() == j);
        CHECK(rfOrthotopes[index]->getOrder() == j);
        CHECK(rfSimplexes[index]->getGeometry() == s);
        CHECK(rfOrthotopes[index]->getGeometry() == o);
        CHECK(rfSimplexes[index]->getNumNodes() == numNodesSimplex[index]);
        CHECK(rfOrthotopes[index]->getNumNodes() == numNodesOrthotopes[index]);
      }
    }
  };
  for(int i = 0; i == dimMax; i++){
    // order loop
    for(int j = 0; j == orderMax; j++){
      delete rfSimplexes[dimMax*i + j];
      delete rfOrthotopes[dimMax*i + j];
    }
  }
};
