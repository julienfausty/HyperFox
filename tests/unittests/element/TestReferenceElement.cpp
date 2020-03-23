#include <catch2/catch.hpp>
#include <vector>
#include <cmath>
#include <string>
#include <boost/math/special_functions/factorials.hpp>
#include "ReferenceElement.h"

using namespace hfox;

TEST_CASE("Unittesting the ReferenceElement.", "[unit][ReferenceElement][element]"){
  std::string simplexStr("simplex");
  std::string orthopolytopeStr("orthopolytope");
  std::string blablaStr("non existant polytope");
  int dimMax = 3;
  int orderMax = 10;
  // dimension loop
  for(int i = 0; i < (dimMax+1); i++){
    // order loop
    for(int j = 0; j < (orderMax+1); j++){
      SECTION("Testing Construction (dim = "+std::to_string(i)+", order = "+std::to_string(j)+")"){
        CHECK_THROWS(ReferenceElement(0, 0, blablaStr));
        CHECK_NOTHROW(ReferenceElement(i, j, simplexStr));
        CHECK_NOTHROW(ReferenceElement(i, j, orthopolytopeStr));
        CHECK_THROWS(ReferenceElement(100, j, simplexStr));
        CHECK_THROWS(ReferenceElement(i, 100, simplexStr));
      };
    }
  }
  elementGeometry s = simplex;
  elementGeometry o = orthotope;
  std::vector<ReferenceElement *> rfSimplexes((dimMax+1)*(orderMax + 1));
  std::vector<ReferenceElement *> rfOrthotopes((dimMax + 1)*(orderMax +1));
  std::vector<int> numNodesSimplex((dimMax + 1)*(orderMax + 1));
  std::vector<int> numNodesOrthotopes((dimMax + 1)*(orderMax + 1));
  // dimension loop
  for(int i = 0; i < (dimMax+1); i++){
    // order loop
    for(int j = 0; j < (orderMax+1); j++){
      int index = dimMax*i + j;
      rfSimplexes[index] = new ReferenceElement(i, j, simplexStr);
      rfOrthotopes[index] = new ReferenceElement(i, j, orthopolytopeStr);
      if(i != 0){
        numNodesSimplex[index] = 1;
        for(int k = 1; k < (j + 1); k++){
          numNodesSimplex[index] += 
            boost::math::factorial<double>(i + k - 1)/(boost::math::factorial<double>(k)*boost::math::factorial<double>(i-1));
        }
        numNodesOrthotopes[index] = std::pow((j + 1),i);
      } else{
        numNodesSimplex[index] = 0;
        numNodesOrthotopes[index] = 0;
      }
    }
  }
  std::vector<double> numFacesSimplex = {0, 2, 3, 4};
  std::vector<double> numFacesOrthotope = {0, 2, 4, 6};
  // dimension loop
  for(int i = 0; i < (dimMax + 1); i++){
    // order loop
    for(int j = 0; j < (orderMax + 1); j++){
      SECTION("Testing easy gets (dim = "+std::to_string(i)+", order = "+std::to_string(j)+")"){
        int index = orderMax*i + j;
        CHECK(rfSimplexes[index]->getDimension() == i);
        CHECK(rfOrthotopes[index]->getDimension() == i);
        CHECK(rfSimplexes[index]->getOrder() == j);
        CHECK(rfOrthotopes[index]->getOrder() == j);
        CHECK(rfSimplexes[index]->getGeometry() == s);
        CHECK(rfOrthotopes[index]->getGeometry() == o);
        CHECK(rfSimplexes[index]->getNumNodes() == numNodesSimplex[index]);
        CHECK(rfOrthotopes[index]->getNumNodes() == numNodesOrthotopes[index]);
        CHECK(rfSimplexes[index]->getNodes()->size() == numNodesSimplex[index]);
        CHECK(rfOrthotopes[index]->getNodes()->size() == numNodesOrthotopes[index]);
        if(i != 0){
          CHECK(rfSimplexes[index]->getFaceElement()->getNumNodes() == numNodesSimplex[orderMax*(i-1) + j]);
          CHECK(rfOrthotopes[index]->getFaceElement()->getNumNodes() == numNodesOrthotopes[orderMax*(i-1)+j]);
        } else{
          CHECK(rfSimplexes[index]->getFaceElement() == NULL);
          CHECK(rfOrthotopes[index]->getFaceElement() == NULL);
        }
        CHECK(rfSimplexes[index]->getNumFaces() == numFacesSimplex[i]);
        CHECK(rfOrthotopes[index]->getNumFaces() == numFacesOrthotope[i]);
        CHECK(rfSimplexes[index]->getFaceNodes()->size() == numFacesSimplex[i]);
        CHECK(rfOrthotopes[index]->getFaceNodes()->size() == numFacesOrthotope[i]);
      };
    }
  }
  // dimension loop
  for(int i = 0; i < (dimMax + 1); i++){
    // order loop
    for(int j = 0; j < (orderMax + 1); j++){
      SECTION("Lagrange property (dim = "+std::to_string(i)+", order = "+std::to_string(j)+")"){
        int index = orderMax*i + j;
        const std::vector< std::vector<double> > * nodes;
        std::vector<double> res;
        nodes = rfSimplexes[index]->getNodes();
        for(int k = 0; k < nodes->size(); k++){
          res = rfSimplexes[index]->interpolate((*nodes)[k]);
          for(int l = 0; l < res.size(); l++){
            if(l!=k){
              CHECK(res[l] == Approx(0.0).margin(1e-12));
            }else{
              CHECK(res[l] == Approx(1.0));
            }
          }
        }
        nodes = rfOrthotopes[index]->getNodes();
        for(int k = 0; k < nodes->size(); k++){
          res = rfOrthotopes[index]->interpolate((*nodes)[k]);
          for(int l = 0; l < res.size(); l++){
            if(l!=k){
              CHECK(res[l] == Approx(0.0).margin(1e-12));
            }else{
              CHECK(res[l] == Approx(1.0));
            }
          }
        }
      };
      SECTION("Monomial interpolation (dim = "+std::to_string(i)+", order = "+std::to_string(j)+")"){
        // x^j + y^j + z^j
        int index = orderMax*i + j;
        const std::vector< std::vector<double> > * nodes;
        std::vector<double> res, contributions;
        std::vector<double> interpPoint(i, -0.5);
        double val, interpVal;
        val = 0;
        for(int d = 0; d < i; d++){
          val += std::pow(interpPoint[d], j);
        }
        nodes = rfSimplexes[index]->getNodes();
        res.resize(nodes->size());
        for(int k = 0; k < nodes->size(); k++){
          res[k] = 0;
          for(int d = 0; d < i; d++){
            res[k] += std::pow((*nodes)[k][d], j);
          }
        }
        interpVal = 0;
        contributions = rfSimplexes[index]->interpolate(interpPoint);
        for(int k = 0; k < contributions.size(); k++){
          interpVal += contributions[k]*res[k];
        }
        CHECK(interpVal == Approx(val));


        nodes = rfOrthotopes[index]->getNodes();
        res.resize(nodes->size());
        for(int k = 0; k < nodes->size(); k++){
          res[k] = 0;
          for(int d = 0; d < i; d++){
            res[k] += std::pow((*nodes)[k][d], j);
          }
        }
        interpVal = 0;
        contributions = rfOrthotopes[index]->interpolate(interpPoint);
        for(int k = 0; k < contributions.size(); k++){
          interpVal += contributions[k]*res[k];
        }
        CHECK(interpVal == Approx(val));


        // {jx^(j-1), jy^(j-1), jz^(j-1)}
        std::vector<double> derivVal, interpDerivVal;
        std::vector< std::vector<double> > derivContributions;
        derivVal.resize(i);
        interpDerivVal.resize(i);
        for(int d = 0; d < i; d++){
          derivVal[d] = j*std::pow(interpPoint[d], j-1);
        }
        nodes = rfSimplexes[index]->getNodes();
        res.resize(nodes->size());
        for(int k = 0; k < nodes->size(); k++){
          res[k] = 0;
          for(int d = 0; d < i; d++){
            res[k] += std::pow((*nodes)[k][d], j);
          }
        }
        derivContributions = rfSimplexes[index]->interpolateDeriv(interpPoint, 1);
        for(int d = 0; d < i; d++){
          interpDerivVal[d] = 0.0;
          for(int k = 0; k < contributions.size(); k++){
            interpDerivVal[d] += derivContributions[k][d]*res[k];
          }
        }


        nodes = rfOrthotopes[index]->getNodes();
        res.resize(nodes->size());
        for(int k = 0; k < nodes->size(); k++){
          res[k] = 0;
          for(int d = 0; d < i; d++){
            res[k] += std::pow((*nodes)[k][d], j);
          }
        }
        derivContributions = rfOrthotopes[index]->interpolateDeriv(interpPoint, 1);
        for(int d = 0; d < i; d++){
          interpDerivVal[d] = 0.0;
          for(int k = 0; k < contributions.size(); k++){
            interpDerivVal[d] += derivContributions[k][d]*res[k];
          }
        }
        for(int d = 0; d < i; d++){
          CHECK(interpDerivVal[d] == Approx(derivVal[d]));
        }
      };
    }
  }
  // dimension loop
  for(int i = 0; i < (dimMax + 1); i++){
    // order loop
    for(int j = 0; j < (orderMax + 1); j++){
      delete rfSimplexes[orderMax*i + j];
      delete rfOrthotopes[orderMax*i + j];
    }
  }
};
