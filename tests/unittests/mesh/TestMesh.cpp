#include <catch2/catch.hpp>
#include <vector>
#include <algorithm>
#include "Mesh.h"

using namespace hfox;

bool arrayEqual4(const std::array<int, 4> & a, const std::array<int, 4> & b){
  bool areEqual = 1;
  for(int i = 0; i < 4; i++){
    areEqual = (a[i] == b[i]);
    if(!areEqual){
      break;
    }
  }
  return areEqual;
}

bool arrayEqual2(const std::array<int, 2> & a, const std::array<int, 2> & b){
  bool areEqual = 1;
  for(int i = 0; i < 2; i++){
    areEqual = (a[i] == b[i]);
    if(!areEqual){
      break;
    }
  }
  return areEqual;
}

TEST_CASE("Unit testing the Mesh class", "[unit][mesh][Mesh]"){
  //setup a 2D test mesh
  std::vector< std::vector<double> > testNodes;
  std::vector< std::vector<int> > testSimplexConnectivity;
  std::vector< std::vector<int> > testOrthotopeConnectivity;
  std::vector<double> bufferNode(2);
  bufferNode[0] = 0.0; bufferNode[1] = 0.0;
  testNodes.push_back(bufferNode);
  bufferNode[0] = 1.0;
  testNodes.push_back(bufferNode);
  bufferNode[1] = 1.0;
  testNodes.push_back(bufferNode);
  bufferNode[0] = 0.0;
  testNodes.push_back(bufferNode);
  bufferNode[0] = 0.5; bufferNode[1] = 0.5;
  testNodes.push_back(bufferNode);
  bufferNode[0] = 0.5; bufferNode[1] = 0.0;
  testNodes.push_back(bufferNode);
  bufferNode[0] = 1.0; bufferNode[1] = 0.5;
  testNodes.push_back(bufferNode);
  bufferNode[0] = 0.5; bufferNode[1] = 1.0;
  testNodes.push_back(bufferNode);
  bufferNode[0] = 0.0; bufferNode[1] = 0.5;
  testNodes.push_back(bufferNode);
  std::vector<int> bufferElement(3);
  bufferElement[0] = 0; bufferElement[1] = 5; bufferElement[2] = 4;
  testSimplexConnectivity.push_back(bufferElement);
  bufferElement[0] = 5; bufferElement[1] = 1;
  testSimplexConnectivity.push_back(bufferElement);
  bufferElement[0] = 1; bufferElement[1] = 6;
  testSimplexConnectivity.push_back(bufferElement);
  bufferElement[0] = 6; bufferElement[1] = 2;
  testSimplexConnectivity.push_back(bufferElement);
  bufferElement[0] = 2; bufferElement[1] = 7;
  testSimplexConnectivity.push_back(bufferElement);
  bufferElement[0] = 7; bufferElement[1] = 3;
  testSimplexConnectivity.push_back(bufferElement);
  bufferElement[0] = 3; bufferElement[1] = 8;
  testSimplexConnectivity.push_back(bufferElement);
  bufferElement[0] = 8; bufferElement[1] = 0;
  testSimplexConnectivity.push_back(bufferElement);
  bufferElement.resize(4);
  bufferElement[0] = 0; bufferElement[1] = 5; bufferElement[2] = 4; bufferElement[3] = 8;
  testOrthotopeConnectivity.push_back(bufferElement);
  bufferElement[0] = 1; bufferElement[1] = 6; bufferElement[2] = 4; bufferElement[3] = 5;
  testOrthotopeConnectivity.push_back(bufferElement);
  bufferElement[0] = 2; bufferElement[1] = 7; bufferElement[2] = 4; bufferElement[3] = 6;
  testOrthotopeConnectivity.push_back(bufferElement);
  bufferElement[0] = 3; bufferElement[1] = 8; bufferElement[2] = 4; bufferElement[3] = 7;
  testOrthotopeConnectivity.push_back(bufferElement);
  SECTION("Testing Construction."){
    CHECK_NOTHROW(Mesh());
    CHECK_NOTHROW(Mesh(2, 1, "simplex", testNodes, testSimplexConnectivity));
    CHECK_NOTHROW(Mesh(2, 1, "orthotope", testNodes, testOrthotopeConnectivity));
  }
  Mesh simplexMesh(2, 1, "simplex", testNodes, testSimplexConnectivity);
  Mesh orthoMesh(2, 1, "orthotope", testNodes, testOrthotopeConnectivity);
  SECTION("Testing easy gets"){
    CHECK(simplexMesh.getDimension() == 2);
    CHECK(simplexMesh.getNumberPoints() == testNodes.size());
    CHECK(simplexMesh.getNumberCells() == testSimplexConnectivity.size());
    CHECK(orthoMesh.getDimension() == 2);
    CHECK(orthoMesh.getNumberPoints() == testNodes.size());
    CHECK(orthoMesh.getNumberCells() == testOrthotopeConnectivity.size());
    std::vector<double> point;
    for(int i = 0; i < simplexMesh.getNumberPoints(); i++){
      simplexMesh.getPoint(i, &point);
      for(int j = 0; j < point.size(); j++){
        CHECK(point[j] == testNodes[i][j]);
      }
    }
    std::vector<int> cell;
    for(int i = 0; i < simplexMesh.getNumberCells(); i++){
      simplexMesh.getCell(i, &cell);
      for(int j = 0; j < cell.size(); j++){
        CHECK(cell[j] == testSimplexConnectivity[i][j]);
      }
    }
    for(int i = 0; i < orthoMesh.getNumberCells(); i++){
      orthoMesh.getCell(i, &cell);
      for(int j = 0; j < cell.size(); j++){
        CHECK(cell[j] == testOrthotopeConnectivity[i][j]);
      }
    }
  }
  std::vector< std::array<int, 4> > innerSimplexFaces;
  std::array<int, 4> bufferFace;
  bufferFace[0] = 0; bufferFace[1] = 1; bufferFace[2] = 1; bufferFace[3] = 2;
  innerSimplexFaces.push_back(bufferFace);
  bufferFace[0] = 1; bufferFace[1] = 1; bufferFace[2] = 2; bufferFace[3] = 2;
  innerSimplexFaces.push_back(bufferFace);
  bufferFace[0] = 2; bufferFace[1] = 1; bufferFace[2] = 3; bufferFace[3] = 2;
  innerSimplexFaces.push_back(bufferFace);
  bufferFace[0] = 3; bufferFace[1] = 1; bufferFace[2] = 4; bufferFace[3] = 2;
  innerSimplexFaces.push_back(bufferFace);
  bufferFace[0] = 4; bufferFace[1] = 1; bufferFace[2] = 5; bufferFace[3] = 2;
  innerSimplexFaces.push_back(bufferFace);
  bufferFace[0] = 5; bufferFace[1] = 1; bufferFace[2] = 6; bufferFace[3] = 2;
  innerSimplexFaces.push_back(bufferFace);
  bufferFace[0] = 6; bufferFace[1] = 1; bufferFace[2] = 7; bufferFace[3] = 2;
  innerSimplexFaces.push_back(bufferFace);
  bufferFace[0] = 7; bufferFace[1] = 1; bufferFace[2] = 0; bufferFace[3] = 2;
  innerSimplexFaces.push_back(bufferFace);
  std::vector< std::array<int, 2> > outerSimplexFaces;
  std::array<int, 2> bufferOFace;
  bufferOFace[0] = 0; bufferOFace[1] = 0;
  outerSimplexFaces.push_back(bufferOFace);
  bufferOFace[0] = 1; bufferOFace[1] = 0;
  outerSimplexFaces.push_back(bufferOFace);
  bufferOFace[0] = 2; bufferOFace[1] = 0;
  outerSimplexFaces.push_back(bufferOFace);
  bufferOFace[0] = 3; bufferOFace[1] = 0;
  outerSimplexFaces.push_back(bufferOFace);
  bufferOFace[0] = 4; bufferOFace[1] = 0;
  outerSimplexFaces.push_back(bufferOFace);
  bufferOFace[0] = 5; bufferOFace[1] = 0;
  outerSimplexFaces.push_back(bufferOFace);
  bufferOFace[0] = 6; bufferOFace[1] = 0;
  outerSimplexFaces.push_back(bufferOFace);
  bufferOFace[0] = 7; bufferOFace[1] = 0;
  outerSimplexFaces.push_back(bufferOFace);
  std::vector< std::array<int, 4> > innerOrthoFaces;
  bufferFace[0] = 0; bufferFace[1] = 1; bufferFace[2] = 1; bufferFace[3] = 2;
  innerOrthoFaces.push_back(bufferFace);
  bufferFace[0] = 1; bufferFace[1] = 1; bufferFace[2] = 2; bufferFace[3] = 2;
  innerOrthoFaces.push_back(bufferFace);
  bufferFace[0] = 2; bufferFace[1] = 1; bufferFace[2] = 3; bufferFace[3] = 2;
  innerOrthoFaces.push_back(bufferFace);
  bufferFace[0] = 3; bufferFace[1] = 1; bufferFace[2] = 0; bufferFace[3] = 2;
  innerOrthoFaces.push_back(bufferFace);
  std::vector< std::array<int, 2> > outerOrthoFaces;
  bufferOFace[0] = 0; bufferOFace[1] = 0;
  outerOrthoFaces.push_back(bufferOFace);
  bufferOFace[0] = 1; bufferOFace[1] = 3;
  outerOrthoFaces.push_back(bufferOFace);
  bufferOFace[0] = 1; bufferOFace[1] = 0;
  outerOrthoFaces.push_back(bufferOFace);
  bufferOFace[0] = 2; bufferOFace[1] = 3;
  outerOrthoFaces.push_back(bufferOFace);
  bufferOFace[0] = 2; bufferOFace[1] = 0;
  outerOrthoFaces.push_back(bufferOFace);
  bufferOFace[0] = 3; bufferOFace[1] = 3;
  outerOrthoFaces.push_back(bufferOFace);
  bufferOFace[0] = 3; bufferOFace[1] = 0;
  outerOrthoFaces.push_back(bufferOFace);
  bufferOFace[0] = 0; bufferOFace[1] = 3;
  outerOrthoFaces.push_back(bufferOFace);
  SECTION("Testing face computation algorithm"){
    std::vector< std::array<int, 4> > innerFaces;
    std::vector< std::array<int, 2> > outerFaces;
    bool areEqual;
    simplexMesh.getInnerFaces(&innerFaces);
    CHECK(innerFaces.size() == innerSimplexFaces.size());
    simplexMesh.getOuterFaces(&outerFaces);
    CHECK(outerFaces.size() == outerSimplexFaces.size());
    for(int i = 0; i < innerFaces.size(); i++){
      bufferFace = innerFaces[i];
      for(int j = 0; j < innerSimplexFaces.size(); j++){
        areEqual = arrayEqual4(bufferFace, innerSimplexFaces[j]);
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }
    for(int i = 0; i < outerFaces.size(); i++){
      bufferOFace = outerFaces[i];
      for(int j = 0; j < outerSimplexFaces.size(); j++){
        areEqual = arrayEqual2(bufferOFace, outerSimplexFaces[j]);
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }
    orthoMesh.getInnerFaces(&innerFaces);
    CHECK(innerFaces.size() == innerOrthoFaces.size());
    orthoMesh.getOuterFaces(&outerFaces);
    CHECK(outerFaces.size() == outerOrthoFaces.size());
    for(int i = 0; i < innerFaces.size(); i++){
      bufferFace = innerFaces[i];
      for(int j = 0; j < innerOrthoFaces.size(); j++){
        areEqual = arrayEqual4(bufferFace, innerOrthoFaces[j]);
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }
    for(int i = 0; i < outerFaces.size(); i++){
      bufferOFace = outerFaces[i];
      for(int j = 0; j < outerOrthoFaces.size(); j++){
        areEqual = arrayEqual2(bufferOFace, outerOrthoFaces[j]);
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }
  }

  testSimplexConnectivity.resize(0);
  testOrthotopeConnectivity.resize(0);
  bufferElement.resize(6);
  bufferElement[0] = 0; bufferElement[1] = 5; bufferElement[2] = 1;
  bufferElement[3] = 3; bufferElement[4] = 4; bufferElement[5] = 8;
  testSimplexConnectivity.push_back(bufferElement);
  bufferElement[0] = 3; bufferElement[1] = 4; bufferElement[2] = 1;
  bufferElement[3] = 2; bufferElement[4] = 6; bufferElement[5] = 7;
  testSimplexConnectivity.push_back(bufferElement);
  bufferElement.resize(9);
  bufferElement[0] = 0; bufferElement[1] = 5; bufferElement[2] = 1;
  bufferElement[3] = 2; bufferElement[4] = 6; bufferElement[5] = 7;
  bufferElement[6] = 3; bufferElement[7] = 8; bufferElement[8] = 4;
  testOrthotopeConnectivity.push_back(bufferElement);

  SECTION("Higher order mesh construction."){
    CHECK_NOTHROW(Mesh(2, 2, "simplex", testNodes, testSimplexConnectivity));
    CHECK_NOTHROW(Mesh(2, 2, "orthotope", testNodes, testOrthotopeConnectivity));
  }

  Mesh simplexMesh2(2, 2, "simplex", testNodes, testSimplexConnectivity);
  Mesh orthoMesh2(2, 2, "orthotope", testNodes, testOrthotopeConnectivity);
  SECTION("Testing easy gets"){
    CHECK(simplexMesh2.getDimension() == 2);
    CHECK(simplexMesh2.getNumberPoints() == testNodes.size());
    CHECK(simplexMesh2.getNumberCells() == testSimplexConnectivity.size());
    CHECK(orthoMesh2.getDimension() == 2);
    CHECK(orthoMesh2.getNumberPoints() == testNodes.size());
    CHECK(orthoMesh2.getNumberCells() == testOrthotopeConnectivity.size());
    std::vector<double> point;
    for(int i = 0; i < simplexMesh2.getNumberPoints(); i++){
      simplexMesh2.getPoint(i, &point);
      for(int j = 0; j < point.size(); j++){
        CHECK(point[j] == testNodes[i][j]);
      }
    }
    std::vector<int> cell;
    for(int i = 0; i < simplexMesh2.getNumberCells(); i++){
      simplexMesh2.getCell(i, &cell);
      for(int j = 0; j < cell.size(); j++){
        CHECK(cell[j] == testSimplexConnectivity[i][j]);
      }
    }
    for(int i = 0; i < orthoMesh2.getNumberCells(); i++){
      orthoMesh2.getCell(i, &cell);
      for(int j = 0; j < cell.size(); j++){
        CHECK(cell[j] == testOrthotopeConnectivity[i][j]);
      }
    }
  }

  innerSimplexFaces.resize(0);
  outerSimplexFaces.resize(0);
  bufferFace[0] = 1; bufferFace[1] = 0; bufferFace[2] = 0; bufferFace[3] = 1;
  innerSimplexFaces.push_back(bufferFace);
  bufferOFace[0] = 0; bufferOFace[1] = 0;
  outerSimplexFaces.push_back(bufferOFace);
  bufferOFace[0] = 1; bufferOFace[1] = 1;
  outerSimplexFaces.push_back(bufferOFace);
  bufferOFace[0] = 1; bufferOFace[1] = 2;
  outerSimplexFaces.push_back(bufferOFace);
  bufferOFace[0] = 0; bufferOFace[1] = 2;
  outerSimplexFaces.push_back(bufferOFace);
  innerOrthoFaces.resize(0);
  outerOrthoFaces.resize(0);
  bufferOFace[0] = 0; bufferOFace[1] = 0;
  outerOrthoFaces.push_back(bufferOFace);
  bufferOFace[0] = 0; bufferOFace[1] = 1;
  outerOrthoFaces.push_back(bufferOFace);
  bufferOFace[0] = 0; bufferOFace[1] = 2;
  outerOrthoFaces.push_back(bufferOFace);
  bufferOFace[0] = 0; bufferOFace[1] = 3;
  outerOrthoFaces.push_back(bufferOFace);
  SECTION("Test higher order face computation"){
    std::vector< std::array<int, 4> > innerFaces;
    std::vector< std::array<int, 2> > outerFaces;
    bool areEqual;
    simplexMesh2.getInnerFaces(&innerFaces);
    CHECK(innerFaces.size() == innerSimplexFaces.size());
    simplexMesh2.getOuterFaces(&outerFaces);
    CHECK(outerFaces.size() == outerSimplexFaces.size());
    for(int i = 0; i < innerFaces.size(); i++){
      bufferFace = innerFaces[i];
      for(int j = 0; j < innerSimplexFaces.size(); j++){
        areEqual = arrayEqual4(bufferFace, innerSimplexFaces[j]);
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }
    for(int i = 0; i < outerFaces.size(); i++){
      bufferOFace = outerFaces[i];
      for(int j = 0; j < outerSimplexFaces.size(); j++){
        areEqual = arrayEqual2(bufferOFace, outerSimplexFaces[j]);
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }
    orthoMesh2.getInnerFaces(&innerFaces);
    CHECK(innerFaces.size() == innerOrthoFaces.size());
    orthoMesh2.getOuterFaces(&outerFaces);
    CHECK(outerFaces.size() == outerOrthoFaces.size());
    for(int i = 0; i < innerFaces.size(); i++){
      bufferFace = innerFaces[i];
      for(int j = 0; j < innerOrthoFaces.size(); j++){
        areEqual = arrayEqual4(bufferFace, innerOrthoFaces[j]);
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }
    for(int i = 0; i < outerFaces.size(); i++){
      bufferOFace = outerFaces[i];
      for(int j = 0; j < outerOrthoFaces.size(); j++){
        areEqual = arrayEqual2(bufferOFace, outerOrthoFaces[j]);
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }
  }
};
