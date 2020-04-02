#include <catch2/catch.hpp>
#include <vector>
#include <algorithm>
#include "Mesh.h"

using namespace hfox;

bool vectorEqual(const std::vector<int> & a, const std::vector<int> & b){
  bool areEqual = 1;
  for(int i = 0; i < a.size(); i++){
    areEqual = (a[i] == b[i]);
    if(!areEqual){
      break;
    }
  }
  return areEqual;
}

//bool arrayEqual4(const std::array<int, 4> & a, const std::array<int, 4> & b){
  //bool areEqual = 1;
  //for(int i = 0; i < 4; i++){
    //areEqual = (a[i] == b[i]);
    //if(!areEqual){
      //break;
    //}
  //}
  //return areEqual;
//}

//bool arrayEqual2(const std::array<int, 2> & a, const std::array<int, 2> & b){
  //bool areEqual = 1;
  //for(int i = 0; i < 2; i++){
    //areEqual = (a[i] == b[i]);
    //if(!areEqual){
      //break;
    //}
  //}
  //return areEqual;
//}

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

  std::vector< std::vector<int> > simplexFaces;
  std::vector< std::vector<int> > orthoFaces;
  bufferElement.resize(2);
  bufferElement[0] = 0; bufferElement[1] = 5;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  bufferElement[0] = 5; bufferElement[1] = 1;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  bufferElement[0] = 1; bufferElement[1] = 6;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  bufferElement[0] = 6; bufferElement[1] = 2;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  bufferElement[0] = 2; bufferElement[1] = 7;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  bufferElement[0] = 7; bufferElement[1] = 3;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  bufferElement[0] = 3; bufferElement[1] = 8;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  bufferElement[0] = 8; bufferElement[1] = 0;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  bufferElement[0] = 4; bufferElement[1] = 0;
  simplexFaces.push_back(bufferElement);
  bufferElement[0] = 4; bufferElement[1] = 5;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  bufferElement[0] = 4; bufferElement[1] = 1;
  simplexFaces.push_back(bufferElement);
  bufferElement[0] = 4; bufferElement[1] = 6;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  bufferElement[0] = 4; bufferElement[1] = 2;
  simplexFaces.push_back(bufferElement);
  bufferElement[0] = 4; bufferElement[1] = 7;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  bufferElement[0] = 4; bufferElement[1] = 3;
  simplexFaces.push_back(bufferElement);
  bufferElement[0] = 4; bufferElement[1] = 8;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  std::vector< std::vector<int> > simplexFaceToCellMap;
  std::vector< std::vector<int> > orthoFaceToCellMap;
  bufferElement.resize(2);
  bufferElement[0] = 0; bufferElement[1] = 1;
  simplexFaceToCellMap.push_back(bufferElement);
  orthoFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 1; bufferElement[1] = 2;
  simplexFaceToCellMap.push_back(bufferElement);
  orthoFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 2; bufferElement[1] = 3;
  simplexFaceToCellMap.push_back(bufferElement);
  orthoFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 3; bufferElement[1] = 4;
  simplexFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 4; bufferElement[1] = 5;
  simplexFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 5; bufferElement[1] = 6;
  simplexFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 6; bufferElement[1] = 7;
  simplexFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 7; bufferElement[1] = 0;
  simplexFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 3; bufferElement[1] = 0;
  orthoFaceToCellMap.push_back(bufferElement);
  bufferElement.resize(1);
  bufferElement[0] = 0;
  simplexFaceToCellMap.push_back(bufferElement);
  orthoFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 1;
  simplexFaceToCellMap.push_back(bufferElement);
  orthoFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 2;
  simplexFaceToCellMap.push_back(bufferElement);
  orthoFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 3;
  simplexFaceToCellMap.push_back(bufferElement);
  orthoFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 4;
  simplexFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 5;
  simplexFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 6;
  simplexFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 7;
  simplexFaceToCellMap.push_back(bufferElement);

  std::vector<int> boundaryNodes = {0, 5, 1, 6, 2, 7, 3, 8};
  SECTION("Testing face computation"){
    CHECK(simplexMesh.getNumberFaces() == simplexFaces.size());
    CHECK(orthoMesh.getNumberFaces() == orthoFaces.size());
    bool areEqual;
    for(int i = 0; i < simplexFaces.size(); i++){
      simplexMesh.getFace(i, &bufferElement);
      for(int j = 0; j < simplexFaces.size(); j++){
        for(int k = 0; k < bufferElement.size(); k++){
          areEqual = (std::find(simplexFaces[j].begin(), simplexFaces[j].end(), bufferElement[k]) != simplexFaces[j].end());
          if(!areEqual){
            break;
          }
        }
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }
    for(int i = 0; i < orthoFaces.size(); i++){
      orthoMesh.getFace(i, &bufferElement);
      for(int j = 0; j < orthoFaces.size(); j++){
        for(int k = 0; k < bufferElement.size(); k++){
          areEqual = (std::find(orthoFaces[j].begin(), orthoFaces[j].end(), bufferElement[k]) != orthoFaces[j].end());
          if(!areEqual){
            break;
          }
        }
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }

    CHECK(simplexMesh.getFace2CellMap()->size() == simplexFaces.size());
    CHECK(orthoMesh.getFace2CellMap()->size() == orthoFaces.size());
    for(int i = 0; i < simplexFaces.size(); i++){
      const std::vector<int> * bufferFace = simplexMesh.getFace2Cell(i);
      for(int j = 0; j < simplexFaceToCellMap.size(); j++){
        if(bufferFace->size() == simplexFaceToCellMap[j].size()){
          for(int k = 0; k < bufferFace->size(); k++){
            areEqual = (std::find(simplexFaceToCellMap[j].begin(), simplexFaceToCellMap[j].end(), (*bufferFace)[k]) != simplexFaceToCellMap[j].end());
            if(!areEqual){
              break;
            }
          }
        } else {
          areEqual = 0;
        }
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }
    for(int i = 0; i < orthoFaces.size(); i++){
      const std::vector<int> * bufferFace = orthoMesh.getFace2Cell(i);
      for(int j = 0; j < orthoFaceToCellMap.size(); j++){
        if(bufferFace->size() == orthoFaceToCellMap[j].size()){
          for(int k = 0; k < bufferFace->size(); k++){
            areEqual = (std::find(orthoFaceToCellMap[j].begin(), orthoFaceToCellMap[j].end(), (*bufferFace)[k]) != orthoFaceToCellMap[j].end());
            if(!areEqual){
              break;
            }
          }
        } else {
          areEqual = 0;
        }
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }
    CHECK(simplexMesh.getCell2FaceMap()->size() == testSimplexConnectivity.size());
    CHECK(orthoMesh.getCell2FaceMap()->size() == testOrthotopeConnectivity.size());
    bool isIn;
    for(int i = 0; i < testSimplexConnectivity.size(); i++){
      const std::vector<int> * buff = simplexMesh.getCell2Face(i); 
      for(int j = 0; j < buff->size(); j++){
        simplexMesh.getFace((*buff)[j], &bufferElement);
        for(int k = 0; k < bufferElement.size(); k++){
          CHECK(std::find(testSimplexConnectivity[i].begin(), testSimplexConnectivity[i].end(), bufferElement[k]) != testSimplexConnectivity[i].end());
        }
      }
    }
    for(int i = 0; i < testOrthotopeConnectivity.size(); i++){
      const std::vector<int> * buff = orthoMesh.getCell2Face(i); 
      for(int j = 0; j < buff->size(); j++){
        orthoMesh.getFace((*buff)[j], &bufferElement);
        for(int k = 0; k < bufferElement.size(); k++){
          CHECK(std::find(testOrthotopeConnectivity[i].begin(), testOrthotopeConnectivity[i].end(), bufferElement[k]) != testOrthotopeConnectivity[i].end());
        }
      }
    }
    CHECK(simplexMesh.getBoundaryFaces()->size() == 8);
    CHECK(orthoMesh.getBoundaryFaces()->size() == 8);
    const std::set<int> * buffBoundary = simplexMesh.getBoundaryFaces();
    for(std::set<int>::iterator it = buffBoundary->begin(); it != buffBoundary->end(); it++){
      simplexMesh.getFace(*it, &bufferElement);
      for(int k = 0; k < bufferElement.size(); k++){
        CHECK((std::find(boundaryNodes.begin(), boundaryNodes.end(), bufferElement[k]) != boundaryNodes.end()));
      }
    }
    buffBoundary = orthoMesh.getBoundaryFaces();
    for(std::set<int>::iterator it = buffBoundary->begin(); it != buffBoundary->end(); it++){
      orthoMesh.getFace(*it, &bufferElement);
      for(int k = 0; k < bufferElement.size(); k++){
        CHECK((std::find(boundaryNodes.begin(), boundaryNodes.end(), bufferElement[k]) != boundaryNodes.end()));
      }
    }
  };

  testSimplexConnectivity.resize(0);
  testOrthotopeConnectivity.resize(0);
  bufferElement.resize(6);
  bufferElement[0] = 0; bufferElement[1] = 1; bufferElement[2] = 3;
  bufferElement[3] = 5; bufferElement[4] = 4; bufferElement[5] = 8;
  testSimplexConnectivity.push_back(bufferElement);
  bufferElement[0] = 3; bufferElement[1] = 1; bufferElement[2] = 2;
  bufferElement[3] = 4; bufferElement[4] = 6; bufferElement[5] = 7;
  testSimplexConnectivity.push_back(bufferElement);
  bufferElement.resize(9);
  bufferElement[0] = 0; bufferElement[1] = 1; bufferElement[2] = 2;
  bufferElement[3] = 3; bufferElement[4] = 5; bufferElement[5] = 6;
  bufferElement[6] = 7; bufferElement[7] = 8; bufferElement[8] = 4;
  testOrthotopeConnectivity.push_back(bufferElement);

  SECTION("Higher order mesh construction."){
    CHECK_NOTHROW(Mesh(2, 2, "simplex", testNodes, testSimplexConnectivity));
    CHECK_NOTHROW(Mesh(2, 2, "orthotope", testNodes, testOrthotopeConnectivity));
  };

  Mesh simplexMesh2(2, 2, "simplex", testNodes, testSimplexConnectivity);
  Mesh orthoMesh2(2, 2, "orthotope", testNodes, testOrthotopeConnectivity);
  SECTION("Testing higher order easy gets"){
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
  };

  simplexFaces.resize(0);
  orthoFaces.resize(0);
  bufferElement.resize(3);
  bufferElement[0] = 0; bufferElement[1] = 1; bufferElement[2] = 5;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  bufferElement[0] = 1; bufferElement[1] = 2; bufferElement[2] = 6;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  bufferElement[0] = 2; bufferElement[1] = 3; bufferElement[2] = 7;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  bufferElement[0] = 3; bufferElement[1] = 0; bufferElement[2] = 8;
  simplexFaces.push_back(bufferElement);
  orthoFaces.push_back(bufferElement);
  bufferElement[0] = 1; bufferElement[1] = 3; bufferElement[2] = 4;
  simplexFaces.push_back(bufferElement);
  simplexFaceToCellMap.resize(0);
  orthoFaceToCellMap.resize(0);
  bufferElement.resize(2);
  bufferElement[0] = 0; bufferElement[1] = 1;
  simplexFaceToCellMap.push_back(bufferElement);
  bufferElement.resize(1);
  bufferElement[0] = 0;
  simplexFaceToCellMap.push_back(bufferElement);
  orthoFaceToCellMap.push_back(bufferElement);
  bufferElement[0] = 1;
  simplexFaceToCellMap.push_back(bufferElement);
  SECTION("Testing higher order face computation"){
    CHECK(simplexMesh2.getNumberFaces() == simplexFaces.size());
    CHECK(orthoMesh2.getNumberFaces() == orthoFaces.size());
    bool areEqual;
    for(int i = 0; i < simplexFaces.size(); i++){
      simplexMesh2.getFace(i, &bufferElement);
      for(int j = 0; j < simplexFaces.size(); j++){
        for(int k = 0; k < bufferElement.size(); k++){
          areEqual = (std::find(simplexFaces[j].begin(), simplexFaces[j].end(), bufferElement[k]) != simplexFaces[j].end());
          if(!areEqual){
            break;
          }
        }
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }
    for(int i = 0; i < orthoFaces.size(); i++){
      orthoMesh2.getFace(i, &bufferElement);
      for(int j = 0; j < orthoFaces.size(); j++){
        for(int k = 0; k < bufferElement.size(); k++){
          areEqual = (std::find(orthoFaces[j].begin(), orthoFaces[j].end(), bufferElement[k]) != orthoFaces[j].end());
          if(!areEqual){
            break;
          }
        }
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }
    CHECK(simplexMesh2.getFace2CellMap()->size() == simplexFaces.size());
    CHECK(orthoMesh2.getFace2CellMap()->size() == orthoFaces.size());
    for(int i = 0; i < simplexFaces.size(); i++){
      const std::vector<int> * bufferFace = simplexMesh2.getFace2Cell(i);
      for(int j = 0; j < simplexFaceToCellMap.size(); j++){
        if(bufferFace->size() == simplexFaceToCellMap[j].size()){
          for(int k = 0; k < bufferFace->size(); k++){
            areEqual = (std::find(simplexFaceToCellMap[j].begin(), simplexFaceToCellMap[j].end(), (*bufferFace)[k]) != simplexFaceToCellMap[j].end());
            if(!areEqual){
              break;
            }
          }
        } else {
          areEqual = 0;
        }
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }
    for(int i = 0; i < orthoFaces.size(); i++){
      const std::vector<int> * bufferFace = orthoMesh2.getFace2Cell(i);
      for(int j = 0; j < orthoFaceToCellMap.size(); j++){
        if(bufferFace->size() == orthoFaceToCellMap[j].size()){
          for(int k = 0; k < bufferFace->size(); k++){
            areEqual = (std::find(orthoFaceToCellMap[j].begin(), orthoFaceToCellMap[j].end(), (*bufferFace)[k]) != orthoFaceToCellMap[j].end());
            if(!areEqual){
              break;
            }
          }
        } else {
          areEqual = 0;
        }
        if(areEqual){
          break;
        }
      }
      CHECK(areEqual);
    }
    CHECK(simplexMesh2.getCell2FaceMap()->size() == testSimplexConnectivity.size());
    CHECK(orthoMesh2.getCell2FaceMap()->size() == testOrthotopeConnectivity.size());
    bool isIn;
    for(int i = 0; i < testSimplexConnectivity.size(); i++){
      const std::vector<int> * buff = simplexMesh2.getCell2Face(i); 
      for(int j = 0; j < buff->size(); j++){
        simplexMesh2.getFace((*buff)[j], &bufferElement);
        for(int k = 0; k < bufferElement.size(); k++){
          CHECK(std::find(testSimplexConnectivity[i].begin(), testSimplexConnectivity[i].end(), bufferElement[k]) != testSimplexConnectivity[i].end());
        }
      }
    }
    for(int i = 0; i < testOrthotopeConnectivity.size(); i++){
      const std::vector<int> * buff = orthoMesh2.getCell2Face(i); 
      for(int j = 0; j < buff->size(); j++){
        orthoMesh2.getFace((*buff)[j], &bufferElement);
        for(int k = 0; k < bufferElement.size(); k++){
          CHECK(std::find(testOrthotopeConnectivity[i].begin(), testOrthotopeConnectivity[i].end(), bufferElement[k]) != testOrthotopeConnectivity[i].end());
        }
      }
    }
    CHECK(simplexMesh2.getBoundaryFaces()->size() == 4);
    CHECK(orthoMesh2.getBoundaryFaces()->size() == 4);
    const std::set<int> * buffBoundary = simplexMesh2.getBoundaryFaces();
    for(std::set<int>::iterator it = buffBoundary->begin(); it != buffBoundary->end(); it++){
      simplexMesh2.getFace(*it, &bufferElement);
      for(int k = 0; k < bufferElement.size(); k++){
        CHECK((std::find(boundaryNodes.begin(), boundaryNodes.end(), bufferElement[k]) != boundaryNodes.end()));
      }
    }
    buffBoundary = orthoMesh2.getBoundaryFaces();
    for(std::set<int>::iterator it = buffBoundary->begin(); it != buffBoundary->end(); it++){
      orthoMesh2.getFace(*it, &bufferElement);
      for(int k = 0; k < bufferElement.size(); k++){
        CHECK((std::find(boundaryNodes.begin(), boundaryNodes.end(), bufferElement[k]) != boundaryNodes.end()));
      }
    }
  };
};
