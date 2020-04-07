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

template <class T>
std::vector<T> unpack(std::vector< std::vector<T> > v){
  std::vector<T> res(v.size()*v[0].size());
  int index = 0;
  for(int i = 0; i < v.size(); i++){
    for(int j = 0; j < v[i].size(); j++){
      res[index] = v[i][j];
      index += 1;
    }
  }
  return res;
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
  std::vector<double> testNodesStream = unpack(testNodes);
  std::vector<int> testSimplexConnectivityStream = unpack(testSimplexConnectivity);
  std::vector<int> testOrthotopeConnectivityStream = unpack(testOrthotopeConnectivity);
  SECTION("Testing Construction."){
    CHECK_NOTHROW(Mesh());
    CHECK_NOTHROW(Mesh(2, 1, "simplex", 2, testNodesStream, testSimplexConnectivityStream));
    CHECK_NOTHROW(Mesh(2, 1, "orthotope", 2, testNodesStream, testOrthotopeConnectivityStream));
  }
  Mesh simplexMesh(2, 1, "simplex", 2, testNodesStream, testSimplexConnectivityStream);
  Mesh orthoMesh(2, 1, "orthotope", 2, testNodesStream, testOrthotopeConnectivityStream);
  SECTION("Testing easy gets"){
    CHECK(simplexMesh.getDimension() == 2);
    CHECK(simplexMesh.getNodeSpaceDimension() == 2);
    CHECK(simplexMesh.getNumberPoints() == testNodes.size());
    CHECK(simplexMesh.getNumberCells() == testSimplexConnectivity.size());
    CHECK(orthoMesh.getDimension() == 2);
    CHECK(orthoMesh.getNodeSpaceDimension() == 2);
    CHECK(orthoMesh.getNumberPoints() == testNodes.size());
    CHECK(orthoMesh.getNumberCells() == testOrthotopeConnectivity.size());
    const double * point;
    for(int i = 0; i < simplexMesh.getNumberPoints(); i++){
      point = simplexMesh.getPoint(i);
      for(int j = 0; j < simplexMesh.getNodeSpaceDimension(); j++){
        CHECK((point[j]) == testNodes[i][j]);
      }
    }
    const int * cell;
    for(int i = 0; i < simplexMesh.getNumberCells(); i++){
      cell = simplexMesh.getCell(i);
      for(int j = 0; j < simplexMesh.getReferenceElement()->getNumNodes(); j++){
        CHECK((cell[j]) == testSimplexConnectivity[i][j]);
      }
    }
    for(int i = 0; i < orthoMesh.getNumberCells(); i++){
      cell = orthoMesh.getCell(i);
      for(int j = 0; j < orthoMesh.getReferenceElement()->getNumNodes(); j++){
        CHECK((cell[j]) == testOrthotopeConnectivity[i][j]);
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
    const int * buffer;
    for(int i = 0; i < simplexFaces.size(); i++){
      buffer = simplexMesh.getFace(i);
      for(int j = 0; j < simplexFaces.size(); j++){
        for(int k = 0; k < simplexMesh.getReferenceElement()->getFaceElement()->getNumNodes(); k++){
          areEqual = (std::find(simplexFaces[j].begin(), simplexFaces[j].end(), (buffer[k])) != simplexFaces[j].end());
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
      buffer = orthoMesh.getFace(i);
      for(int j = 0; j < orthoFaces.size(); j++){
        for(int k = 0; k < orthoMesh.getReferenceElement()->getFaceElement()->getNumNodes(); k++){
          areEqual = (std::find(orthoFaces[j].begin(), orthoFaces[j].end(), (buffer[k])) != orthoFaces[j].end());
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
    //for(int i = 0; i < simplexFaces.size(); i++){
      //std::vector<const int> bufferFace(2);
      //bufferFace.data() = simplexMesh.getFace2Cell(i);
      //for(int j = 0; j < simplexFaceToCellMap.size(); j++){
        //if(bufferFace->size() == simplexFaceToCellMap[j].size()){
          //for(int k = 0; k < bufferFace->size(); k++){
            //areEqual = (std::find(simplexFaceToCellMap[j].begin(), simplexFaceToCellMap[j].end(), (*bufferFace)[k]) != simplexFaceToCellMap[j].end());
            //if(!areEqual){
              //break;
            //}
          //}
        //} else {
          //areEqual = 0;
        //}
        //if(areEqual){
          //break;
        //}
      //}
      //CHECK(areEqual);
    //}
    //for(int i = 0; i < orthoFaces.size(); i++){
      //const std::vector<int> * bufferFace = orthoMesh.getFace2Cell(i);
      //for(int j = 0; j < orthoFaceToCellMap.size(); j++){
        //if(bufferFace->size() == orthoFaceToCellMap[j].size()){
          //for(int k = 0; k < bufferFace->size(); k++){
            //areEqual = (std::find(orthoFaceToCellMap[j].begin(), orthoFaceToCellMap[j].end(), (*bufferFace)[k]) != orthoFaceToCellMap[j].end());
            //if(!areEqual){
              //break;
            //}
          //}
        //} else {
          //areEqual = 0;
        //}
        //if(areEqual){
          //break;
        //}
      //}
      //CHECK(areEqual);
    //}
    CHECK(simplexMesh.getCell2FaceMap()->size() == testSimplexConnectivity.size());
    CHECK(orthoMesh.getCell2FaceMap()->size() == testOrthotopeConnectivity.size());
    //bool isIn;
    //for(int i = 0; i < testSimplexConnectivity.size(); i++){
      //const std::vector<int> * buff = simplexMesh.getCell2Face(i); 
      //for(int j = 0; j < buff->size(); j++){
        //buffer = simplexMesh.getFace((*buff)[j]);
        //for(int k = 0; k < buffer->size(); k++){
          //CHECK(std::find(testSimplexConnectivity[i].begin(), testSimplexConnectivity[i].end(), (*buffer)[k]) != testSimplexConnectivity[i].end());
        //}
      //}
    //}
    //for(int i = 0; i < testOrthotopeConnectivity.size(); i++){
      //const std::vector<int> * buff = orthoMesh.getCell2Face(i); 
      //for(int j = 0; j < buff->size(); j++){
        //buffer = orthoMesh.getFace((*buff)[j]);
        //for(int k = 0; k < buffer->size(); k++){
          //CHECK(std::find(testOrthotopeConnectivity[i].begin(), testOrthotopeConnectivity[i].end(), (*buffer)[k]) != testOrthotopeConnectivity[i].end());
        //}
      //}
    //}
    CHECK(simplexMesh.getBoundaryFaces()->size() == 8);
    CHECK(orthoMesh.getBoundaryFaces()->size() == 8);
    //const std::set<int> * buffBoundary = simplexMesh.getBoundaryFaces();
    //for(std::set<int>::iterator it = buffBoundary->begin(); it != buffBoundary->end(); it++){
      //buffer = simplexMesh.getFace(*it);
      //for(int k = 0; k < buffer->size(); k++){
        //CHECK((std::find(boundaryNodes.begin(), boundaryNodes.end(), (*buffer)[k]) != boundaryNodes.end()));
      //}
    //}
    //buffBoundary = orthoMesh.getBoundaryFaces();
    //for(std::set<int>::iterator it = buffBoundary->begin(); it != buffBoundary->end(); it++){
      //buffer = orthoMesh.getFace(*it);
      //for(int k = 0; k < buffer->size(); k++){
        //CHECK((std::find(boundaryNodes.begin(), boundaryNodes.end(), (*buffer)[k]) != boundaryNodes.end()));
      //}
    //}
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
  testSimplexConnectivityStream = unpack(testSimplexConnectivity);
  testOrthotopeConnectivityStream = unpack(testOrthotopeConnectivity);
  SECTION("Higher order mesh construction."){
    CHECK_NOTHROW(Mesh(2, 2, "simplex", 2, testNodesStream, testSimplexConnectivityStream));
    CHECK_NOTHROW(Mesh(2, 2, "orthotope", 2, testNodesStream, testOrthotopeConnectivityStream));
  };

  Mesh simplexMesh2(2, 2, "simplex", 2, testNodesStream, testSimplexConnectivityStream);
  Mesh orthoMesh2(2, 2, "orthotope", 2, testNodesStream, testOrthotopeConnectivityStream);
  SECTION("Testing higher order easy gets"){
    CHECK(simplexMesh2.getDimension() == 2);
    CHECK(simplexMesh2.getNumberPoints() == testNodes.size());
    CHECK(simplexMesh2.getNumberCells() == testSimplexConnectivity.size());
    CHECK(orthoMesh2.getDimension() == 2);
    CHECK(orthoMesh2.getNumberPoints() == testNodes.size());
    CHECK(orthoMesh2.getNumberCells() == testOrthotopeConnectivity.size());
    const double * point;
    for(int i = 0; i < simplexMesh2.getNumberPoints(); i++){
      point = simplexMesh2.getPoint(i);
      for(int j = 0; j < simplexMesh2.getNodeSpaceDimension(); j++){
        CHECK((point)[j] == testNodes[i][j]);
      }
    }
    const int * cell;
    for(int i = 0; i < simplexMesh2.getNumberCells(); i++){
      cell = simplexMesh2.getCell(i);
      for(int j = 0; j < simplexMesh2.getReferenceElement()->getNumNodes(); j++){
        CHECK((cell)[j] == testSimplexConnectivity[i][j]);
      }
    }
    for(int i = 0; i < orthoMesh2.getNumberCells(); i++){
      cell = orthoMesh2.getCell(i);
      for(int j = 0; j < orthoMesh2.getReferenceElement()->getNumNodes(); j++){
        CHECK((cell)[j] == testOrthotopeConnectivity[i][j]);
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
    //bool areEqual;
    //const std::vector<int> * buffer;
    //for(int i = 0; i < simplexFaces.size(); i++){
      //buffer = simplexMesh2.getFace(i);
      //for(int j = 0; j < simplexFaces.size(); j++){
        //for(int k = 0; k < buffer->size(); k++){
          //areEqual = (std::find(simplexFaces[j].begin(), simplexFaces[j].end(), (*buffer)[k]) != simplexFaces[j].end());
          //if(!areEqual){
            //break;
          //}
        //}
        //if(areEqual){
          //break;
        //}
      //}
      //CHECK(areEqual);
    //}
    //for(int i = 0; i < orthoFaces.size(); i++){
      //buffer = orthoMesh2.getFace(i);
      //for(int j = 0; j < orthoFaces.size(); j++){
        //for(int k = 0; k < buffer->size(); k++){
          //areEqual = (std::find(orthoFaces[j].begin(), orthoFaces[j].end(), (*buffer)[k]) != orthoFaces[j].end());
          //if(!areEqual){
            //break;
          //}
        //}
        //if(areEqual){
          //break;
        //}
      //}
      //CHECK(areEqual);
    //}
    CHECK(simplexMesh2.getFace2CellMap()->size() == simplexFaces.size());
    CHECK(orthoMesh2.getFace2CellMap()->size() == orthoFaces.size());
    //for(int i = 0; i < simplexFaces.size(); i++){
      //const std::vector<int> * bufferFace = simplexMesh2.getFace2Cell(i);
      //for(int j = 0; j < simplexFaceToCellMap.size(); j++){
        //if(bufferFace->size() == simplexFaceToCellMap[j].size()){
          //for(int k = 0; k < bufferFace->size(); k++){
            //areEqual = (std::find(simplexFaceToCellMap[j].begin(), simplexFaceToCellMap[j].end(), (*bufferFace)[k]) != simplexFaceToCellMap[j].end());
            //if(!areEqual){
              //break;
            //}
          //}
        //} else {
          //areEqual = 0;
        //}
        //if(areEqual){
          //break;
        //}
      //}
      //CHECK(areEqual);
    //}
    //for(int i = 0; i < orthoFaces.size(); i++){
      //const std::vector<int> * bufferFace = orthoMesh2.getFace2Cell(i);
      //for(int j = 0; j < orthoFaceToCellMap.size(); j++){
        //if(bufferFace->size() == orthoFaceToCellMap[j].size()){
          //for(int k = 0; k < bufferFace->size(); k++){
            //areEqual = (std::find(orthoFaceToCellMap[j].begin(), orthoFaceToCellMap[j].end(), (*bufferFace)[k]) != orthoFaceToCellMap[j].end());
            //if(!areEqual){
              //break;
            //}
          //}
        //} else {
          //areEqual = 0;
        //}
        //if(areEqual){
          //break;
        //}
      //}
      //CHECK(areEqual);
    //}
    CHECK(simplexMesh2.getCell2FaceMap()->size() == testSimplexConnectivity.size());
    CHECK(orthoMesh2.getCell2FaceMap()->size() == testOrthotopeConnectivity.size());
    //bool isIn;
    //for(int i = 0; i < testSimplexConnectivity.size(); i++){
      //const std::vector<int> * buff = simplexMesh2.getCell2Face(i); 
      //for(int j = 0; j < buff->size(); j++){
        //buffer = simplexMesh2.getFace((*buff)[j]);
        //for(int k = 0; k < buffer->size(); k++){
          //CHECK(std::find(testSimplexConnectivity[i].begin(), testSimplexConnectivity[i].end(), (*buffer)[k]) != testSimplexConnectivity[i].end());
        //}
      //}
    //}
    //for(int i = 0; i < testOrthotopeConnectivity.size(); i++){
      //const std::vector<int> * buff = orthoMesh2.getCell2Face(i); 
      //for(int j = 0; j < buff->size(); j++){
        //buffer = orthoMesh2.getFace((*buff)[j]);
        //for(int k = 0; k < bufferElement.size(); k++){
          //CHECK(std::find(testOrthotopeConnectivity[i].begin(), testOrthotopeConnectivity[i].end(), (*buffer)[k]) != testOrthotopeConnectivity[i].end());
        //}
      //}
    //}
    CHECK(simplexMesh2.getBoundaryFaces()->size() == 4);
    CHECK(orthoMesh2.getBoundaryFaces()->size() == 4);
    //const std::set<int> * buffBoundary = simplexMesh2.getBoundaryFaces();
    //for(std::set<int>::iterator it = buffBoundary->begin(); it != buffBoundary->end(); it++){
      //buffer = simplexMesh2.getFace(*it);
      //for(int k = 0; k < bufferElement.size(); k++){
        //CHECK((std::find(boundaryNodes.begin(), boundaryNodes.end(), (*buffer)[k]) != boundaryNodes.end()));
      //}
    //}
    //buffBoundary = orthoMesh2.getBoundaryFaces();
    //for(std::set<int>::iterator it = buffBoundary->begin(); it != buffBoundary->end(); it++){
      //buffer = orthoMesh2.getFace(*it);
      //for(int k = 0; k < bufferElement.size(); k++){
        //CHECK((std::find(boundaryNodes.begin(), boundaryNodes.end(), (*buffer)[k]) != boundaryNodes.end()));
      //}
    //}
  };

  std::vector<double> nodesOrd3({0, 0,
1, 0,
1, 1,
0, 1,
0.499999999998694, 0,
0.1666666666663331, 0,
0.3333333333325136, 0,
0.666666666665796, 0,
0.8333333333328979, 0,
1, 0.499999999998694,
1, 0.1666666666663331,
1, 0.3333333333325136,
1, 0.666666666665796,
1, 0.8333333333328979,
0.5000000000020591, 1,
0.8333333333339644, 1,
0.6666666666680346, 1,
0.3333333333347061, 1,
0.166666666667353, 1,
0, 0.5000000000020591,
0, 0.8333333333339644,
0, 0.6666666666680346,
0, 0.3333333333347061,
0, 0.166666666667353,
0.5, 0.5,
0.2500000000010295, 0.7500000000010295,
0.7500000000005148, 0.7499999999996735,
0.749999999999432, 0.249999999999432,
0.2499999999989704, 0.2500000000010296,
0.2499999999996568, 0.4166666666676962,
0.2500000000003432, 0.5833333333343629,
0.166666666667353, 0.6666666666680394,
0.08333333333367651, 0.5833333333350492,
0.08333333333299014, 0.4166666666683826,
0.1666666666659803, 0.3333333333347061,
0.1666666666666666, 0.5000000000013726,
0.4166666666670099, 0.5833333333336765,
0.3333333333340197, 0.666666666667353,
0.3333333333326469, 0.3333333333340197,
0.4166666666663235, 0.4166666666670099,
0.3333333333333331, 0.5000000000006861,
0.0833333333336765, 0.9166666666670098,
0.166666666667353, 0.8333333333340197,
0.3333333333347061, 0.8333333333340197,
0.4166666666683826, 0.9166666666670098,
0.2500000000010294, 0.9166666666670094,
0.08333333333299014, 0.08333333333367651,
0.1666666666659803, 0.166666666667353,
0.08333333333299012, 0.2500000000010294,
0.9166666666668383, 0.9166666666665578,
0.8333333333336765, 0.8333333333331157,
0.8333333333336765, 0.6666666666660137,
0.9166666666668383, 0.5833333333323538,
0.9166666666668377, 0.7499999999994553,
0.9166666666664773, 0.08333333333314398,
0.8333333333329547, 0.166666666666288,
0.6666666666658527, 0.166666666666288,
0.5833333333322733, 0.08333333333314399,
0.7499999999993748, 0.08333333333314397,
0.5000000000006863, 0.6666666666666666,
0.5000000000013728, 0.8333333333333333,
0.416666666667696, 0.7500000000003427,
0.5833333333335049, 0.5833333333332245,
0.6666666666670098, 0.6666666666664489,
0.6666666666676963, 0.8333333333331157,
0.5833333333348777, 0.9166666666665578,
0.5833333333341911, 0.7499999999998906,
0.6666666666666666, 0.4999999999995647,
0.8333333333333333, 0.4999999999991293,
0.7500000000001711, 0.5833333333327889,
0.4999999999991293, 0.1666666666666667,
0.4999999999995647, 0.3333333333333333,
0.3333333333322116, 0.1666666666673531,
0.4166666666654528, 0.08333333333367654,
0.4166666666658879, 0.2500000000003432,
0.666666666666288 ,0.3333333333329546,
0.5833333333331441, 0.4166666666664773,
0.5833333333327084, 0.2499999999998106,
0.8333333333329547, 0.3333333333325193,
0.9166666666664773, 0.4166666666656066,
0.7499999999998102, 0.4166666666660417,
0.2499999999992977, 0.08333333333367651,
0.9166666666664769, 0.2499999999994515,
0.08333333333367648, 0.7500000000009838,
0.7500000000008122, 0.9166666666665574
  });

  std::vector<int> connectOrd3({
29, 26, 20, 30, 31, 32, 33, 34, 35, 36,
25, 26, 29, 37, 38, 31, 30, 39, 40, 41,
4, 26, 15, 42, 43, 44, 45, 18, 19, 46,
1, 29, 20, 47, 48, 35, 34, 23, 24, 49,
3, 27, 10, 50, 51, 52, 53, 13, 14, 54,
2, 28, 5, 55, 56, 57, 58, 8, 9, 59,
15, 26, 25, 45, 44, 38, 37, 60, 61, 62,
15, 25, 27, 61, 60, 63, 64, 65, 66, 67,
10, 27, 25, 53, 52, 64, 63, 68, 69, 70,
5, 25, 29, 71, 72, 40, 39, 73, 74, 75,
5, 28, 25, 58, 57, 76, 77, 72, 71, 78,
10, 25, 28, 69, 68, 77, 76, 79, 80, 81,
1, 5, 29, 6, 7, 74, 73, 48, 47, 82,
2, 10, 28, 11, 12, 80, 79, 56, 55, 83,
4, 20, 26, 21, 22, 33, 32, 43, 42, 84,
3, 15, 27, 16, 17, 66, 65, 51, 50, 85
      });


  SECTION("Order 3 tests"){
    Mesh m(2, 3, "simplex", 2, nodesOrd3, connectOrd3);
    CHECK(m.getNumberPoints() == nodesOrd3.size());
    CHECK(m.getNumberCells() == connectOrd3.size());
    for(int i = 0; i < nodesOrd3.size(); i++){
      CHECK((*m.getPoints())[i] == nodesOrd3[i]);
    }
    for(int i = 0; i < connectOrd3.size(); i++){
      CHECK((*m.getCells())[i] == connectOrd3[i]);
    }
  };
};
