#include <catch2/catch.hpp>
#include <iostream>
#include <cstdio>
#include "configRessources.h"
#include "HDF5Io.h"

using namespace hfox;

TEST_CASE("Testing HDF5Io", "[unit][io][HDF5Io]"){
  std::string meshDirPath = RessourceConfig::getRessourcePath() + "/meshes/";
  SECTION("Testing construction"){
    CHECK_NOTHROW(HDF5Io());
    Mesh mesh(1, 1, "simplex");
    CHECK_NOTHROW(HDF5Io(&mesh));
    Mesh empty;
    CHECK_THROWS(HDF5Io(&empty));
  };


  Io * hdfio = new HDF5Io();
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
  SECTION("Test load mesh file"){
    std::string inputFile = meshDirPath + "lightTri.h5";
    CHECK_THROWS(hdfio->load(inputFile));
    Mesh simplexMesh(2, 1, "simplex");
    hdfio->setMesh(&simplexMesh);
    CHECK_NOTHROW(hdfio->load(inputFile));
    hdfio->load(inputFile);
    CHECK(simplexMesh.getDimension() == 2);
    CHECK(simplexMesh.getNumberPoints() == testNodes.size());
    CHECK(simplexMesh.getNumberCells() == testSimplexConnectivity.size());
    const std::vector<double> * point;
    for(int i = 0; i < simplexMesh.getNumberPoints(); i++){
      point = simplexMesh.getPoint(i);
      for(int j = 0; j < point->size(); j++){
        CHECK((*point)[j] == testNodes[i][j]);
      }
    }
    const std::vector<int> * cell;
    for(int i = 0; i < simplexMesh.getNumberCells(); i++){
      cell = simplexMesh.getCell(i);
      for(int j = 0; j < cell->size(); j++){
        CHECK((*cell)[j] == testSimplexConnectivity[i][j]);
      }
    }
    inputFile = meshDirPath + "lightQuad.h5";
    Mesh orthoMesh(2, 1, "orthotope");
    hdfio->setMesh(&orthoMesh);
    CHECK_NOTHROW(hdfio->load(inputFile));
    hdfio->load(inputFile);
    CHECK(orthoMesh.getDimension() == 2);
    CHECK(orthoMesh.getNumberPoints() == testNodes.size());
    CHECK(orthoMesh.getNumberCells() == testOrthotopeConnectivity.size());
    for(int i = 0; i < orthoMesh.getNumberPoints(); i++){
      point = orthoMesh.getPoint(i);
      for(int j = 0; j < point->size(); j++){
        CHECK((*point)[j] == testNodes[i][j]);
      }
    }
    for(int i = 0; i < orthoMesh.getNumberCells(); i++){
      cell = orthoMesh.getCell(i);
      for(int j = 0; j < cell->size(); j++){
        CHECK((*cell)[j] == testSimplexConnectivity[i][j]);
      }
    }
  };

  SECTION("Test write mesh file"){
    Mesh simplexMesh(2, 1, "simplex", testNodes, testSimplexConnectivity);
    hdfio->setMesh(&simplexMesh);
    std::string checkFile = meshDirPath + "lightTri.h5";
    std::string outputFile = meshDirPath + "tmp.h5";
    hdfio->write(outputFile);
    CHECK(RessourceConfig::checkFilesEqual(checkFile, outputFile));
    std::remove(outputFile.c_str());
    Mesh orthoMesh(2, 1, "orthotope", testNodes, testSimplexConnectivity);
    hdfio->setMesh(&orthoMesh);
    checkFile = meshDirPath + "lightQuad.h5";
    outputFile = meshDirPath + "tmp.h5";
    hdfio->write(outputFile);
    CHECK(RessourceConfig::checkFilesEqual(checkFile, outputFile));
    std::remove(outputFile.c_str());
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
  SECTION("Test load mesh file (order 2)"){
    std::string inputFile = meshDirPath + "lightTri2.h5";
    CHECK_THROWS(hdfio->load(inputFile));
    Mesh simplexMesh(2, 2, "simplex");
    hdfio->setMesh(&simplexMesh);
    CHECK_NOTHROW(hdfio->load(inputFile));
    hdfio->load(inputFile);
    CHECK(simplexMesh.getDimension() == 2);
    CHECK(simplexMesh.getNumberPoints() == testNodes.size());
    CHECK(simplexMesh.getNumberCells() == testSimplexConnectivity.size());
    const std::vector<double> * point;
    for(int i = 0; i < simplexMesh.getNumberPoints(); i++){
      point = simplexMesh.getPoint(i);
      for(int j = 0; j < point->size(); j++){
        CHECK((*point)[j] == testNodes[i][j]);
      }
    }
    const std::vector<int> * cell;
    for(int i = 0; i < simplexMesh.getNumberCells(); i++){
      cell = simplexMesh.getCell(i);
      for(int j = 0; j < cell->size(); j++){
        CHECK((*cell)[j] == testSimplexConnectivity[i][j]);
      }
    }
    inputFile = meshDirPath + "lightQuad2.h5";
    Mesh orthoMesh(2, 2, "orthotope");
    hdfio->setMesh(&orthoMesh);
    CHECK_NOTHROW(hdfio->load(inputFile));
    hdfio->load(inputFile);
    CHECK(orthoMesh.getDimension() == 2);
    CHECK(orthoMesh.getNumberPoints() == testNodes.size());
    CHECK(orthoMesh.getNumberCells() == testOrthotopeConnectivity.size());
    for(int i = 0; i < orthoMesh.getNumberPoints(); i++){
      point = orthoMesh.getPoint(i);
      for(int j = 0; j < point->size(); j++){
        CHECK((*point)[j] == testNodes[i][j]);
      }
    }
    for(int i = 0; i < orthoMesh.getNumberCells(); i++){
      cell = orthoMesh.getCell(i);
      for(int j = 0; j < cell->size(); j++){
        CHECK((*cell)[j] == testSimplexConnectivity[i][j]);
      }
    }
  };

  SECTION("Test write mesh file (order 2)"){
    Mesh simplexMesh(2, 2, "simplex", testNodes, testSimplexConnectivity);
    hdfio->setMesh(&simplexMesh);
    std::string checkFile = meshDirPath + "lightTri2.h5";
    std::string outputFile = meshDirPath + "tmp.h5";
    hdfio->write(outputFile);
    CHECK(RessourceConfig::checkFilesEqual(checkFile, outputFile));
    std::remove(outputFile.c_str());
    Mesh orthoMesh(2, 2, "orthotope", testNodes, testSimplexConnectivity);
    hdfio->setMesh(&orthoMesh);
    checkFile = meshDirPath + "lightQuad2.h5";
    outputFile = meshDirPath + "tmp.h5";
    hdfio->write(outputFile);
    CHECK(RessourceConfig::checkFilesEqual(checkFile, outputFile));
    std::remove(outputFile.c_str());
  };

  delete hdfio;
};
