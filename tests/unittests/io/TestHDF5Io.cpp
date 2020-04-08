#include <catch2/catch.hpp>
#include <iostream>
#include <cstdio>
#include <algorithm>
#include "TestUtils.h"
#include "HDF5Io.h"
#include "Mesh.h"
#include "Field.h"

using namespace hfox;

TEST_CASE("Testing HDF5Io", "[unit][io][HDF5Io]"){
  std::string meshDirPath = TestUtils::getRessourcePath() + "/meshes/";
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
  bufferElement[0] = 5; bufferElement[1] = 1; bufferElement[2] = 6; bufferElement[3] = 4;
  testOrthotopeConnectivity.push_back(bufferElement);
  bufferElement[0] = 4; bufferElement[1] = 6; bufferElement[2] = 2; bufferElement[3] = 7;
  testOrthotopeConnectivity.push_back(bufferElement);
  bufferElement[0] = 8; bufferElement[1] = 4; bufferElement[2] = 7; bufferElement[3] = 3;
  testOrthotopeConnectivity.push_back(bufferElement);

  std::vector<double> testNodesStream = TestUtils::unpack(testNodes);
  std::vector<int> testSimplexConnectivityStream = TestUtils::unpack(testSimplexConnectivity);
  std::vector<int> testOrthotopeConnectivityStream = TestUtils::unpack(testOrthotopeConnectivity);
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
    std::vector<double> point;
    for(int i = 0; i < simplexMesh.getNumberPoints(); i++){
      simplexMesh.getPoint(i, &point);
      for(int j = 0; j < point.size(); j++){
        CHECK((point)[j] == testNodes[i][j]);
      }
    }
    std::vector<int> cell;
    for(int i = 0; i < simplexMesh.getNumberCells(); i++){
      simplexMesh.getCell(i, &cell);
      for(int j = 0; j < cell.size(); j++){
        CHECK((cell)[j] == testSimplexConnectivity[i][j]);
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
      orthoMesh.getPoint(i, &point);
      for(int j = 0; j < point.size(); j++){
        CHECK((point)[j] == testNodes[i][j]);
      }
    }
    for(int i = 0; i < orthoMesh.getNumberCells(); i++){
      orthoMesh.getCell(i, &cell);
      for(int j = 0; j < cell.size(); j++){
        CHECK((cell)[j] == testOrthotopeConnectivity[i][j]);
      }
    }
  };

  SECTION("Test write mesh file"){
    Mesh simplexMesh(2, 1, "simplex", 2, testNodesStream, testSimplexConnectivityStream);
    hdfio->setMesh(&simplexMesh);
    std::string outputFile = meshDirPath + "tmp.h5";
    hdfio->write(outputFile);
    hdfio->load(outputFile);
    std::vector<double> point;
    for(int i = 0; i < simplexMesh.getNumberPoints(); i++){
      simplexMesh.getPoint(i, &point);
      for(int j = 0; j < point.size(); j++){
        CHECK((point)[j] == testNodes[i][j]);
      }
    }
    std::vector<int> cell;
    for(int i = 0; i < simplexMesh.getNumberCells(); i++){
      simplexMesh.getCell(i, &cell);
      for(int j = 0; j < cell.size(); j++){
        CHECK((cell)[j] == testSimplexConnectivity[i][j]);
      }
    }
    std::remove(outputFile.c_str());
    Mesh orthoMesh(2, 1, "orthotope", 2, testNodesStream, testOrthotopeConnectivityStream);
    hdfio->setMesh(&orthoMesh);
    outputFile = meshDirPath + "tmp.h5";
    hdfio->write(outputFile);
    hdfio->load(outputFile);
    for(int i = 0; i < orthoMesh.getNumberPoints(); i++){
      orthoMesh.getPoint(i, &point);
      for(int j = 0; j < point.size(); j++){
        CHECK((point)[j] == testNodes[i][j]);
      }
    }
    for(int i = 0; i < orthoMesh.getNumberCells(); i++){
      orthoMesh.getCell(i, &cell);
      for(int j = 0; j < cell.size(); j++){
        CHECK((cell)[j] == testOrthotopeConnectivity[i][j]);
      }
    }
    std::remove(outputFile.c_str());
  };


  std::vector<double> nodeFieldData(9);
  std::iota(nodeFieldData.begin(), nodeFieldData.end(), 0);
  std::vector<double> cellFieldData(8*2);
  for(int i = 0; i < 8; i++){
    cellFieldData[i*2] = i;
    cellFieldData[i*2 + 1] = -i;
  }
  SECTION("Load test field"){
    Mesh simplexMesh(2, 1, "simplex", 2, testNodesStream, testSimplexConnectivityStream);
    Io * tmpIo = new HDF5Io(&simplexMesh);
    Field nodeField(&simplexMesh, Node, 1, 1);
    tmpIo->setField("NodeField", &nodeField);
    Field cellField(&simplexMesh, Cell, 2, 1);
    tmpIo->setField("CellField", &cellField);
    std::string loadFile = meshDirPath + "fieldTest.h5";
    try{
    tmpIo->load(loadFile);
    } catch(std::exception){
      std::cout << "there was an error in load" << std::endl;
    }
    for(int i = 0; i < nodeFieldData.size(); i++){
      CHECK((*(nodeField.getValues()))[i] == nodeFieldData[i]);
    }
    for(int i = 0; i < cellFieldData.size(); i++){
      CHECK((*(cellField.getValues()))[i] == cellFieldData[i]);
    }
    delete tmpIo;
  };

  SECTION("Write test field"){
    Mesh simplexMesh(2, 1, "simplex", 2, testNodesStream, testSimplexConnectivityStream);
    Io * tmpIo = new HDF5Io(&simplexMesh);
    Field nodeField(&simplexMesh, Node, 1, 1);
    *(nodeField.getValues()) = nodeFieldData;
    tmpIo->setField("NodeField", &nodeField);
    Field cellField(&simplexMesh, Cell, 2, 1);
    *(cellField.getValues()) = cellFieldData;
    tmpIo->setField("CellField", &cellField);
    std::string outputFile = meshDirPath + "tmp.h5";
    try{
    tmpIo->write(outputFile);
    tmpIo->load(outputFile);
    }catch(std::exception){
      std::cout << "there was an error in io" << std::endl;
    }
    for(int i = 0; i < nodeFieldData.size(); i++){
      CHECK((*(nodeField.getValues()))[i] == nodeFieldData[i]);
    }
    for(int i = 0; i < cellFieldData.size(); i++){
      CHECK((*(cellField.getValues()))[i] == cellFieldData[i]);
    }
    std::remove(outputFile.c_str());
    delete tmpIo;
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
  testSimplexConnectivityStream = TestUtils::unpack(testSimplexConnectivity);
  testOrthotopeConnectivityStream = TestUtils::unpack(testOrthotopeConnectivity);
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
    std::vector<double> point;
    for(int i = 0; i < simplexMesh.getNumberPoints(); i++){
      simplexMesh.getPoint(i, &point);
      for(int j = 0; j < point.size(); j++){
        CHECK((point)[j] == testNodes[i][j]);
      }
    }
    std::vector<int> cell;
    for(int i = 0; i < simplexMesh.getNumberCells(); i++){
      simplexMesh.getCell(i, &cell);
      for(int j = 0; j < cell.size(); j++){
        CHECK((cell)[j] == testSimplexConnectivity[i][j]);
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
      orthoMesh.getPoint(i, &point);
      for(int j = 0; j < point.size(); j++){
        CHECK((point)[j] == testNodes[i][j]);
      }
    }
    for(int i = 0; i < orthoMesh.getNumberCells(); i++){
      orthoMesh.getCell(i, &cell);
      for(int j = 0; j < cell.size(); j++){
        CHECK((cell)[j] == testOrthotopeConnectivity[i][j]);
      }
    }
  };

  SECTION("Test write mesh file (order 2)"){
    Mesh simplexMesh(2, 2, "simplex", 2, testNodesStream, testSimplexConnectivityStream);
    hdfio->setMesh(&simplexMesh);
    std::string outputFile = meshDirPath + "tmp.h5";
    hdfio->write(outputFile);
    hdfio->load(outputFile);
    std::vector<double> point;
    for(int i = 0; i < simplexMesh.getNumberPoints(); i++){
      simplexMesh.getPoint(i, &point);
      for(int j = 0; j < point.size(); j++){
        CHECK((point)[j] == testNodes[i][j]);
      }
    }
    std::vector<int> cell;
    for(int i = 0; i < simplexMesh.getNumberCells(); i++){
      simplexMesh.getCell(i, &cell);
      for(int j = 0; j < cell.size(); j++){
        CHECK((cell)[j] == testSimplexConnectivity[i][j]);
      }
    }
    std::remove(outputFile.c_str());
    Mesh orthoMesh(2, 2, "orthotope", 2, testNodesStream, testOrthotopeConnectivityStream);
    hdfio->setMesh(&orthoMesh);
    outputFile = meshDirPath + "tmp.h5";
    hdfio->write(outputFile);
    hdfio->load(outputFile);
    for(int i = 0; i < orthoMesh.getNumberPoints(); i++){
      orthoMesh.getPoint(i, &point);
      for(int j = 0; j < point.size(); j++){
        CHECK((point)[j] == testNodes[i][j]);
      }
    }
    for(int i = 0; i < orthoMesh.getNumberCells(); i++){
      orthoMesh.getCell(i, &cell);
      for(int j = 0; j < cell.size(); j++){
        CHECK((cell)[j] == testOrthotopeConnectivity[i][j]);
      }
    }
    std::remove(outputFile.c_str());
  };

  delete hdfio;
};
