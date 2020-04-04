#include <catch2/catch.hpp>
#include <iostream>
#include <cstdio>
#include "configRessources.h"
#include "MoabMeshIo.h"

using namespace hfox;

TEST_CASE("Unit testing the MoabMeshIo class", "[unit][io][MoabMeshIo]"){
  std::string meshDirPath = getRessourcePath() + "/meshes/";
  SECTION("Testing construction"){
    CHECK_NOTHROW(MoabMeshIo());
    Mesh refElemMesh;
    refElemMesh.setReferenceElement(1, 1, "simplex");
    CHECK_NOTHROW(MoabMeshIo(&refElemMesh));
    Mesh emptyMesh;
    CHECK_THROWS(MoabMeshIo(&emptyMesh));
  };

  MoabMeshIo * moabio;
  moabio = new MoabMeshIo();

  SECTION("Testing setMesh"){
    Mesh refElemMesh;
    refElemMesh.setReferenceElement(1, 1, "simplex");
    CHECK_NOTHROW(moabio->setMesh(&refElemMesh));
    Mesh emptyMesh;
    CHECK_THROWS(moabio->setMesh(&emptyMesh));
  };

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
  SECTION("Testing load VTK"){
    std::string lightTriVtk = "lightTestTri.vtk";
    CHECK_THROWS(moabio->load(meshDirPath + lightTriVtk));
    Mesh simplexMesh;
    simplexMesh.setReferenceElement(2, 1, "simplex");
    moabio->setMesh(&simplexMesh);
    CHECK_NOTHROW(moabio->load(meshDirPath + lightTriVtk));
    moabio->load(meshDirPath + lightTriVtk);
    CHECK(simplexMesh.getDimension() == 2);
    CHECK(simplexMesh.getNumberPoints() == testNodes.size());
    CHECK(simplexMesh.getNumberCells() == testSimplexConnectivity.size());
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
    std::string lightQuadVtk = "lightTestQuad.vtk";
    Mesh orthoMesh;
    orthoMesh.setReferenceElement(2, 1, "orthotope");
    moabio->setMesh(&orthoMesh);
    moabio->load(meshDirPath + lightQuadVtk);
    CHECK(orthoMesh.getDimension() == 2);
    CHECK(orthoMesh.getNumberPoints() == testNodes.size());
    CHECK(orthoMesh.getNumberCells() == testOrthotopeConnectivity.size());
    for(int i = 0; i < orthoMesh.getNumberPoints(); i++){
      orthoMesh.getPoint(i, &point);
      for(int j = 0; j < point.size(); j++){
        CHECK(point[j] == testNodes[i][j]);
      }
    }
    for(int i = 0; i < orthoMesh.getNumberCells(); i++){
      orthoMesh.getCell(i, &cell);
      for(int j = 0; j < cell.size(); j++){
        CHECK(cell[j] == testOrthotopeConnectivity[i][j]);
      }
    }
  };

  SECTION("Testing write VTK"){
    Mesh simplexMesh(2, 1, "simplex", testNodes, testSimplexConnectivity);
    moabio->setMesh(&simplexMesh);
    std::string checkFile = meshDirPath + "generatedLightTri.vtk";
    std::string outputFile = meshDirPath + "tmpoutput.vtk";
    moabio->write(outputFile);
    CHECK(checkFilesEqual(checkFile, outputFile));
    std::remove(outputFile.c_str());
    Mesh orthoMesh(2, 1, "orthotope", testNodes, testOrthotopeConnectivity);
    moabio->setMesh(&orthoMesh);
    checkFile = meshDirPath + "generatedLightQuad.vtk";
    outputFile = meshDirPath + "tmpoutput.vtk";
    moabio->write(outputFile);
    CHECK(checkFilesEqual(checkFile, outputFile));
    std::remove(outputFile.c_str());
  };
    Mesh simplexMesh(2, 1, "simplex", testNodes, testSimplexConnectivity);
    moabio->setMesh(&simplexMesh);
    std::string outputFile = meshDirPath + "generatedLightTri.h5m";
    moabio->write(outputFile);
    moabio->load(outputFile);
    std::cout << "Num Nodes: " << simplexMesh.getNumberPoints() << std::endl;
    std::cout << "Num Cells: " << simplexMesh.getNumberCells() << std::endl;
    std::cout << "Num Faces: " << simplexMesh.getNumberFaces() << std::endl;

  delete moabio;
};
