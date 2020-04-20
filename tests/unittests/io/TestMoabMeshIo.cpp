#include <catch2/catch.hpp>
#include <iostream>
#include <cstdio>
#include "TestUtils.h"
#include "MoabMeshIo.h"

using namespace hfox;

TEST_CASE("Unit testing the MoabMeshIo class", "[unit][io][MoabMeshIo]"){
  std::string meshDirPath = TestUtils::getRessourcePath() + "/meshes/";
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
  std::vector<double> testNodesStream = TestUtils::unpack(testNodes);
  std::vector<int> testSimplexConnectivityStream = TestUtils::unpack(testSimplexConnectivity);
  std::vector<int> testOrthotopeConnectivityStream = TestUtils::unpack(testOrthotopeConnectivity);
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

  //SECTION("Testing write VTK"){
    //Mesh simplexMesh(2, 1, "simplex", 2, testNodesStream, testSimplexConnectivityStream);
    //moabio->setMesh(&simplexMesh);
    //std::string checkFile = meshDirPath + "generatedLightTri.vtk";
    //std::string outputFile = meshDirPath + "tmpoutput.vtk";
    //moabio->write(outputFile);
    //CHECK(TestUtils::checkFilesEqual(checkFile, outputFile));
    //std::remove(outputFile.c_str());
    //Mesh orthoMesh(2, 1, "orthotope", 2, testNodesStream, testOrthotopeConnectivityStream);
    //moabio->setMesh(&orthoMesh);
    //checkFile = meshDirPath + "generatedLightQuad.vtk";
    //outputFile = meshDirPath + "test.vtk";
    //moabio->write(outputFile);
    //CHECK(TestUtils::checkFilesEqual(checkFile, outputFile));
    //std::remove(outputFile.c_str());
  //};

  SECTION("Testing load H5M"){
    std::string lightTriH5M = "generatedLightTri.h5m";
    Mesh simplexMesh;
    simplexMesh.setReferenceElement(2, 1, "simplex");
    moabio->setMesh(&simplexMesh);
    CHECK_NOTHROW(moabio->load(meshDirPath + lightTriH5M));
    moabio->load(meshDirPath + lightTriH5M);
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
    std::string lightQuadH5M = "generatedLightQuad.h5m";
    Mesh orthoMesh;
    orthoMesh.setReferenceElement(2, 1, "orthotope");
    moabio->setMesh(&orthoMesh);
    moabio->load(meshDirPath + lightQuadH5M);
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

  SECTION("Testing write H5M"){
    Mesh simplexMesh(2, 1, "simplex", 2, testNodesStream, testSimplexConnectivityStream);
    moabio->setMesh(&simplexMesh);
    std::string outputFile = meshDirPath + "tmpoutput.h5m";
    moabio->write(outputFile);
    moabio->load(outputFile);
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
    std::remove(outputFile.c_str());
    Mesh orthoMesh(2, 1, "orthotope", 2, testNodesStream, testOrthotopeConnectivityStream);
    moabio->setMesh(&orthoMesh);
    outputFile = meshDirPath + "tmpoutput.h5m";
    moabio->write(outputFile);
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
  testSimplexConnectivityStream = TestUtils::unpack(testSimplexConnectivity);
  testOrthotopeConnectivityStream = TestUtils::unpack(testOrthotopeConnectivity);
  SECTION("Testing higher order load VTK"){
    std::string lightTriVtk = "generatedLightTri2.vtk";
    CHECK_THROWS(moabio->load(meshDirPath + lightTriVtk));
    Mesh simplexMesh;
    simplexMesh.setReferenceElement(2, 2, "simplex");
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
    std::string lightQuadVtk = "generatedLightQuad2.vtk";
    Mesh orthoMesh;
    orthoMesh.setReferenceElement(2, 2, "orthotope");
    moabio->setMesh(&orthoMesh);
    moabio->load(meshDirPath + lightQuadVtk);
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

  //SECTION("Testing higher order write VTK"){
    //Mesh simplexMesh(2, 2, "simplex", 2, testNodesStream, testSimplexConnectivityStream);
    //moabio->setMesh(&simplexMesh);
    //std::string checkFile = meshDirPath + "generatedLightTri2.vtk";
    //std::string outputFile = meshDirPath + "tmpoutput.vtk";
    //moabio->write(outputFile);
    //CHECK(TestUtils::checkFilesEqual(checkFile, outputFile));
    //std::remove(outputFile.c_str());
    //Mesh orthoMesh(2, 2, "orthotope", 2, testNodesStream, testOrthotopeConnectivityStream);
    //moabio->setMesh(&orthoMesh);
    //checkFile = meshDirPath + "generatedLightQuad2.vtk";
    //outputFile = meshDirPath + "tmpoutput.vtk";
    //moabio->write(outputFile);
    //CHECK(TestUtils::checkFilesEqual(checkFile, outputFile));
    //std::remove(outputFile.c_str());
  //};

  SECTION("Testing higher order load H5M"){
    std::string lightTriH5M = "generatedLightTri2.h5m";
    Mesh simplexMesh;
    simplexMesh.setReferenceElement(2, 2, "simplex");
    moabio->setMesh(&simplexMesh);
    CHECK_NOTHROW(moabio->load(meshDirPath + lightTriH5M));
    moabio->load(meshDirPath + lightTriH5M);
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
    std::string lightQuadH5M = "generatedLightQuad2.h5m";
    Mesh orthoMesh;
    orthoMesh.setReferenceElement(2, 2, "orthotope");
    moabio->setMesh(&orthoMesh);
    moabio->load(meshDirPath + lightQuadH5M);
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

  SECTION("Testing higher order write H5M"){
    Mesh simplexMesh(2, 2, "simplex", 2, testNodesStream, testSimplexConnectivityStream);
    moabio->setMesh(&simplexMesh);
    std::string outputFile = meshDirPath + "tmpoutput.h5m";
    moabio->write(outputFile);
    moabio->load(outputFile);
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
    std::remove(outputFile.c_str());
    Mesh orthoMesh(2, 2, "orthotope", 2, testNodesStream, testOrthotopeConnectivityStream);
    moabio->setMesh(&orthoMesh);
    outputFile = meshDirPath + "tmpoutput.h5m";
    moabio->write(outputFile);
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
    std::remove(outputFile.c_str());
  };


  std::vector<int> nNodesOrdsGmsh = {13, 41};
  SECTION("Testing gmsh load"){
    // higher order not supported by MOAB
    int nOrders = 2;
    for(int i = 0; i < nOrders; i++){
      int order = i + 1;
      Mesh m;
      m.setReferenceElement(2, order, "simplex");
      moabio->setMesh(&m);
      std::string loadFile = meshDirPath + "lightSquareOrd" + std::to_string(order) + ".msh";
      moabio->load(loadFile);
      CHECK(m.getNumberPoints() == nNodesOrdsGmsh[i]);
      CHECK(m.getNumberCells() == 16);
      CHECK(m.getNumberFaces() == 28);
      std::vector<int> cell;
      m.getCell(0, &cell);
      CHECK(cell.size() == m.getReferenceElement()->getNumNodes());
    }
  };

  delete moabio;
};
