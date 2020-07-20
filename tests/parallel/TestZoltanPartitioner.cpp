#include <catch2/catch.hpp>
#include "ZoltanPartitioner.h"
#include "Mesh.h"
#include "HDF5Io.h"
#include "TestUtils.h"

using namespace hfox;

TEST_CASE("Testing ZoltanPartitioner class", "[par][unit][parallel][ZoltanPartitioner]"){

  SECTION("Testing setup"){
    Mesh seqMesh(2, 3, "simplex");
    CHECK_NOTHROW(ZoltanPartitioner(&seqMesh));
    ZoltanOpts opts;
    CHECK_NOTHROW(ZoltanPartitioner(&seqMesh, opts));
    ZoltanPartitioner zPart(&seqMesh);
    CHECK_NOTHROW(zPart.initialize());
  };

  SECTION("Test initialization"){
    Mesh seqMesh(2, 3, "simplex");
    ZoltanPartitioner zPart(&seqMesh);
    zPart.initialize();
    int worldSize;
    MPI_Comm_size(MPI_COMM_WORLD, &worldSize);
    CHECK(zPart.getNumPartitions() == worldSize);
    int worldRank;
    MPI_Comm_rank(MPI_COMM_WORLD, &worldRank);
    CHECK(zPart.getRank() == worldRank);
  };

  SECTION("Test a mesh partition"){
    std::string meshLocation = TestUtils::getRessourcePath() + "/meshes/regression/regression_dim-2_h-5e-2_ord-3.h5";
    Mesh seqMesh(2, 3, "simplex");
    HDF5Io hdf5io(&seqMesh);
    hdf5io.load(meshLocation);
    Mesh parMesh(2, 3, "simplex");
    HDF5Io parhdf5io(&parMesh);
    parhdf5io.load(meshLocation);
    ZoltanOpts myOpts;
    myOpts.debugLevel = "1";
    ZoltanPartitioner zPart(&parMesh, myOpts);
    zPart.initialize();
    zPart.computePartition();
    zPart.update();
    int totNodes = zPart.getTotalNumberNodes();
    int seqNodes = seqMesh.getNumberPoints();
    MPI_Bcast(&seqNodes, 1, MPI_INT, 0, MPI_COMM_WORLD);
    CHECK(totNodes == seqNodes);
    int totCells = zPart.getTotalNumberEls();
    int seqCells = seqMesh.getNumberCells();
    MPI_Bcast(&seqCells, 1, MPI_INT, 0, MPI_COMM_WORLD);
    CHECK(totCells == seqCells);
    int totFaces = zPart.getTotalNumberFaces();
    int seqFaces = seqMesh.getNumberFaces();
    MPI_Bcast(&seqFaces, 1, MPI_INT, 0, MPI_COMM_WORLD);
    CHECK(totFaces == seqFaces);
    std::vector<double> nodes(totNodes*2);
    std::copy(seqMesh.getPoints()->begin(), seqMesh.getPoints()->end(), nodes.begin());
    MPI_Bcast(nodes.data(), nodes.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    std::vector<int> cells(totCells*(parMesh.getReferenceElement()->getNumNodes()));
    std::copy(seqMesh.getCells()->begin(), seqMesh.getCells()->end(), cells.begin());
    MPI_Bcast(cells.data(), cells.size(), MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> faces(totFaces*(parMesh.getReferenceElement()->getFaceElement()->getNumNodes()));
    std::copy(seqMesh.getFaces()->begin(), seqMesh.getFaces()->end(), faces.begin());
    MPI_Bcast(faces.data(), faces.size(), MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> face2CellMap(totFaces*2);
    std::copy(seqMesh.getFace2CellMap()->begin(), seqMesh.getFace2CellMap()->end(), face2CellMap.begin());
    MPI_Bcast(face2CellMap.data(), face2CellMap.size(), MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> cell2FaceMap(totCells*3);
    std::copy(seqMesh.getCell2FaceMap()->begin(), seqMesh.getCell2FaceMap()->end(), cell2FaceMap.begin());
    MPI_Bcast(cell2FaceMap.data(), cell2FaceMap.size(), MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<double> point;
    std::vector<int> cell;
    //local checking
    int globIndex;
    int nNodesPCell = parMesh.getReferenceElement()->getNumNodes();
    int nNodesPFace = parMesh.getReferenceElement()->getFaceElement()->getNumNodes();
    for(int i = 0; i < parMesh.getNumberPoints(); i++){
      parMesh.getPoint(i, &point);
      globIndex = zPart.local2GlobalNode(i);
      for(int k = 0; k < 2; k++){
        CHECK(point[k] == nodes[globIndex*2 + k]);
      }
    }
    for(int i = 0; i < parMesh.getNumberCells(); i++){
      parMesh.getCell(i, &cell);
      globIndex = zPart.local2GlobalEl(i);
      for(int k = 0; k < nNodesPCell; k++){
        CHECK(cell[k] == cells[globIndex*nNodesPCell + k]);
      }
    }
    for(int i = 0; i < parMesh.getNumberFaces(); i++){
      parMesh.getFace(i, &cell);
      globIndex = zPart.local2GlobalFace(i);
      for(int k = 0; k < nNodesPFace; k++){
        CHECK(cell[k] == faces[globIndex*nNodesPFace + k]);
      }
    }
    for(int i = 0; i < parMesh.getNumberFaces(); i++){
      parMesh.getFace2Cell(i, &cell);
      globIndex = zPart.local2GlobalFace(i);
      for(int k = 0; k < cell.size(); k++){
        CHECK(cell[k] == face2CellMap[globIndex*2 + k]);
      }
    }
    for(int i = 0; i < parMesh.getNumberCells(); i++){
      parMesh.getCell2Face(i, &cell);
      globIndex = zPart.local2GlobalEl(i);
      for(int k = 0; k < 3; k++){
        CHECK(cell[k] == cell2FaceMap[globIndex*3 + k]);
      }
    }
    //global checking
    int locIndex;
    for(int i = 0; i < nodes.size()/2; i++){
      locIndex = zPart.global2LocalNode(i);
      if(locIndex != -1){
        parMesh.getPoint(locIndex, &point);
        for(int k = 0; k < 2; k++){
          CHECK(point[k] == nodes[i*2 + k]);
        }
      }
    }
    for(int i = 0; i < cells.size()/nNodesPCell; i++){
      locIndex = zPart.global2LocalElement(i);
      if(locIndex != -1){
        parMesh.getCell(locIndex, &cell);
        for(int k = 0; k < nNodesPCell; k++){
          CHECK(cell[k] == cells[i*nNodesPCell + k]);
        }
        parMesh.getCell2Face(locIndex, &cell);
        for(int k = 0; k < 3; k++){
          CHECK(cell[k] == cell2FaceMap[i*3 + k]);
        }
      }
    }
    for(int i = 0; i < faces.size()/nNodesPFace; i++){
      locIndex = zPart.global2LocalFace(i);
      if(locIndex != -1){
        parMesh.getFace(locIndex, &cell);
        for(int k = 0; k < nNodesPFace; k++){
          CHECK(cell[k] == faces[i*nNodesPFace + k]);
        }
        parMesh.getFace2Cell(locIndex, &cell);
        for(int k = 0; k < cell.size(); k++){
          CHECK(cell[k] == face2CellMap[i*2 + k]);
        }
      }
    }
  };

  SECTION("Testing a field partition"){
    std::string meshLocation = TestUtils::getRessourcePath() + "/meshes/regression/regression_dim-2_h-5e-2_ord-3.h5";
    Mesh seqMesh(2, 3, "simplex");
    HDF5Io hdf5io(&seqMesh);
    hdf5io.load(meshLocation);
    Mesh parMesh(2, 3, "simplex");
    HDF5Io parhdf5io(&parMesh);
    parhdf5io.load(meshLocation);
    Field nodeField(&parMesh, Node, 1, 1);
    std::iota(nodeField.getValues()->begin(), nodeField.getValues()->end(), 0);
    Field cellField(&parMesh, Cell, 1, 1);
    std::iota(cellField.getValues()->begin(), cellField.getValues()->end(), 0);
    Field faceField(&parMesh, Face, 1, 1);
    std::iota(faceField.getValues()->begin(), faceField.getValues()->end(), 0);
    std::vector<Field*> fieldList = {&nodeField, &cellField, &faceField};
    ZoltanPartitioner zPart(&parMesh);
    zPart.initialize();
    zPart.setFields(fieldList);
    zPart.computePartition();
    zPart.update();
    int totNodes = zPart.getTotalNumberNodes();
    int totFaces = zPart.getTotalNumberFaces();
    int totCells = zPart.getTotalNumberEls();
    CHECK(*(nodeField.getNumEntities()) == parMesh.getNumberPoints());
    CHECK(*(faceField.getNumEntities()) == parMesh.getNumberFaces());
    CHECK(*(cellField.getNumEntities()) == parMesh.getNumberCells());
    CHECK((nodeField.getValues()->size()) == parMesh.getNumberPoints());
    CHECK((faceField.getValues()->size()) == parMesh.getNumberFaces());
    CHECK((cellField.getValues()->size()) == parMesh.getNumberCells());
    int globIndex;
    for(int i = 0; i < parMesh.getNumberPoints(); i++){
      globIndex = zPart.local2GlobalNode(i);
      CHECK(globIndex == nodeField.getValues()->at(i));
    }
    for(int i = 0; i < parMesh.getNumberCells(); i++){
      globIndex = zPart.local2GlobalEl(i);
      CHECK(globIndex == cellField.getValues()->at(i));
    }
    for(int i = 0; i < parMesh.getNumberFaces(); i++){
      globIndex = zPart.local2GlobalFace(i);
      CHECK(globIndex == faceField.getValues()->at(i));
    }
    int locIndex;
    for(int i = 0; i < totNodes; i++){
      locIndex = zPart.global2LocalNode(i);
      if(locIndex != -1){
        CHECK(i == nodeField.getValues()->at(locIndex));
      }
    }
    for(int i = 0; i < totFaces; i++){
      locIndex = zPart.global2LocalFace(i);
      if(locIndex != -1){
        CHECK(i == faceField.getValues()->at(locIndex));
      }
    }
    for(int i = 0; i < totCells; i++){
      locIndex = zPart.global2LocalElement(i);
      if(locIndex != -1){
        CHECK(i == cellField.getValues()->at(locIndex));
      }
    }
  };

  SECTION("Testing shared information"){
    std::string meshLocation = TestUtils::getRessourcePath() + "/meshes/regression/regression_dim-2_h-5e-2_ord-3.h5";
    Mesh seqMesh(2, 3, "simplex");
    HDF5Io hdf5io(&seqMesh);
    hdf5io.load(meshLocation);
    Mesh parMesh(2, 3, "simplex");
    HDF5Io parhdf5io(&parMesh);
    parhdf5io.load(meshLocation);
    Field nodeField(&parMesh, Node, 1, 1);
    std::iota(nodeField.getValues()->begin(), nodeField.getValues()->end(), 0);
    Field cellField(&parMesh, Cell, 1, 1);
    std::iota(cellField.getValues()->begin(), cellField.getValues()->end(), 0);
    Field faceField(&parMesh, Face, 1, 1);
    std::iota(faceField.getValues()->begin(), faceField.getValues()->end(), 0);
    std::vector<Field*> fieldList = {&nodeField, &cellField, &faceField};
    ZoltanPartitioner zPart(&parMesh);
    zPart.initialize();
    zPart.setFields(fieldList);
    zPart.computePartition();
    zPart.update();
    int totNodes = zPart.getTotalNumberNodes();
    int totFaces = zPart.getTotalNumberFaces();
    int totCells = zPart.getTotalNumberEls();
    std::vector<double> nodes(totNodes*2);
    std::copy(seqMesh.getPoints()->begin(), seqMesh.getPoints()->end(), nodes.begin());
    MPI_Bcast(nodes.data(), nodes.size(), MPI_DOUBLE, 0, MPI_COMM_WORLD);
    std::vector<int> cells(totCells*(parMesh.getReferenceElement()->getNumNodes()));
    std::copy(seqMesh.getCells()->begin(), seqMesh.getCells()->end(), cells.begin());
    MPI_Bcast(cells.data(), cells.size(), MPI_INT, 0, MPI_COMM_WORLD);
    std::vector<int> faces(totFaces*(parMesh.getReferenceElement()->getFaceElement()->getNumNodes()));
    std::copy(seqMesh.getFaces()->begin(), seqMesh.getFaces()->end(), faces.begin());
    MPI_Bcast(faces.data(), faces.size(), MPI_INT, 0, MPI_COMM_WORLD);
    int dimSpace = parMesh.getNodeSpaceDimension();
    int nNodesPFace = parMesh.getReferenceElement()->getFaceElement()->getNumNodes();
    int globIndex, locIndex;
    std::vector<int> cell;
    std::vector<int> cell2Face;
    std::vector<int> face2Cell;
    std::vector<double> point;
    std::vector<double> fVals;
    for(int iEl = 0; iEl < parMesh.getNumberCells(); iEl++){
      parMesh.getCell2Face(iEl, &cell2Face);
      for(int i = 0; i < cell2Face.size(); i++){
        locIndex = zPart.global2LocalFace(cell2Face[i]);
        if(locIndex == -1){
          CHECK_NOTHROW(parMesh.getGhostFace(cell2Face[i], &cell));
          for(int k = 0; k < cell.size(); k++){
            CHECK(cell[k] == faces[cell2Face[i]*nNodesPFace + k]);
            CHECK_NOTHROW(parMesh.getGhostPoint(cell[k], &point));
            for(int j = 0; j < dimSpace; j++){
              CHECK(point[j] == nodes[cell[k]*dimSpace + j]);
            }
          } 
          CHECK_NOTHROW(faceField.getParValues(cell2Face[i], &fVals));
          CHECK(fVals[0] == cell2Face[i]);
          for(int k = 0; k < cell.size(); k++){
            CHECK_NOTHROW(nodeField.getParValues(cell[k], &fVals));
            CHECK(fVals[0] == cell[k]);
          }
        }
      }
    }
  };
};
