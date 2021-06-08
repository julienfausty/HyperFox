#include <string>
#include <vector>
#include <cstdio>
#include <exception>
#include <tuple>
#include <cmath>
#include <moab/Core.hpp>
#include "ErrorHandle.h"
#include "DenseEigen.h"
#include "Mesh.h"
#include "ReferenceElement.h"
#include "HDF5Io.h"

/*!
 * \brief a conversion utility from linear simplicial meshes generated with gmsh to the native h5 format of higher user defined order
 */

using namespace hfox;

std::string usage();
void generateHigherOrderMesh(Mesh * hoMesh, int spaceDim, std::string input);
void readMesh(std::string input, int spaceDim, std::vector<double> * linNodes, std::vector< std::vector<int> > * linCells,
    std::vector< std::vector<int> > * elementAdjacencies);
void generateCellNodes(const ReferenceElement * refElement, std::vector<double> * linCellNodes, 
    std::vector<double> * cellNodes);

int main(int argc, char **argv){
  MPI_Init(&argc, &argv);
  std::string inputFile, outputFile;
  int dimension = 0;
  int spaceDim = 0;
  int order = 0;
  std::string optstr("i:d:s:r:o:vh");
  int opt;
  while((opt = getopt(argc, argv, optstr.c_str())) != -1){
    switch(opt){
      case 'i':
        {
          inputFile.assign(optarg);
          break;
        }
      case 'd':
        {
          dimension = std::stoi(optarg);
          break;
        }
      case 's':
        {
          spaceDim = std::stoi(optarg);
          break;
        }
      case 'r':
        {
          order = std::stoi(optarg);
          break;
        }
      case 'o':
        {
          outputFile.assign(optarg);
          break;
        }
      case 'h':
        {
          std::cout << usage() << std::endl;
          return EXIT_SUCCESS;
          break;
        }
    }
  }
  std::string addUsage = "\n " + usage();

  if(inputFile.empty()){
    std::cout << "convertGmsh2H5HO: main: no inputFile provided" + addUsage << std::endl;
    return EXIT_FAILURE;
  }
  std::string inputExt = inputFile.substr(inputFile.find_last_of('.'), inputFile.size());
  if(inputExt != ".msh"){
    std::cout << "convertGmsh2H5HO: main: inputFile extension must be .msh" + addUsage << std::endl;
    return EXIT_FAILURE;
  }
  if(dimension == 0){
    std::cout << "convertGmsh2H5HO: main: dimension variable cannot be null and must be specified" + addUsage << std::endl;
    return EXIT_FAILURE;
  }
  if(spaceDim == 0){
    spaceDim = dimension;
  }
  if(order == 0){
    std::cout << "convertGmsh2H5HO: main: requiredOrder variable cannot be null and must be specified" + addUsage << std::endl;
    return EXIT_FAILURE;
  }
  if(outputFile.empty()){
    outputFile = inputFile.substr(0, inputFile.find_last_of('.')) + ".h5";
    std::cout << "No output file provided, so setting output to be: " + outputFile << std::endl;
  }
  Mesh hoMesh(dimension, order, "simplex");
  try{
    generateHigherOrderMesh(&hoMesh, spaceDim, inputFile);
  } catch(const std::exception & e){
    std::cout << "The mesh generation process failed with an exception:" << std::endl;
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "New mesh characteristics" << std::endl;
  std::cout << "Dimension: " << hoMesh.getDimension() << std::endl;
  std::cout << "Space dimension: " << hoMesh.getNodeSpaceDimension() << std::endl;
  std::cout << "Number of nodes: " << hoMesh.getNumberPoints() << std::endl;
  std::cout << "Number of cells: " << hoMesh.getNumberCells() << std::endl;
  std::cout << "Number of faces: " << hoMesh.getNumberFaces() << std::endl;
  std::cout << "Number of boundary faces: " << hoMesh.getBoundaryFaces()->size() << std::endl; 
  try{
    HDF5Io hdfio(&hoMesh);
    hdfio.write(outputFile);
  }catch(const std::exception & e){
    std::cout << "The mesh writing process failed with an exception:" << std::endl;
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  MPI_Finalize();
  return EXIT_SUCCESS;
};//main

std::string usage(){
  return std::string("Usage is: convertGmsh2H5HO [-i inputFile] [-d dimension] [-r requiredOrder] [-o outputFile] [-s space dimension] [-h]");
};

void generateHigherOrderMesh(Mesh * hoMesh, int spaceDim, std::string input){
  int dim = hoMesh->getDimension();
  std::vector<double> linNodes;
  std::vector< std::vector<int> > linCells(dim);
  std::vector< std::vector<int> > elementAdjacencies(dim);

  readMesh(input, spaceDim, &linNodes, &linCells, &elementAdjacencies);

  std::cout << "Linear mesh characteristics" << std::endl;
  std::cout << "Number of nodes: " << linNodes.size()/spaceDim << std::endl;
  for(int i = 0; i < dim; i++){
    std::cout << "Number of " + std::to_string(i+1) + " cells: " << linCells[i].size()/(i+2) << std::endl;
  }

  int nCells = linCells[dim-1].size()/(dim+1);
  int nNodesPerEl = hoMesh->getReferenceElement()->getNumNodes();
  std::vector<double> hoNodes;
  std::vector< std::vector<int> > hoCells(dim+1);
  std::vector<int> hoElements(nCells*nNodesPerEl);

  std::vector<int> nCellsPerTopDim(dim);
  std::vector<int> nInnerNodesPerCell(dim);
  std::vector<ReferenceElement *> refElements(dim);
  std::vector< std::vector<bool> > tracking(dim+1);
  tracking[0].resize(linNodes.size()/spaceDim , false);
  hoCells[0].resize(linNodes.size()/spaceDim);
  for(int i = 0; i < dim; i++){
    nCellsPerTopDim[i] = elementAdjacencies[i].size()/nCells;
    refElements[i] = new ReferenceElement(i+1, hoMesh->getReferenceElement()->getOrder(), "simplex");
    tracking[i+1].resize(linCells[i].size()/(i+2), false);
    nInnerNodesPerCell[i] = refElements[i]->getInnerNodes()->size();
    hoCells[i+1].resize(nInnerNodesPerCell[i]*((linCells[i].size())/(i+2)));
  }

  std::vector<double> cellNodes, totalCellNodes, elNodes, linCellNodes, linElNodes((dim+1)*spaceDim);
  std::vector<int> cell;
  std::vector<double> cellNode(spaceDim);
  int cellId, nodeId;
  int hoNodesIndex = 0;
  //big element loop
  for(int iEl = 0; iEl < nCells; iEl++){
    cell.resize(dim+1);
    for(int i = 0; i < dim+1; i++){
      nodeId = linCells[dim-1][iEl*(dim+1) + i];
      for(int j = 0; j < spaceDim; j++){
        linElNodes[i*spaceDim + j] = linNodes[nodeId*spaceDim + j];
      }
      if(!tracking[0][nodeId]){
        for(int j = 0; j < spaceDim ; j++){
          hoNodes.push_back(linElNodes[i*spaceDim + j]);
        }
        tracking[0][nodeId] = true;
        hoCells[0][nodeId] = hoNodesIndex;
        cell[i] = hoNodesIndex;
        hoNodesIndex += 1;
      } else{
        cell[i] = hoCells[0][nodeId];
      }
    }
    //add linNodes
    for(int i = 0; i < cell.size(); i++){
      hoElements[iEl*nNodesPerEl + i] = cell[i];
    }
    generateCellNodes(hoMesh->getReferenceElement(), &linElNodes, &elNodes);
    //topological dim loop
    for(int i = 0; i < dim; i++){
      linCellNodes.resize((i+2)*spaceDim);
      //cells in element loop
      for(int j = 0; j < nCellsPerTopDim[i]; j++){
        cellId = elementAdjacencies[i][iEl*nCellsPerTopDim[i] + j];
        cell.resize(nInnerNodesPerCell[i]);
        //did we already generate the cell or not?
        if(!tracking[i+1][cellId]){
          //generate the nodes and the cell
          for(int l = 0; l < (i+2); l++){
            nodeId = linCells[i][cellId*(i+2) + l];
            for(int k = 0; k < spaceDim ; k++){
              linCellNodes[l*spaceDim + k] = linNodes[nodeId*spaceDim + k];
            }
          }
          generateCellNodes(refElements[i], &linCellNodes, &totalCellNodes);
          //generate the new cell
          std::iota(cell.begin(), cell.end(), hoNodesIndex);
          //put these new nodes in the hoNodes (only inner nodes though)
          const std::vector<int> * innerNodes = refElements[i]->getInnerNodes();
          cellNodes.resize(nInnerNodesPerCell[i]*spaceDim);
          for(int l = 0; l < innerNodes->size(); l++){
            for(int k = 0; k < spaceDim; k++){
              cellNodes[l*spaceDim + k] = totalCellNodes[(*innerNodes)[l]*spaceDim + k];
              hoNodes.push_back(cellNodes[l*spaceDim+k]);
            }
            hoNodesIndex += 1;
          }
          //put that in hoCells
          for(int l = 0; l < cell.size(); l++){
            hoCells[i+1][cellId*nInnerNodesPerCell[i] + l] = cell[l];
          }
          //set tracking to true
          tracking[i+1][cellId] = true;
        }else{
          //get the right cell and corresponding nodes
          cell.assign(hoCells[i+1].begin() + cellId*nInnerNodesPerCell[i], 
              hoCells[i+1].begin() + (cellId+1)*nInnerNodesPerCell[i]);
          //get the cellNodes from hoNodes
          cellNodes.resize(cell.size()*spaceDim);
          for(int l = 0; l < cell.size(); l++){
            for(int k = 0; k < spaceDim; k++){
              cellNodes[l*spaceDim + k] = hoNodes[cell[l]*spaceDim + k];
            }
          }
        }
        //put the cell in hoElements in the right order using elNodes
        bool areEqual;
        for(int l = 0; l < cell.size(); l++){
          for(int k = 0; k < spaceDim; k++){
            cellNode[k] = cellNodes[l*spaceDim + k];
          }
          for(int k = 0; k < nNodesPerEl; k++){
            for(int p = 0; p < spaceDim; p++){
              areEqual = (std::abs((cellNode[p] - elNodes[k*spaceDim + p])) < 1e-8);
              if(!areEqual){
                break;
              }
            }
            if(areEqual){
              hoElements[iEl*nNodesPerEl + k] = cell[l];
              break;
            }
          }
          if(!areEqual){
            throw(ErrorHandle("convertGmsh2H5HO", "generateHigherOrderMesh", "one of the cell nodes could not be found in element"));
          }
        }
      }
    }
  }
  hoMesh->setMesh(spaceDim, hoNodes, hoElements);
  for(int i = 0; i < refElements.size(); i++){
    delete refElements[i];
  }
};

void readMesh(std::string input, int spaceDim, std::vector<double> * linNodes, std::vector< std::vector<int> > * linCells,
    std::vector< std::vector<int> > * elementAdjacencies){
  int dim = linCells->size();
  moab::ErrorCode mbErr;
  moab::Interface * mbIFace = new (std::nothrow) moab::Core;
  if(mbIFace == NULL){
    throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not construct new moab mesh."));
  }
  moab::EntityHandle meshset;
  mbErr = mbIFace->create_meshset(moab::MESHSET_SET, meshset);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not create meshset"));
  }
  mbErr = mbIFace->load_file(input.c_str(), &meshset);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not load mesh file: " + input));
  }
  //fill nodes
  moab::Range vertexes;
  mbErr = mbIFace->get_entities_by_dimension(meshset, 0, vertexes);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not get vertexes from meshset"));
  }
  linNodes->resize(vertexes.size()*spaceDim);
  std::vector<double> vertex(3);
  for(moab::Range::iterator it = vertexes.begin(); it != vertexes.end(); it++){
    mbErr = mbIFace->get_coords(&(*it), 1, vertex.data());
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not get vertex coords from list"));
    }
    int index = mbIFace->id_from_handle(*it) -1;
    for(int i = 0; i < spaceDim ; i++){
      (*linNodes)[index*spaceDim + i] = vertex[i];
    }
  }
  //fill cells
  moab::Range elems;
  mbErr = mbIFace->get_entities_by_dimension(meshset, dim, elems);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not get elements from meshset"));
  }
  for(int i = 0; i < dim; i++){
    moab::Range tmp0;
    mbErr = mbIFace->get_adjacencies(elems, i+1, 1, tmp0, moab::Interface::UNION);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not develop face adjacencies from cells (moab)"));
    }
    mbErr = mbIFace->add_entities(meshset, tmp0);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "computeFaces", "could not add faces to meshset (moab)"));
    }
    moab::Range tmp1;
    mbErr = mbIFace->get_entities_by_dimension(meshset, i+1, tmp1);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("convertGmsh2H5HO", "readMesh", 
            "could not get " + std::to_string(i) + " cells from meshset"));
    }
    int nNodesPerCell = i+2;
    (*linCells)[i].resize(nNodesPerCell*tmp1.size());
    for(moab::Range::iterator it = tmp1.begin(); it != tmp1.end(); it++){
      std::vector<moab::EntityHandle> verts;
      mbErr = mbIFace->get_connectivity(&(*it), 1, verts);
      if(mbErr != moab::MB_SUCCESS){
        throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not get vertices from element"));
      }
      int index = mbIFace->id_from_handle(*it) - 1;
      for(int j = 0; j < nNodesPerCell; j++){
        (*linCells)[i][index*nNodesPerCell + j] = mbIFace->id_from_handle(verts[j])-1;
      }
    }
    moab::Range tmpAdj;
    mbErr = mbIFace->get_adjacencies(&(*(elems.begin())), 1, i+1, 1, tmpAdj);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not " + std::to_string(i+1) + " adjacencies"));
    }
    int nCellsPerElement = tmpAdj.size();
    (*elementAdjacencies)[i].resize(elems.size()*nCellsPerElement);
    for(moab::Range::iterator itEl = elems.begin(); itEl != elems.end(); itEl++){
      moab::Range adj;
      mbErr = mbIFace->get_adjacencies(&(*itEl), 1, i+1, 1, adj);
      if(mbErr != moab::MB_SUCCESS){
        throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not " + std::to_string(i+1) + " adjacencies"));
      }
      int index = mbIFace->id_from_handle(*itEl) - 1;
      int indexCell = 0;
      for(moab::Range::iterator itadj = adj.begin(); itadj != adj.end(); itadj++){
        (*elementAdjacencies)[i][index*nCellsPerElement + indexCell] = mbIFace->id_from_handle(*itadj) - 1;
        indexCell += 1;
      }
    }
  }
  moab::Range allEnts;
  mbErr = mbIFace->get_entities_by_handle(meshset, allEnts);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not get all entities"));
  }
  mbErr = mbIFace->delete_entities(allEnts);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not delete all entities"));
  }
  mbErr = mbIFace->delete_entities(&meshset, 1);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not clear meshset."));
  }
  delete mbIFace;
};


void generateCellNodes(const ReferenceElement * refElement, std::vector<double> * linCellNodes, 
    std::vector<double> * cellNodes){
  const std::vector< std::vector<double> > * refNodes = refElement->getNodes();
  int topDim = refElement->getDimension();
  int dim = linCellNodes->size()/(topDim+1);
  EMatrix locT(topDim, topDim);
  for(int i = 0; i < topDim; i++){
    for(int j = 0; j < topDim; j++){
      locT(j, i) = (*refNodes)[(i+1)][j] - (*refNodes)[0][j];
    }
  }
  EMatrix T(dim, topDim);
  for(int i = 0; i < topDim; i++){
    for(int j = 0; j < dim; j++){
      T(j, i) = (*linCellNodes)[(i+1)*dim + j] - (*linCellNodes)[j];
    }
  }
  T = T*locT.inverse();
  const EVector locV0 = Eigen::Map<const EVector>((*refNodes)[0].data(), topDim);
  const EVector v0 = Eigen::Map<const EVector>(linCellNodes->data(), dim);
  EVector v(topDim), node(dim);
  cellNodes->resize(dim*(refNodes->size()));
  for(int i = 0; i < refNodes->size(); i++){
    for(int j = 0; j < topDim; j++){
      v[j] = (*refNodes)[i][j];
    }
    node = T*(v-locV0) + v0;
    for(int j = 0; j < dim; j++){
      (*cellNodes)[i*dim + j] = node[j];
    }
  }
};
