#include <string>
#include <vector>
#include <cstdio>
#include <exception>
#include <tuple>
#include <moab/Core.hpp>
#include "ErrorHandle.h"
#include "Mesh.h"
#include "ReferenceElement.h"
#include "HDF5Io.h"

/*!
 * \brief a conversion utility from linear simplicial meshes generated with gmsh to the native h5 format of higher user defined order
 */

using namespace hfox;

std::string usage();
void generateHigherOrderMesh(Mesh * hoMesh, std::string input);
void readMesh(std::string input, std::vector<double> * linNodes, std::vector< std::vector<int> > * linCells,
    std::vector< std::vector<int> > * elementAdjacencies);
void generateCellNodes(const ReferenceElement * refElement, std::vector<double> * linCellNodes, 
    std::vector<double> * cellNodes);

int main(int argc, char **argv){
  std::string inputFile, outputFile;
  int dimension = 0;
  int order = 0;
  bool verbose;
  std::string optstr("i:d:r:o:vh");
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
      case 'v':
        {
          verbose = 1;
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
    generateHigherOrderMesh(&hoMesh, inputFile);
  } catch(const std::exception & e){
    std::cout << "The mesh generation process failed with an exception:" << std::endl;
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  if(verbose){
    std::cout << "New mesh characteristics" << std::endl;
    std::cout << "Dimension: " << hoMesh.getDimension() << std::endl;
    std::cout << "Number of nodes: " << hoMesh.getNumberPoints() << std::endl;
    std::cout << "Number of cells: " << hoMesh.getNumberCells() << std::endl;
    std::cout << "Number of faces: " << hoMesh.getNumberFaces() << std::endl;
    std::cout << "Number of boundary faces: " << hoMesh.getBoundaryFaces()->size() << std::endl;
  }
  try{
    HDF5Io hdfio(&hoMesh);
    hdfio.write(outputFile);
  }catch(const std::exception & e){
    std::cout << "The mesh writing process failed with an exception:" << std::endl;
    std::cout << e.what() << std::endl;
    return EXIT_FAILURE;
  }
  return EXIT_SUCCESS;
};//main

std::string usage(){
  return std::string("Usage is: convertGmsh2H5HO [-i inputFile] [-d dimension] [-r requiredOrder] [-o outputFile] [-h]");
};

void generateHigherOrderMesh(Mesh * hoMesh, std::string input){
  int dim = hoMesh->getDimension();
  std::vector<double> linNodes;
  std::vector< std::vector<int> > linCells(dim);
  std::vector< std::vector<int> > elementAdjacencies(dim);
  try{
    readMesh(input, &linNodes, &linCells, &elementAdjacencies);
  } catch(const std::exception & e){
    throw(e);
  }

  int nCells = linCells[dim-1].size()/(dim+1);
  int nNodesPerEl = hoMesh->getReferenceElement()->getNumNodes();
  std::vector<double> hoNodes;
  std::vector< std::vector<int> > hoCells(dim+1);
  std::vector<int> hoElements;

  std::vector<int> nCellsPerTopDim(dim);
  std::vector<ReferenceElement *> refElements(dim);
  std::vector< std::vector<bool> > tracking(dim);
  for(int i = 0; i < dim; i++){
    nCellsPerTopDim[i] = elementAdjacencies[i].size()/nCells;
    refElements[i] = new ReferenceElement(i+1, hoMesh->getReferenceElement()->getOrder(), "simplex");
    tracking[i].resize(linCells[i].size()/nCellsPerTopDim[i], 0);
  }

  std::vector<double> cellNodes, elNodes, linCellNodes, linElNodes;
  std::vector<int> cell;
  int cellId;
  for(int iEl = 0; iEl < nCells; iEl++){
    linElNodes.resize((dim+1)*dim);
    for(int i = 0; i < dim+1; i++){
      for(int j = 0; j < dim; j++){
        linElNodes[i*dim + j] = linNodes[linCells[dim-1][iEl*(dim+1) + i]*dim + j];
      }
    }
    generateCellNodes(hoMesh->getReferenceElement(), &linElNodes, &elNodes);
    for(int i = 0; i < dim; i++){
      linCellNodes.resize((i+2)*dim);
      for(int j = 0; j < nCellsPerTopDim[i]; j++){
        cellId = elementAdjacencies[i][iEl*nCellsPerTopDim[i] + j];
        if(!tracking[i][cellId]){
          for(int l = 0; l < (i+2); l++){
            for(int k = 0; k < dim ; k++){
              linCellNodes[l*dim + k] = linNodes[linCells[i][cellId*(i+2) + l]*dim + k];
            }
          }
          generateCellNodes(refElements[i], &linCellNodes, &cellNodes);
          //put these new nodes in the hoNodes
          //generate the new cell
          //put that in hoCells
          //set tracking to true
        }else{
          int nNodesPerCell = refElements[i]->getNumNodes();
          cell.resize(nNodesPerCell);
          cell.assign(hoCells[i].begin() + cellId*nNodesPerCell, hoCells[i].begin() + (cellId+1)*nNodesPerCell);
          //get the cellNodes from hoNodes
        }
        //put the cell in hoElements in the right order using elNodes (only inner nodes though)
      }
    }
  }
  hoMesh->setMesh(dim, hoNodes, hoElements);
  for(int i = 0; i < refElements.size(); i++){
    delete refElements[i];
  }
};

void readMesh(std::string input, std::vector<double> * linNodes, std::vector< std::vector<int> > * linCells,
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
  linNodes->resize(vertexes.size()*dim);
  std::vector<double> vertex(3);
  for(moab::Range::iterator it = vertexes.begin(); it != vertexes.end(); it++){
    mbErr = mbIFace->get_coords(&(*it), 1, vertex.data());
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not get vertex coords from list"));
    }
    int index = mbIFace->id_from_handle(*it) -1;
    for(int i = 0; i < dim ; i++){
      (*linNodes)[index*dim + i] = vertex[i];
    }
  }
  //fill cells
  moab::Range elems;
  mbErr = mbIFace->get_entities_by_dimension(meshset, dim, elems);
  if(mbErr != moab::MB_SUCCESS){
    throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not get elements from meshset"));
  }
  for(int i = 0; i < dim; i++){
    moab::Range tmp;
    mbErr = mbIFace->get_adjacencies(elems, i+1, 1, tmp, moab::Interface::UNION);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not develop face adjacencies from cells (moab)"));
    }
    mbErr = mbIFace->add_entities(meshset, tmp);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("Mesh", "computeFaces", "could not add faces to meshset (moab)"));
    }
    mbErr = mbIFace->get_entities_by_dimension(meshset, i+1, tmp);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("convertGmsh2H5HO", "readMesh", 
            "could not get " + std::to_string(i) + " cells from meshset"));
    }
    std::vector<moab::EntityHandle> verts;
    int nNodesPerCell = i+2;
    (*linCells)[i].resize(nNodesPerCell*tmp.size());
    for(moab::Range::iterator it = tmp.begin(); it != tmp.end(); it++){
      mbErr = mbIFace->get_connectivity(&(*it), 1, verts);
      if(mbErr != moab::MB_SUCCESS){
        throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not get vertices from element"));
      }
      int index = mbIFace->id_from_handle(*it) - 1;
      for(int j = 0; j < nNodesPerCell; j++){
        (*linCells)[i][index*nNodesPerCell + j] = mbIFace->id_from_handle(verts[j])-1;
      }
    }
    moab::Range adj;
    mbErr = mbIFace->get_adjacencies(&(*(elems.begin())), 1, i+1, 1, adj);
    if(mbErr != moab::MB_SUCCESS){
      throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not " + std::to_string(i) + " adjacencies"));
    }
    int nCellsPerElement = adj.size();
    (*elementAdjacencies)[i].resize(elems.size()*nCellsPerElement);
    for(moab::Range::iterator itEl = elems.begin(); itEl != elems.end(); itEl++){
      mbErr = mbIFace->get_adjacencies(&(*itEl), 1, i+1, 1, adj);
      if(mbErr != moab::MB_SUCCESS){
        throw(ErrorHandle("convertGmsh2H5HO", "readMesh", "could not " + std::to_string(i) + " adjacencies"));
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
};
