#include "CGSolver.h"

namespace hfox{

void CGSolver::allocate(){
  if(!initialized){
    throw(ErrorHandle("CGSolver", "allocate", "must initialize the solver before allocating."));
  }
  if(myMesh == NULL){
    throw(ErrorHandle("CGSolver", "allocate", "must set the Mesh before allocating."));
  }
  if(linSystem == NULL){
    throw(ErrorHandle("CGSolver", "allocate", "must set the linear system before allocating."));
  }
  if(model == NULL){
    throw(ErrorHandle("CGSolver", "allocate", "must set the model before allocating."));
  }
  if(boundaryModel == NULL){
    throw(ErrorHandle("CGSolver", "allocate", "must set the boundary model before allocating."));
  }
  if(fieldMap->size() == 0){
    throw(ErrorHandle("CGSolver", "allocate", "must set the fields before allocating."));
  }
  std::map<std::string, Field * >::iterator it = fieldMap->find("Solution");
  if(it == fieldMap->end()){
    throw(ErrorHandle("CGSolver", "allocate", "the field map must have a Solution field."));
  }
  if(*(it->second->getFieldType()) != Node){
    throw(ErrorHandle("CGSolver", "allocate", "the Solution field must be a nodal field."));
  }
  nDOFsPerNode = (*(it->second->getNumObjPerEnt()))*(*(it->second->getNumValsPerObj()));
  calcSparsityPattern();
  int nNodes = myMesh->getNumberPoints();
  linSystem->allocate(nNodes*nDOFsPerNode, &diagSparsePattern, &offSparsePattern);
  model->allocate(nDOFsPerNode);
  boundaryModel->allocate(nDOFsPerNode);
  allocated = 1;
};//allocate

void CGSolver::assemble(){
  if(!(initialized and allocated)){
    throw(ErrorHandle("CGSolver", "assemble", "must initialize and allocate the solver before allocating."));
  }
  if(assembled){
    linSystem->clearSystem();
  }
  const ReferenceElement * refEl = myMesh->getReferenceElement();
  std::map<std::string, std::vector<double> > nodalFieldMap = prepareLocalFieldMap(Node);
  std::map<std::string, std::vector<double> > faceFieldMap = prepareLocalFieldMap(Face);

  std::vector<int> cell(refEl->getNumNodes());
  std::vector< std::vector<double> > nodes(refEl->getNumNodes(), std::vector<double>(myMesh->getNodeSpaceDimension(), 0.0));
  const AssemblyType * modAssemble = model->getAssemblyType();
  //element loop
  int iEl = 0;
  ProgressBar pb;
  if(verbose){
    pb.setIterIndex(&iEl);
    pb.setNumIterations(myMesh->getNumberCells());
    std::cout << "Assembly - Element loop (nElements = " << myMesh->getNumberCells() << ", nNodes/El = " << cell.size() << "):" << std::endl;
    pb.update();
  }
  for(iEl = 0; iEl < myMesh->getNumberCells(); iEl++){
    myMesh->getCell(iEl, &cell);
    myMesh->getSlicePoints(cell, &nodes);
    constructLocalFields(cell, &nodalFieldMap);
    model->setElementNodes(&nodes);
    model->setFieldMap(&nodalFieldMap);
    model->compute();
    switch(modAssemble->matrix){
      case Add:{
                 linSystem->addValsMatrix(cell, cell, model->getLocalMatrix()->transpose().data());
                 break;
               }
      case Set:{
                 linSystem->setValsMatrix(cell, cell, model->getLocalMatrix()->transpose().data());
                 break;
               }
    }
    switch(modAssemble->rhs){
      case Add:{
                 linSystem->addValsRHS(cell, model->getLocalRHS()->data());
                 break;
               }
      case Set:{
                 linSystem->setValsRHS(cell, model->getLocalRHS()->data());
                 break;
               }
    }
    if(verbose){
      pb.update();
    }
  }
  linSystem->assemble();

  cell.resize(refEl->getFaceElement()->getNumNodes());
  nodes.resize(0);
  nodes.resize(cell.size(), std::vector<double>(myMesh->getNodeSpaceDimension(), 0.0));
  modAssemble = boundaryModel->getAssemblyType();
  std::vector<int> slice(1, 0);
  const std::set<int> * boundaryFaces = myMesh->getBoundaryFaces();
  std::set<int>::const_iterator itFace;
  int index = 0;
  if(verbose){
    pb.setIterIndex(&index);
    pb.setNumIterations(boundaryFaces->size());
  }
  if(modAssemble->matrix == Set){
    if(verbose){
      std::cout << "Assembly - Zero loop (nFaces = " << boundaryFaces->size() << ", nNodes/Face = " << cell.size() << "):" << std::endl;
      pb.update();
    }
    for(itFace = boundaryFaces->begin(); itFace != boundaryFaces->end(); itFace++){
      myMesh->getFace(*itFace, &cell);
      linSystem->zeroOutRows(cell);
      if(verbose){
        index += 1;
        pb.update();
      }
    }
    linSystem->assemble();
  }
  if(verbose){
    index = 0;
    std::cout << "Assembly - Boundary loop (nFaces = " << boundaryFaces->size() << ", nNodes/Face = " << cell.size() << "):" << std::endl;
    pb.update();
  }
  for(itFace = boundaryFaces->begin(); itFace != boundaryFaces->end(); itFace++){
    myMesh->getFace(*itFace, &cell);
    myMesh->getSlicePoints(cell, &nodes);
    slice[0] = *itFace;
    constructLocalFields(slice, &faceFieldMap);
    boundaryModel->setElementNodes(&nodes);
    boundaryModel->setFieldMap(&faceFieldMap);
    boundaryModel->compute();
    switch(modAssemble->matrix){
      case Add:{
                 linSystem->addValsMatrix(cell, cell, boundaryModel->getLocalMatrix()->transpose().data());
                 break;
               }
      case Set:{
                 linSystem->setValsMatrix(cell, cell, boundaryModel->getLocalMatrix()->transpose().data());
                 break;
               }
    }
    switch(modAssemble->rhs){
      case Add:{
                 linSystem->addValsRHS(cell, boundaryModel->getLocalRHS()->data());
                 break;
               }
      case Set:{
                 linSystem->setValsRHS(cell, boundaryModel->getLocalRHS()->data());
                 break;
               }
    }
    if(verbose){
      index += 1;
      pb.update();
    }
  }
  linSystem->assemble();
  assembled = 1;
};//assemble

void CGSolver::solve(){
  if(!assembled){
    throw(ErrorHandle("CGSolver", "solve", "system must be assembled before solving"));
  }
  linSystem->solve((*fieldMap)["Solution"]->getValues());
};//solve

void CGSolver::calcSparsityPattern(){
  diagSparsePattern.resize(myMesh->getNumberPoints()*nDOFsPerNode, 0);
  offSparsePattern.resize(myMesh->getNumberPoints()*nDOFsPerNode, 0);
  std::vector< std::set<int> > sparsity(myMesh->getNumberPoints());
  int nNodesEl = myMesh->getReferenceElement()->getNumNodes();
  std::vector<int> cell(nNodesEl, 0);
  for(int iEl = 0; iEl < myMesh->getNumberCells(); iEl++){
    myMesh->getCell(iEl, &cell);
    for(int iCell = 0; iCell < nNodesEl; iCell++){
      std::set<int> * nodeSet = &(sparsity[cell[iCell]]);
      for(int i = 0; i < nNodesEl; i++){
        nodeSet->insert(cell[i]);
      }
    }
  }
  for(int iRow = 0; iRow < sparsity.size(); iRow++){
    for(int kDof = 0; kDof < nDOFsPerNode; kDof++){
      diagSparsePattern[iRow*nDOFsPerNode + kDof] = sparsity[iRow].size()*nDOFsPerNode;
    }
    //in parallel we have to divide the set into diagonal and off diagonal parts
  }
};//calcSparsityPattern

}//hfox
