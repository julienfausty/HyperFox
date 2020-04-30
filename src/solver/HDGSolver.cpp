#include "HDGSolver.h"

namespace hfox{

void HDGSolver::allocate(){
  if(!initialized){
    throw(ErrorHandle("HDGSolver", "allocate", "must initialize the solver before allocating."));
  }
  if(myMesh == NULL){
    throw(ErrorHandle("HDGSolver", "allocate", "must set the Mesh before allocating."));
  }
  if(linSystem == NULL){
    throw(ErrorHandle("HDGSolver", "allocate", "must set the linear system before allocating."));
  }
  if(model == NULL){
    throw(ErrorHandle("HDGSolver", "allocate", "must set the model before allocating."));
  }
  if(boundaryModel == NULL){
    throw(ErrorHandle("HDGSolver", "allocate", "must set the boundary model before allocating."));
  }
  if(fieldMap->size() == 0){
    throw(ErrorHandle("HDGSolver", "allocate", "must set the fields before allocating."));
  }
  std::map<std::string, Field * >::iterator it = fieldMap->find("Solution");
  if(it == fieldMap->end()){
    throw(ErrorHandle("HDGSolver", "allocate", "the field map must have a Solution field."));
  }
  if(*(it->second->getFieldType()) != Cell){
    throw(ErrorHandle("HDGSolver", "allocate", "the Solution field must be a cell field."));
  }
  if(*(it->second->getNumObjPerEnt()) != myMesh->getReferenceElement()->getNumNodes()){
    throw(ErrorHandle("HDGSolver", "allocate", "the Solution field must have an object per element node."));
  }
  nDOFsPerNode = *(fieldMap->at("Solution")->getNumValsPerObj());
  it = fieldMap->find("Flux");
  if(it == fieldMap->end()){
    throw(ErrorHandle("HDGSolver", "allocate", "the field map must have a Flux field."));
  }
  if(*(it->second->getFieldType()) != Cell){
    throw(ErrorHandle("HDGSolver", "allocate", "the Flux field must be a cell field."));
  }
  if(*(it->second->getNumObjPerEnt()) != myMesh->getReferenceElement()->getNumNodes()){
    throw(ErrorHandle("HDGSolver", "allocate", "the Flux field must have an object per element node."));
  }
  if(*(it->second->getNumValsPerObj()) != nDOFsPerNode*(myMesh->getNodeSpaceDimension())){
    throw(ErrorHandle("HDGSolver", "allocate", "the Flux field must represent a spatial derivative of the Solution field."));
  }
  it = fieldMap->find("Trace");
  if(it == fieldMap->end()){
    throw(ErrorHandle("HDGSolver", "allocate", "the field map must have a Trace field."));
  }
  if(*(it->second->getFieldType()) != Face){
    throw(ErrorHandle("HDGSolver", "allocate", "the Trace field must be a face field."));
  }
  if(*(it->second->getNumObjPerEnt()) != myMesh->getReferenceElement()->getFaceElement()->getNumNodes()){
    throw(ErrorHandle("HDGSolver", "allocate", "the Trace field must have an object per element node."));
  }
  if(*(it->second->getNumValsPerObj()) != nDOFsPerNode){
    throw(ErrorHandle("HDGSolver", "allocate", "the Trace field must have the same number of values per object as the Solution field."));
  }
  calcSparsityPattern();
  int dim = myMesh->getNodeSpaceDimension();
  int nNodesPEl = myMesh->getReferenceElement()->getNumNodes();
  int nFacesPEl = myMesh->getReferenceElement()->getNumFaces();
  int nFaces = myMesh->getNumberFaces();
  int nNodesPFace = myMesh->getReferenceElement()->getFaceElement()->getNumNodes();
  linSystem->allocate(nDOFsPerNode*nNodesPFace*nFaces, &diagSparsePattern, &offSparsePattern);
  model->allocate(nDOFsPerNode);
  boundaryModel->allocate(nDOFsPerNode);
  if(U != NULL){delete U;}
  U = new Field(myMesh, Cell, 1, nDOFsPerNode*nNodesPEl*nDOFsPerNode*nFacesPEl*nNodesPFace);
  if(Q != NULL){delete Q;}
  Q = new Field(myMesh, Cell, 1, dim*nDOFsPerNode*nNodesPEl*dim*nDOFsPerNode*nFacesPEl*nNodesPFace);
  if(U0 != NULL){delete U0;}
  U0 = new Field(myMesh, Cell, 1, nDOFsPerNode*nNodesPEl);
  if(Q0 != NULL){delete Q0;}
  Q0 = new Field(myMesh, Cell, 1, dim*nDOFsPerNode*nNodesPEl);
  allocated = 1;
};//allocate

HDGSolver::~HDGSolver(){
  delete U;
  delete Q;
  delete U0;
  delete Q0;
}

void HDGSolver::calcSparsityPattern(){
  const ReferenceElement * refEl = myMesh->getReferenceElement();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int nNodesPFace = fEl->getNumNodes();
  int nFaces = myMesh->getNumberFaces();
  int nFacesPEl = refEl->getNumFaces();
  diagSparsePattern.resize(nFaces*nNodesPFace*nDOFsPerNode, 0);
  offSparsePattern.resize(nFaces*nNodesPFace*nDOFsPerNode, 0);
  std::vector< std::set<int> > sparsity(nFaces*nNodesPFace);
  std::set<int> * nodeSet;
  std::vector<int> facesInEl(nFacesPEl, 0);
  for(int iEl = 0; iEl < myMesh->getNumberCells(); iEl++){
    myMesh->getCell2Face(iEl, &facesInEl);
    for(int locFc = 0; locFc < nFacesPEl; locFc++){
      for(int iNode = 0; iNode < nNodesPFace; iNode++){
        nodeSet = &(sparsity[facesInEl[locFc]*nNodesPFace + iNode]);
        for(int i = 0; i < nFacesPEl; i++){
          for(int j = 0; j < nNodesPFace; j++){
            nodeSet->insert(facesInEl[i]*nNodesPFace + j);
          }
        }
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

void HDGSolver::assemble(){

};//assemble

void HDGSolver::solve(){

};//solve

}//hfox
