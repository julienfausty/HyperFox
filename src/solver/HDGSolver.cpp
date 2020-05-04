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
  it = fieldMap->find("Tau");
  if(it == fieldMap->end()){
    throw(ErrorHandle("HDGSolver", "allocate", "the field map must have a Tau field."));
  }
  if(*(it->second->getFieldType()) != Face){
    throw(ErrorHandle("HDGSolver", "allocate", "the Tau field must be a face field."));
  }
  if(*(it->second->getNumObjPerEnt()) != myMesh->getReferenceElement()->getFaceElement()->getNumNodes()){
    throw(ErrorHandle("HDGSolver", "allocate", "the Tau field must have an object per element node."));
  }
  if(*(it->second->getNumValsPerObj()) != nDOFsPerNode){
    throw(ErrorHandle("HDGSolver", "allocate", "the Tau field must have the same number of values per object as the Solution field."));
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
  Q = new Field(myMesh, Cell, 1, dim*nDOFsPerNode*nNodesPEl*nDOFsPerNode*nFacesPEl*nNodesPFace);
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
  if(!(initialized and allocated)){
    throw(ErrorHandle("HDGSolver", "assemble", "the solver must be initialized and allocated before assembling."));
  }
  if(assembled){
    linSystem->clearSystem();
  }
  const ReferenceElement * refEl = myMesh->getReferenceElement();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int dim = myMesh->getNodeSpaceDimension();
  int nFacesPEl = refEl->getNumFaces();
  int nNodesPEl = refEl->getNumNodes();
  int nNodesPFc = fEl->getNumNodes();
  int startU = 0;
  int lenU = nNodesPEl*nDOFsPerNode;
  int startQ = lenU;
  int lenQ = lenU*dim;
  int startL = lenU + lenQ;
  int lenL = nFacesPEl*nNodesPFc*nDOFsPerNode;
  const std::vector< std::vector<int> > * nodeMap = refEl->getFaceNodes();
  const AssemblyType * modAssembly = model->getAssemblyType();
  std::map<std::string, std::vector<double> > locFieldMap;
  std::map<std::string, std::vector<double> > faceFieldMap = prepareLocalFieldMap(Face);
  std::map<std::string, std::vector<double> > cellFieldMap = prepareLocalFieldMap(Cell);
  std::map<std::string, std::vector<double> > nodalFieldMap = prepareLocalFieldMap(Node);
  std::map<std::string, std::vector<double> >::iterator itfm;
  std::vector< std::vector<double> > nodes(nNodesPEl, std::vector<double>(dim, 0.0));
  std::vector<int> cell(nNodesPEl, 0);
  std::vector<int> face(nNodesPFc, 0);
  std::vector<int> facesInCell(nFacesPEl, 0);
  std::vector<int> face2Cells(2, 0);
  std::vector<int> unitCell(1, 0);
  std::vector<int> matRowCols(lenL, 0);
  int sgn = 1;
  int offset = 0;
  EMatrix locS(lenL, lenL);
  EVector locS0(lenL);
  EMatrix invSqqSqu(lenQ, lenU);
  EMatrix invSqqSql(lenQ, lenL);
  Eigen::HouseholderQR<EMatrix> invSuuMinSuqinvSqqSqu(lenU, lenU);
  Eigen::HouseholderQR<EMatrix> invSqq(lenQ, lenQ);
  const EMatrix * locMat;
  const EVector * locVec;
  //element loop
  int iEl = 0;
  ProgressBar pb;
  if(verbose){
    pb.setIterIndex(&iEl);
    pb.setNumIterations(myMesh->getNumberCells());
    std::cout << "Assembly - Element loop (nElements = " << myMesh->getNumberCells() << ", nNodes/El = " << nNodesPEl << "):" << std::endl;
    pb.update();
  }
  for(iEl = 0; iEl < myMesh->getNumberCells(); iEl++){
    myMesh->getCell(iEl, &cell);
    myMesh->getSlicePoints(cell, &nodes);
    model->setElementNodes(&nodes);
    myMesh->getCell2Face(iEl, &facesInCell);
    unitCell[0] = iEl;
    constructLocalFields(facesInCell, &faceFieldMap);
    constructLocalFields(unitCell, &cellFieldMap);
    constructLocalFields(cell, &nodalFieldMap);
    for(itfm = nodalFieldMap.begin(); itfm != nodalFieldMap.end(); itfm++){locFieldMap[itfm->first] = itfm->second;}
    for(itfm = cellFieldMap.begin(); itfm != cellFieldMap.end(); itfm++){locFieldMap[itfm->first] = itfm->second;}
    for(int i = 0; i < nFacesPEl; i++){
      myMesh->getFace(facesInCell[i], &face);
      for(int j = 0; j < nNodesPFc; j++){
        matRowCols[i*nNodesPFc + j] = facesInCell[i] * nNodesPFc + std::distance(face.begin(), std::find(face.begin(), face.end(), cell[nodeMap->at(i)[j]]));
      }
      myMesh->getFace2Cell(facesInCell[i], &face2Cells);
    }
    //re order the face fields correctly for local ordering
    for(itfm = faceFieldMap.begin(); itfm != faceFieldMap.end(); itfm++){
      locFieldMap[itfm->first].resize((itfm->second).size());
      for(int i = 0; i < nFacesPEl; i++){
        offset = i*nNodesPFc;
        for(int j = 0; j <nNodesPFc; j++){
          locFieldMap[itfm->first][offset + j] = itfm->second[matRowCols[offset + j] - facesInCell[i] * nNodesPFc];
        }
      }
    }
    //orienting faces here using the sign of tau
    for(int i = 0; i < nFacesPEl; i++){
      myMesh->getFace2Cell(facesInCell[i], &face2Cells);
      sgn = 1;
      if(iEl == face2Cells[1]){
        sgn = -1;
      }
      offset = i*nNodesPFc;
      for(int j = 0; j < nNodesPFc; j++){
        locFieldMap["Tau"][offset + j] *= sgn;
      }
    }
    model->setFieldMap(&locFieldMap);
    model->compute();
    locMat = model->getLocalMatrix();
    locVec = model->getLocalRHS();
    invSqq.compute(locMat->block(startQ, startQ, lenQ, lenQ));
    invSqqSqu = invSqq.solve(locMat->block(startQ, startU, lenQ, lenU));
    invSqqSql = invSqq.solve(locMat->block(startQ, startL, lenQ, lenL));
    //(Suu - Suq(Sqq^-1)Squ)^-1
    invSuuMinSuqinvSqqSqu.compute(locMat->block(startU, startU, lenU, lenU) - locMat->block(startU, startQ, lenU, lenQ)*invSqqSqu);
    EMap<EMatrix> locU(U->getValues()->data() + iEl*lenU*lenL, lenU, lenL);
    EMap<EMatrix> locQ(Q->getValues()->data() + iEl*lenQ*lenL, lenQ, lenL);
    EMap<EVector> locU0(U0->getValues()->data() + iEl*lenU, lenU);
    EMap<EVector> locQ0(Q0->getValues()->data() + iEl*lenQ, lenQ);
    locU = -invSuuMinSuqinvSqqSqu.solve(locMat->block(startU, startL, lenU, lenL) - locMat->block(startU, startQ, lenU, lenQ)*invSqqSql);
    locU0 = invSuuMinSuqinvSqqSqu.solve(locVec->segment(startU, lenU));
    locQ = -invSqqSqu*locU - invSqqSql;
    locQ0 = -invSqqSqu*locU0;
    locS = locMat->block(startL, startU, lenL, lenU)*locU + locMat->block(startL, startQ, lenL, lenQ)*locQ + locMat->block(startL, startL, lenL, lenL);
    locS0 = locVec->segment(startL, lenL) - locMat->block(startL, startU, lenL, lenU)*locU0 - locMat->block(startL, startQ, lenL, lenQ)*locQ0;
    switch(modAssembly->matrix){
      case Add:{
                 linSystem->addValsMatrix(matRowCols, matRowCols, locS.transpose().data());
                 break;
               }
      case Set:{
                 linSystem->setValsMatrix(matRowCols, matRowCols, locS.transpose().data());
                 break;
               }
    }
    switch(modAssembly->rhs){
      case Add:{
                 linSystem->addValsRHS(matRowCols, locS0.data());
                 break;
               }
      case Set:{
                 linSystem->setValsRHS(matRowCols, locS0.data());
                 break;
               }
    }
    if(verbose){
      pb.update();
    }
  }
  linSystem->assemble();

  matRowCols.resize(nNodesPFc*nDOFsPerNode);
  cell.resize(nNodesPFc);
  nodes.resize(nNodesPFc);
  modAssembly = boundaryModel->getAssemblyType();
  const std::set<int> * boundaryFaces = myMesh->getBoundaryFaces();
  std::set<int>::const_iterator itFace;
  int index = 0;
  if(verbose){
    pb.setIterIndex(&index);
    pb.setNumIterations(boundaryFaces->size());
  }
  if(modAssembly->matrix == Set){
    if(verbose){
      std::cout << "Assembly - Zero loop (nFaces = " << boundaryFaces->size() << ", nNodes/Face = " << nNodesPFc << "):" << std::endl;
      pb.update();
    }
    for(itFace = boundaryFaces->begin(); itFace != boundaryFaces->end(); itFace++){
      std::iota(matRowCols.begin(), matRowCols.end(), (*itFace)*(matRowCols.size()));
      linSystem->zeroOutRows(matRowCols);
      if(verbose){
        index += 1;
        pb.update();
      }
    }
    linSystem->assemble();
  }
  if(verbose){
    index = 0;
    std::cout << "Assembly - Boundary loop (nFaces = " << boundaryFaces->size() << ", nNodes/Face = " << nNodesPFc << "):" << std::endl;
    pb.update();
  }
  for(itFace = boundaryFaces->begin(); itFace != boundaryFaces->end(); itFace++){
    myMesh->getFace(*itFace, &cell);
    myMesh->getSlicePoints(cell, &nodes);
    unitCell[0] = *itFace;
    constructLocalFields(unitCell, &faceFieldMap);
    boundaryModel->setElementNodes(&nodes);
    boundaryModel->setFieldMap(&faceFieldMap);
    boundaryModel->compute();
    std::iota(matRowCols.begin(), matRowCols.end(), (*itFace)*(matRowCols.size()));
    switch(modAssembly->matrix){
      case Add:{
                 linSystem->addValsMatrix(matRowCols, matRowCols, boundaryModel->getLocalMatrix()->transpose().data());
                 break;
               }
      case Set:{
                 linSystem->setValsMatrix(matRowCols, matRowCols, boundaryModel->getLocalMatrix()->transpose().data());
                 break;
               }
    }
    switch(modAssembly->rhs){
      case Add:{
                 linSystem->addValsRHS(matRowCols, boundaryModel->getLocalRHS()->data());
                 break;
               }
      case Set:{
                 linSystem->setValsRHS(matRowCols, boundaryModel->getLocalRHS()->data());
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

void HDGSolver::solve(){
  if(!assembled){
    throw(ErrorHandle("HDGSolver", "solve", "system must be assembled before solving"));
  }
  linSystem->solve((*fieldMap)["Trace"]->getValues());
  const ReferenceElement * refEl = myMesh->getReferenceElement();
  const std::vector< std::vector<int> > * nodeMap = refEl->getFaceNodes();
  const ReferenceElement * fEl = refEl->getFaceElement();
  int dim = myMesh->getNodeSpaceDimension();
  int nFacesPEl = refEl->getNumFaces();
  int nNodesPEl = refEl->getNumNodes();
  int nNodesPFc = fEl->getNumNodes();
  int lenU = nNodesPEl*nDOFsPerNode;
  int lenQ = lenU*dim;
  int lenT = nNodesPFc * nDOFsPerNode;
  int lenL = lenT * nFacesPEl;
  std::vector<int> facesInCell(nFacesPEl, 0);
  std::vector<int> cell(nNodesPEl, 0);
  std::vector<int> face(nNodesPFc, 0);
  EVector locL(lenL);
  int iEl = 0;
  ProgressBar pb;
  if(verbose){
    pb.setIterIndex(&iEl);
    pb.setNumIterations(myMesh->getNumberCells());
    std::cout << "Computing Solution and Flux - Element loop (nElements = " << myMesh->getNumberCells() << ", nNodes/El = " << nNodesPEl << "):" << std::endl;
    pb.update();
  }
  for(iEl = 0; iEl < myMesh->getNumberCells(); iEl++){
    myMesh->getCell(iEl, &cell);
    myMesh->getCell2Face(iEl, &facesInCell);
    EMap<EMatrix> locU(U->getValues()->data() + iEl*lenU*lenL, lenU, lenL);
    EMap<EMatrix> locQ(Q->getValues()->data() + iEl*lenQ*lenL, lenQ, lenL);
    EMap<EVector> locU0(U0->getValues()->data() + iEl*lenU, lenU);
    EMap<EVector> locQ0(Q0->getValues()->data() + iEl*lenQ, lenQ);
    for(int i = 0; i < nFacesPEl; i++){
      EMap<EVector> buff(fieldMap->at("Trace")->getValues()->data()+facesInCell[i]*lenT, lenT);
      myMesh->getFace(facesInCell[i], &face);
      for(int j = 0; j < nNodesPFc; j++){
        locL[i*lenT + j] = buff[std::distance(face.begin(), std::find(face.begin(), face.end(), cell[nodeMap->at(i)[j]]))];
      }
    }
    EMap<EVector> locSol(fieldMap->at("Solution")->getValues()->data()+iEl*lenU, lenU);
    EMap<EVector> locFlux(fieldMap->at("Flux")->getValues()->data()+iEl*lenQ, lenQ);
    locSol = locU*locL + locU0;
    locFlux = locQ*locL + locQ0;
    if(verbose){
      pb.update();
    }
  }
};//solve

}//hfox
