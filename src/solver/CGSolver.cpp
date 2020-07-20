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
  Partitioner * part = myMesh->getPartitioner();
  const ReferenceElement * refEl = myMesh->getReferenceElement();
  std::map<std::string, std::vector<double> > nodalFieldMap = prepareLocalFieldMap(Node);
  std::map<std::string, std::vector<double> > faceFieldMap = prepareLocalFieldMap(Face);

  std::vector<int> cell(refEl->getNumNodes());
  std::vector<int> locCell(refEl->getNumNodes());
  std::vector< std::vector<double> > nodes(refEl->getNumNodes(), std::vector<double>(myMesh->getNodeSpaceDimension(), 0.0));
  const AssemblyType * modAssemble = model->getAssemblyType();
  //element loop
  int iEl = 0;
  ProgressBar pb;
  if(verbose){
    pb.setIterIndex(&iEl);
    pb.setNumIterations(myMesh->getNumberCells());
    int totCells = myMesh->getNumberCells();
    bool speak = 1;
    if(part != NULL){
      totCells = part->getTotalNumberEls();
      if(part->getRank() != 0){
        speak = 0;
      }
    }
    if(speak){
      std::cout << "Assembly - Element loop (nElements = " << totCells << ", nNodes/El = " << cell.size() << "):" << std::endl;
    }
    pb.update();
  }
  EMatrix modelT(model->getLocalMatrix()->rows(), model->getLocalMatrix()->cols());
  for(iEl = 0; iEl < myMesh->getNumberCells(); iEl++){
    std::cout << "iEl/totEl: " << iEl << "/" << myMesh->getNumberCells() << std::endl;
    myMesh->getCell(iEl, &cell);
    if(part == NULL){
      constructLocalFields(cell, &nodalFieldMap);
      myMesh->getSlicePoints(cell, &nodes);
    } else {
      part->global2LocalNodeSlice(cell, &locCell);
      constructLocalFields(cell, locCell, &nodalFieldMap);
      for(int i = 0; i < locCell.size(); i++){
        if(locCell[i] != -1){
          myMesh->getPoint(locCell[i], &(nodes[i]));
        } else {
          myMesh->getGhostPoint(cell[i], &(nodes[i]));
        }
      }
    }
    model->setElementNodes(&nodes);
    model->setFieldMap(&nodalFieldMap);
    model->compute();
    modelT = model->getLocalMatrix()->transpose();
    switch(modAssemble->matrix){
      case Add:{
                 linSystem->addValsMatrix(cell, cell, modelT.data());
                 break;
               }
      case Set:{
                 linSystem->setValsMatrix(cell, cell, modelT.data());
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
  locCell.resize(refEl->getFaceElement()->getNumNodes());
  nodes.resize(0);
  nodes.resize(cell.size(), std::vector<double>(myMesh->getNodeSpaceDimension(), 0.0));
  modelT.resize(boundaryModel->getLocalMatrix()->rows(), boundaryModel->getLocalMatrix()->cols());
  modAssemble = boundaryModel->getAssemblyType();
  std::vector<int> slice(1, 0);
  std::vector<int> locSlice(1, 0);
  const std::set<int> * boundaryFaces = myMesh->getBoundaryFaces();
  std::set<int>::const_iterator itFace;
  int index = 0;
  if(verbose){
    pb.setIterIndex(&index);
    pb.setNumIterations(boundaryFaces->size());
  }
  if(modAssemble->matrix == Set){
    if(verbose){
      int totCells = boundaryFaces->size();
      bool speak = 1;
      if(part != NULL){
        int locBs = boundaryFaces->size();
        MPI_Reduce(&locBs, &totCells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
        if(part->getRank() != 0){
          speak = 0;
        }
      }
      if(speak){
        std::cout << "Assembly - Zero loop (nFaces = " << totCells << ", nNodes/Face = " << cell.size() << "):" << std::endl;
      }
      pb.update();
    }
    std::set<int> zeroSet;
    for(itFace = boundaryFaces->begin(); itFace != boundaryFaces->end(); itFace++){
      if(part == NULL){
        myMesh->getFace(*itFace, &cell);
      } else {
        myMesh->getFace(part->global2LocalFace(*itFace), &cell);
      }
      for(int i = 0; i < cell.size(); i++){ 
        zeroSet.insert(cell[i]);
      }
      if(verbose){
        index += 1;
        pb.update();
      }
    }
    std::vector<int> zeroRows(zeroSet.size(), 0);
    std::copy(zeroSet.begin(), zeroSet.end(), zeroRows.begin());
    linSystem->zeroOutRows(zeroRows);
    linSystem->assemble();
  }
  if(verbose){
    index = 0;
    int totCells = boundaryFaces->size();
    bool speak = 1;
    if(part != NULL){
      int locBs = boundaryFaces->size();
      MPI_Reduce(&locBs, &totCells, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
      if(part->getRank() != 0){
        speak = 0;
      }
    }
    if(speak){
      std::cout << "Assembly - Boundary loop (nFaces = " << totCells << ", nNodes/Face = " << cell.size() << "):" << std::endl;
    }
    pb.update();
  }
  std::cout << "starting boundary assembly" << std::endl;
  for(itFace = boundaryFaces->begin(); itFace != boundaryFaces->end(); itFace++){
    if(part == NULL){
      myMesh->getFace(*itFace, &cell);
      myMesh->getSlicePoints(cell, &nodes);
      slice[0] = *itFace;
      constructLocalFields(slice, &faceFieldMap);
    } else {
      locSlice[0] = part->global2LocalFace(*itFace);
      myMesh->getFace(locSlice[0], &cell);
      part->global2LocalNodeSlice(cell, &locCell);
      slice[0] = *itFace;
      constructLocalFields(slice, locSlice, &faceFieldMap);
      for(int i = 0; i < locCell.size(); i++){
        if(locCell[i] != -1){
          myMesh->getPoint(locCell[i], &(nodes[i]));
        } else {
          myMesh->getGhostPoint(cell[i], &(nodes[i]));
        }
      }
    }
    boundaryModel->setElementNodes(&nodes);
    boundaryModel->setFieldMap(&faceFieldMap);
    boundaryModel->compute();
    modelT = boundaryModel->getLocalMatrix()->transpose();
    switch(modAssemble->matrix){
      case Add:{
                 linSystem->addValsMatrix(cell, cell, modelT.data());
                 break;
               }
      case Set:{
                 linSystem->setValsMatrix(cell, cell, modelT.data());
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
  std::vector<double> * solutionVec = (*fieldMap)["Solution"]->getValues();
  linSystem->solve(solutionVec);
  Partitioner * part = myMesh->getPartitioner();
  if(part != NULL){
    std::vector<int> thisRange;
    linSystem->getSolutionOwnership(&thisRange);
    //int nSendDofs = thisRange.size();
    //std::vector<int> nDofs(part->getNumPartitions());
    //std::vector<int> offsets(part->getNumPartitions());
    //MPI_Allgather(&nSendDofs, 1, MPI_INT, nDofs.data(), 1, MPI_INT, MPI_COMM_WORLD);
    //for(int i = 0; i < nDofs.size(); i++){
      //offsets[i] = std::accumulate(nDofs.begin(), nDofs.begin() + i, 0);
    //}
    //std::vector<int> allDofIds(std::accumulate(nDofs.begin(), nDofs.end(), 0));
    //std::vector<double> allDofVals(allDofIds.size());
    //MPI_Allgatherv(thisRange.data(), thisRange.size(), MPI_INT, allDofIds.data(), nDofs.data(), offsets.data(), MPI_INT, MPI_COMM_WORLD);
    //MPI_Allgatherv(solutionVec->data(), solutionVec->size(), MPI_DOUBLE, allDofVals.data(), nDofs.data(), offsets.data(), MPI_DOUBLE, MPI_COMM_WORLD);
    //std::vector<int>::iterator it;
    //for(int i = 0; i < part->getNodeIds()->size(); i++){
      //for(int k = 0; k < nDOFsPerNode; k++){
        //it = std::find(allDofIds.begin(), allDofIds.end(), part->getNodeIds()->at(i)*nDOFsPerNode + k);
        //if(it != allDofIds.end()){
          //(*solutionVec)[i*nDOFsPerNode + k] = allDofVals[std::distance(allDofIds.begin(), it)];
        //} else {
          //throw(ErrorHandle("CGSolver", "solve", "one of the dof indexes of the partition could not be found in the communication structure."));
        //}
      //}
    //}
    std::vector<int> dofsWeNeed(myMesh->getNumberPoints() * nDOFsPerNode, 0);
    Utils::multiplyIndexes(nDOFsPerNode, part->getNodeIds(), &dofsWeNeed);
    std::vector<int> currentMap, sendIds;
    std::vector<double> sendVals;
    std::vector<int>::iterator it;
    for(int i = 0; i < thisRange.size(); i++){
      it = std::find(dofsWeNeed.begin(), dofsWeNeed.end(), thisRange[i]);
      if(it != dofsWeNeed.end()){
        currentMap.push_back(i);
        currentMap.push_back(std::distance(dofsWeNeed.begin(), it));
      } else {
        sendIds.push_back(thisRange[i]);
        sendVals.push_back(solutionVec->at(i));
      }
    }
    //do the reordering we can currently do on this partition
    std::vector<double> buffer(solutionVec->size(), 0.0);
    for(int i = 0; i < currentMap.size()/2; i++){
      buffer[currentMap[2*i+1]] = solutionVec->at(currentMap[2*i]);
    }
    std::copy(buffer.begin(), buffer.end(), solutionVec->begin());
    buffer.clear();
    std::vector<int> recvIds;
    for(int i = 0; i < dofsWeNeed.size(); i++){
      it = std::find(thisRange.begin(), thisRange.end(), dofsWeNeed[i]);
      if(it == thisRange.end()){
        recvIds.push_back(dofsWeNeed[i]);
        recvIds.push_back(i);
      }
    }
    //gather the mismatched dofs to all nodes
    int nSendDofs = sendIds.size();
    std::vector<int> nDofs(part->getNumPartitions());
    std::vector<int> offsets(part->getNumPartitions());
    MPI_Allgather(&nSendDofs, 1, MPI_INT, nDofs.data(), 1, MPI_INT, MPI_COMM_WORLD);
    for(int i = 0; i < nDofs.size(); i++){
      offsets[i] = std::accumulate(nDofs.begin(), nDofs.begin() + i, 0);
    }
    std::vector<int> allDofIds(std::accumulate(nDofs.begin(), nDofs.end(), 0));
    std::vector<double> allDofVals(allDofIds.size());
    MPI_Allgatherv(sendIds.data(), sendIds.size(), MPI_INT, allDofIds.data(), nDofs.data(), offsets.data(), MPI_INT, MPI_COMM_WORLD);
    MPI_Allgatherv(sendVals.data(), sendVals.size(), MPI_DOUBLE, allDofVals.data(), nDofs.data(), offsets.data(), MPI_DOUBLE, MPI_COMM_WORLD);
    for(int i = 0; i < recvIds.size()/2; i++){
      it = std::find(allDofIds.begin(), allDofIds.end(), recvIds[2*i]);
      if(it != allDofIds.end()){
        (*solutionVec)[recvIds[2*i + 1]] = allDofVals[std::distance(allDofIds.begin(), it)];
      } else {
        throw(ErrorHandle("CGSolver", "solve", "one of the dof indexes of the partition could not be found in the communication structure."));
      }
    }
  }
};//solve

void CGSolver::calcSparsityPattern(){
  Partitioner * part = myMesh->getPartitioner();
  diagSparsePattern.resize(myMesh->getNumberPoints()*nDOFsPerNode, 0);
  offSparsePattern.resize(myMesh->getNumberPoints()*nDOFsPerNode, 0);
  std::map<int, std::vector<int> > ghostSparse;
  std::vector< std::set<int> > sparsity(myMesh->getNumberPoints());
  int nNodesEl = myMesh->getReferenceElement()->getNumNodes();
  std::vector<int> cell(nNodesEl, 0);
  std::vector<int> locCell(nNodesEl, 0);
  for(int iEl = 0; iEl < myMesh->getNumberCells(); iEl++){
    myMesh->getCell(iEl, &cell);
    if(part != NULL){
      part->global2LocalNodeSlice(cell, &locCell);
    } else {
      locCell = cell;
    }
    for(int iCell = 0; iCell < nNodesEl; iCell++){
      if(locCell[iCell] != -1){
        std::set<int> * nodeSet = &(sparsity[locCell[iCell]]);
        for(int i = 0; i < nNodesEl; i++){
          nodeSet->insert(cell[i]);
        }
      } 
    }
  }
  //iterate through ghost cells to add to sparsity pattern
  int nVals = 2 + nNodesEl;
  const std::vector<int> * pGhostCells = myMesh->getGhostCells();
  for(int i = 0; i < pGhostCells->size()/nVals; i++){
    std::copy(pGhostCells->begin() + nVals*i + 2, pGhostCells->begin() + nVals*(i+1), cell.begin());
    part->global2LocalNodeSlice(cell, &locCell);
    for(int iCell = 0; iCell < nNodesEl; iCell++){
      if(locCell[iCell] != -1){
        std::set<int> * nodeSet = &(sparsity[locCell[iCell]]);
        for(int i = 0; i < nNodesEl; i++){
          nodeSet->insert(cell[i]);
        }
      } 
    }
  }
  std::vector<int> locBuff;
  for(int iRow = 0; iRow < sparsity.size(); iRow++){
    int rowSize = sparsity[iRow].size()*nDOFsPerNode;
    cell.resize(sparsity[iRow].size());
    std::copy(sparsity[iRow].begin(), sparsity[iRow].end(), cell.begin());
    locBuff.resize(cell.size());
    if(part != NULL){
      part->global2LocalNodeSlice(cell, &locBuff);
    } else {
      locBuff = cell;
    }
    int offDiag = 0;
    for(int k = 0; k < locBuff.size(); k++){
      if(locBuff[k] == -1){
        offDiag += nDOFsPerNode;
      }
    }
    rowSize -= offDiag;
    for(int kDof = 0; kDof < nDOFsPerNode; kDof++){
      diagSparsePattern[iRow*nDOFsPerNode + kDof] = rowSize;
      offSparsePattern[iRow*nDOFsPerNode + kDof] = offDiag;
    }
  }
};//calcSparsityPattern

}//hfox
