#include "HDGSolver.h"

namespace hfox{

void HDGSolver::allocate(){
  if(!initialized){
    throw(ErrorHandle("HDGSolver", "allocate", "must initialize the solver before allocating."));
  }
  if(myMesh == NULL){
    throw(ErrorHandle("HDGSolver", "allocate", "must set the Mesh before allocating."));
  }
  if((linSystem == NULL) and (myOpts.type != SEXPLICIT)){
    throw(ErrorHandle("HDGSolver", "allocate", "must set the linear system before allocating."));
  }
  if(model == NULL){
    throw(ErrorHandle("HDGSolver", "allocate", "must set the model before allocating."));
  }
  if(bSystem.getBoundaryList()->empty()){
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
  if(((*(it->second->getNumValsPerObj()) != nDOFsPerNode*nDOFsPerNode) and (*(it->second->getNumValsPerObj()) != 2*nDOFsPerNode*nDOFsPerNode))){
    throw(ErrorHandle("HDGSolver", "allocate", "the Tau field must have the same or twice the number of values per object as the Solution field."));
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
  int dim = myMesh->getNodeSpaceDimension();
  int nNodesPEl = myMesh->getReferenceElement()->getNumNodes();
  int nFacesPEl = myMesh->getReferenceElement()->getNumFaces();
  int nFaces = myMesh->getNumberFaces();
  int nNodesPFace = myMesh->getReferenceElement()->getFaceElement()->getNumNodes();
  if(myOpts.type != SEXPLICIT){
    calcSparsityPattern();
    linSystem->allocate(nDOFsPerNode*nNodesPFace*nFaces, &diagSparsePattern, &offSparsePattern);
  } else{
    int nParts;
    MPI_Comm_size(MPI_COMM_WORLD, &nParts);
    if(nParts != 1){
      throw(ErrorHandle("HDGSolver", "allocate", "the SEXPLICIT mode is currently unsuported in parallel, please use the WEXPLICIT mode instead."));
    }
    if(S != NULL){delete S;}
    S = new Field(myMesh, Face, 1, std::pow(nDOFsPerNode*nNodesPFace, 2));
    if(S0 != NULL){delete S0;}
    S0 = new Field(myMesh, Face, 1, nDOFsPerNode*nNodesPFace);
  }
  model->allocate(nDOFsPerNode);
  for(int iBoundary = 0; iBoundary < bSystem.getBoundaryList()->size(); iBoundary++){
    std::get<0>(bSystem.getBoundaryList()->at(iBoundary))->allocate(nDOFsPerNode);
  }
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
  if(U != NULL){delete U;}
  if(Q != NULL){delete Q;}
  if(S != NULL){delete S;}
  if(U0 != NULL){delete U0;}
  if(Q0 != NULL){delete Q0;}
  if(S0 != NULL){delete S0;}
};//destructor

void HDGSolver::calcSparsityPattern(){
  Partitioner * part = myMesh->getPartitioner();
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
  std::vector<int> locCellIndexes;
  for(int iEl = 0; iEl < myMesh->getNumberCells(); iEl++){
    myMesh->getCell2Face(iEl, &facesInEl);
    locCellIndexes = facesInEl;
    if(part != NULL){
      part->global2LocalFaceSlice(facesInEl, &locCellIndexes);
    }
    for(int locFc = 0; locFc < nFacesPEl; locFc++){
      if(locCellIndexes[locFc] != -1){
        for(int iNode = 0; iNode < nNodesPFace; iNode++){
          nodeSet = &(sparsity[locCellIndexes[locFc]*nNodesPFace + iNode]);
          for(int i = 0; i < nFacesPEl; i++){
            for(int j = 0; j < nNodesPFace; j++){
              nodeSet->insert(locCellIndexes[i]*nNodesPFace + j);
            }
          }
        }
      }
    }
  }
  std::set<int>::iterator it;
  int diagSize, offDiagSize;
  for(int iRow = 0; iRow < sparsity.size(); iRow++){
    offDiagSize = 0;
    for(it = sparsity[iRow].begin(); it != sparsity[iRow].end(); it++){
      if(*it < 0){
        offDiagSize += 1;
      }
    }
    diagSize = sparsity[iRow].size() - offDiagSize;
    for(int kDof = 0; kDof < nDOFsPerNode; kDof++){
      diagSparsePattern[iRow*nDOFsPerNode + kDof] = diagSize*nDOFsPerNode;
      offSparsePattern[iRow*nDOFsPerNode + kDof] = offDiagSize*nDOFsPerNode;
    }
  }
};//calcSparsityPattern

void HDGSolver::assemble(){
  if(!(initialized and allocated)){
    throw(ErrorHandle("HDGSolver", "assemble", "the solver must be initialized and allocated before assembling."));
  }
  if(myOpts.type != SEXPLICIT){
    if(assembled){
      linSystem->clearSystem();
    }
  } else {
    std::fill(S->getValues()->begin(), S->getValues()->end(), 0.0);
    std::fill(S0->getValues()->begin(), S0->getValues()->end(), 0.0);
  }
  Partitioner * part = myMesh->getPartitioner();
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
  int matLen = nDOFsPerNode*nDOFsPerNode;
  const std::vector< std::vector<int> > * nodeMap = refEl->getFaceNodes();
  const AssemblyType * modAssembly = model->getAssemblyType();
  std::map<std::string, std::vector<double> > locFieldMap;
  std::map<std::string, std::vector<double> > faceFieldMap = prepareLocalFieldMap(Face);
  std::map<std::string, std::vector<double> > cellFieldMap = prepareLocalFieldMap(Cell);
  std::map<std::string, std::vector<double> > nodalFieldMap = prepareLocalFieldMap(Node);
  std::map<std::string, std::vector<double> >::iterator itfm;
  std::vector< std::vector<double> > nodes(nNodesPEl, std::vector<double>(dim, 0.0));
  std::vector<int> cell(nNodesPEl, 0);
  std::vector<int> locCell(nNodesPEl, 0);
  std::vector<int> face(nNodesPFc, 0);
  std::vector<int> locFace(nNodesPFc, 0);
  std::vector<int> facesInCell(nFacesPEl, 0);
  std::vector<int> locFacesInCell(nFacesPEl, 0);
  std::vector<int> face2Cells(2, 0);
  std::vector<int> locFace2Cells(2, 0);
  std::vector<int> unitCell(1, 0);
  std::vector<int> locUnitCell(1, 0);
  std::vector<int> face2CellMap(nFacesPEl*nNodesPFc, 0);
  std::vector<int> matRowCols(lenL, 0);
  std::vector<double> fieldBuffer;
  std::vector<int>::iterator intIt;
  int buff = 0;
  int offset = 0;
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> locS(lenL, lenL);
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
    myMesh->getCell2Face(iEl, &facesInCell);
    unitCell[0] = iEl;
    if(part == NULL){
      myMesh->getSlicePoints(cell, &nodes);
      constructLocalFields(facesInCell, &faceFieldMap);
      constructLocalFields(unitCell, &cellFieldMap);
      constructLocalFields(cell, &nodalFieldMap);
    } else {
      part->global2LocalFaceSlice(facesInCell, &locFacesInCell);
      constructLocalFields(facesInCell, locFacesInCell, &faceFieldMap);
      unitCell[0] = part->local2GlobalEl(iEl);
      locUnitCell[0] = iEl;
      constructLocalFields(unitCell, locUnitCell, &cellFieldMap);
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
    for(itfm = nodalFieldMap.begin(); itfm != nodalFieldMap.end(); itfm++){locFieldMap[itfm->first] = itfm->second;}
    for(itfm = cellFieldMap.begin(); itfm != cellFieldMap.end(); itfm++){locFieldMap[itfm->first] = itfm->second;}
    for(int i = 0; i < nFacesPEl; i++){
      if(part == NULL){
        myMesh->getFace(facesInCell[i], &face);
      } else {
        if(locFacesInCell[i] != -1){
          myMesh->getFace(locFacesInCell[i], &face);
        } else {
          myMesh->getGhostFace(facesInCell[i], &face);
        }
      }
      for(int j = 0; j < nNodesPFc; j++){
        intIt = std::find(face.begin(), face.end(), cell[nodeMap->at(i)[j]]);
        if(intIt == face.end()){
          throw(ErrorHandle("HDGSolver", "assemble", "couldn't find cell node in face."));
        }
        face2CellMap[i*nNodesPFc + j] = std::distance(face.begin(), intIt);
        for(int k = 0; k < nDOFsPerNode; k++){
          matRowCols[(i*nNodesPFc + j)*nDOFsPerNode + k] = (facesInCell[i] * nNodesPFc + face2CellMap[i*nNodesPFc + j])*nDOFsPerNode + k;
        }
      }
    }
    //check if Tau field is supposed to be double valued and deal with it
    if(myOpts.doubleValuedTau){
      fieldBuffer.resize(faceFieldMap["Tau"].size()/2, 0.0);
      for(int i = 0; i < nFacesPEl; i++){
        if(part == NULL){
          myMesh->getFace2Cell(facesInCell[i], &face2Cells);
        } else {
          if(locFacesInCell[i] != -1){
            myMesh->getFace2Cell(locFacesInCell[i], &face2Cells);
          } else {
            face2Cells[0] = unitCell[0];
          }
        }
        offset = 1;
        if(unitCell[0] == face2Cells[0]){
          offset = 0;
        }
        for(int j = 0; j < nNodesPFc; j++){
          buff = (i*nNodesPFc + j)*matLen;
          for(int k = 0; k < matLen; k++){
            fieldBuffer[buff + k] = faceFieldMap["Tau"][buff*2 + offset*matLen + k];
          }
        }
      }
      faceFieldMap["Tau"] = fieldBuffer;
    }
    //re order the face fields correctly for local ordering
    for(itfm = faceFieldMap.begin(); itfm != faceFieldMap.end(); itfm++){
      int len = itfm->second.size()/(nFacesPEl * nNodesPFc);
      locFieldMap[itfm->first].resize((itfm->second).size());
      for(int i = 0; i < nFacesPEl; i++){
        if(part == NULL){
          myMesh->getFace(facesInCell[i], &face);
        } else {
          if(locFacesInCell[i] != -1){
            myMesh->getFace(locFacesInCell[i], &face);
          } else {
            myMesh->getGhostFace(facesInCell[i], &face);
          }
        }
        offset = i*nNodesPFc;
        for(int j = 0; j <nNodesPFc; j++){
          for(int k = 0; k < len; k++){
            locFieldMap[itfm->first][(offset + j)*len + k] = itfm->second[(offset + face2CellMap[offset + j])*len + k];
          }
        }
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
    if(myOpts.type == IMPLICIT){
      locS = locMat->block(startL, startU, lenL, lenU)*locU + locMat->block(startL, startQ, lenL, lenQ)*locQ + locMat->block(startL, startL, lenL, lenL);
      locS0 = locVec->segment(startL, lenL) - locMat->block(startL, startU, lenL, lenU)*locU0 - locMat->block(startL, startQ, lenL, lenQ)*locQ0;
    } else {
      EMap<const EVector> sol(locFieldMap["Solution"].data(), locFieldMap["Solution"].size());
      EMap<const EVector> flux(locFieldMap["Flux"].data(), locFieldMap["Flux"].size());
      locS0 = locVec->segment(startL, lenL) - locMat->block(startL, startU, lenL, lenU)*sol - locMat->block(startL, startQ, lenL, lenQ)*flux;
      locS = locMat->block(startL, startL, lenL, lenL);
    }
    if((myOpts.type == IMPLICIT) or (myOpts.type == WEXPLICIT)){
      switch(modAssembly->matrix){
        case Add:{
                   linSystem->addValsMatrix(matRowCols, matRowCols, locS.data());
                   break;
                 }
        case Set:{
                   linSystem->setValsMatrix(matRowCols, matRowCols, locS.data());
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
    } else {
      buff = nDOFsPerNode*nNodesPFc;
      for(int i = 0; i < nFacesPEl; i++){
        EMap<EVector> S0map(S0->getValues()->data() + facesInCell[i]*buff, buff);
        for(int j = 0; j < nNodesPFc; j++){
          for(int k = 0; k < nDOFsPerNode; k++){
            switch(modAssembly->rhs){
              case Add:{
                         S0map[face2CellMap[i*nNodesPFc + j]*nDOFsPerNode + k] += locS0[(i*nNodesPFc + j)*nDOFsPerNode + k];
                         break;
                       }
              case Set:{
                         S0map[face2CellMap[i*nNodesPFc + j]*nDOFsPerNode + k] = locS0[(i*nNodesPFc + j)*nDOFsPerNode + k];
                         break;
                       }
            }
          }
        }
      }
      int buff2 = std::pow(buff,2);
      for(int i = 0; i < nFacesPEl; i++){
        EMap<EMatrix> Smap(S->getValues()->data() + facesInCell[i]*buff2, buff, buff);
        for(int j = 0; j < nNodesPFc; j++){
          for(int k = 0; k < nDOFsPerNode; k++){
            for(int l = 0; l < nNodesPFc; l++){
              for(int m = 0; m < nDOFsPerNode; m++){
                switch(modAssembly->rhs){
                  case Add:{
                             Smap(face2CellMap[i*nNodesPFc + j]*nDOFsPerNode + k, face2CellMap[i*nNodesPFc + l]*nDOFsPerNode + m) += locS((i*nNodesPFc + j)*nDOFsPerNode + k, (i*nNodesPFc + l)*nDOFsPerNode + m);
                             break;
                           }
                  case Set:{
                             Smap(face2CellMap[i*nNodesPFc + j]*nDOFsPerNode + k, face2CellMap[i*nNodesPFc + l]*nDOFsPerNode + m) = locS((i*nNodesPFc + j)*nDOFsPerNode + k, (i*nNodesPFc + l)*nDOFsPerNode + m);
                             break;
                           }
                }
              }
            }
          }
        }
      }
    }
    if(verbose){
      pb.update();
    }
  }
  if(myOpts.type != SEXPLICIT){
    linSystem->assemble();
  }
  //boundary
  matRowCols.resize(nNodesPFc*nDOFsPerNode);
  cell.resize(nNodesPFc);
  locCell.resize(nNodesPFc);
  nodes.resize(nNodesPFc);
  for(int iBoundary = 0; iBoundary < bSystem.getBoundaryList()->size(); iBoundary++){
    FEModel * boundaryModel = std::get<0>(bSystem.getBoundaryList()->at(iBoundary));
    modAssembly = boundaryModel->getAssemblyType();
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> boundaryT(
        boundaryModel->getLocalMatrix()->rows(), 
        boundaryModel->getLocalMatrix()->cols());
    const std::set<int> * boundaryFaces = std::get<1>(bSystem.getBoundaryList()->at(iBoundary));
    std::set<int>::const_iterator itFace;
    int locFaceInd;
    int index = 0;
    if(verbose){
      pb.setIterIndex(&index);
      pb.setNumIterations(boundaryFaces->size());
    }
    if(myOpts.type != SEXPLICIT){
      if(modAssembly->matrix == Set){
        if(verbose){
          std::cout << "Assembly - Zero loop (nFaces = " << boundaryFaces->size() << ", nNodes/Face = " << nNodesPFc << "):" << std::endl;
          pb.update();
        }
        std::set<int> zeroSet;
        for(itFace = boundaryFaces->begin(); itFace != boundaryFaces->end(); itFace++){
          std::iota(matRowCols.begin(), matRowCols.end(), (*itFace)*(matRowCols.size()));
          for(int i = 0; i < matRowCols.size(); i++){
            zeroSet.insert(matRowCols[i]);
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
    }
    if(verbose){
      index = 0;
      std::cout << "Assembly - Boundary loop (nFaces = " << boundaryFaces->size() << ", nNodes/Face = " << nNodesPFc << "):" << std::endl;
      pb.update();
    }
    for(itFace = boundaryFaces->begin(); itFace != boundaryFaces->end(); itFace++){
      locFaceInd = *itFace;
      if(part != NULL){
        locFaceInd = part->global2LocalFace(*itFace);
      }
      myMesh->getFace(locFaceInd, &cell);
      unitCell[0] = *itFace;
      if(part == NULL){
        myMesh->getSlicePoints(cell, &nodes);
        constructLocalFields(unitCell, &faceFieldMap);    
      } else {
        locUnitCell[0] = locFaceInd;
        constructLocalFields(unitCell, locUnitCell, &faceFieldMap);
        part->global2LocalNodeSlice(cell, &locCell);
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
      if(myOpts.type != SEXPLICIT){
        boundaryT = *(boundaryModel->getLocalMatrix());
        std::iota(matRowCols.begin(), matRowCols.end(), (*itFace)*(matRowCols.size()));
        switch(modAssembly->matrix){
          case Add:{
                     linSystem->addValsMatrix(matRowCols, matRowCols, boundaryT.data());
                     break;
                   }
          case Set:{
                     linSystem->setValsMatrix(matRowCols, matRowCols, boundaryT.data());
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
      } else {
        buff = nDOFsPerNode*nNodesPFc;
        EMap<EVector> S0map(S0->getValues()->data() + (*itFace)*buff, buff);
        switch(modAssembly->rhs){
          case Add:{
                     S0map += *(boundaryModel->getLocalRHS());
                     break;
                   }
          case Set:{
                     S0map = *(boundaryModel->getLocalRHS());
                     break;
                   }
        }
        int buff2 = std::pow(buff, 2);
        EMap<EMatrix> Smap(S->getValues()->data() + (*itFace)*buff2, buff, buff);
        switch(modAssembly->rhs){
          case Add:{
                     Smap += *(boundaryModel->getLocalMatrix());
                     break;
                   }
          case Set:{
                     Smap = *(boundaryModel->getLocalMatrix());
                     break;
                   }
        }
      }
      if(verbose){
        index += 1;
        pb.update();
      }
    }
    if(myOpts.type != SEXPLICIT){
      linSystem->assemble();
    }
  }
  assembled = 1;
};//assemble

void HDGSolver::solve(){
  if(!assembled){
    throw(ErrorHandle("HDGSolver", "solve", "system must be assembled before solving"));
  }
  Partitioner * part = myMesh->getPartitioner();
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
  std::vector<int> locFacesInCell(nFacesPEl, 0);
  std::vector<int> cell(nNodesPEl, 0);
  std::vector<int> face(nNodesPFc, 0);
  EVector locL(lenL);
  ProgressBar pb;
  if(myOpts.type != SEXPLICIT){
    std::vector<double> * solutionVec = (*fieldMap)["Trace"]->getValues();
    linSystem->solve(solutionVec);
    if(part != NULL){
      std::vector<int> thisRange;
      linSystem->getSolutionOwnership(&thisRange);
      std::vector<int> dofsWeNeed(myMesh->getNumberFaces() * nDOFsPerNode * nNodesPFc, 0);
      Utils::multiplyIndexes(nDOFsPerNode*nNodesPFc, part->getFaceIds(), &dofsWeNeed);
      reorder(solutionVec, &thisRange, &dofsWeNeed);
    }
  } else {
    int lenT2 = std::pow(lenT, 2);
    Eigen::HouseholderQR<EMatrix> invS(lenT, lenT);
    int iFace = 0;
    if(verbose){
      pb.setIterIndex(&iFace);
      pb.setNumIterations(myMesh->getNumberFaces());
      std::cout << "Computing Trace - Face loop (nFaces = " << myMesh->getNumberFaces() << ", nNodes/Fc = " << nNodesPFc << "):" << std::endl;
      pb.update();
    }
    for(iFace = 0; iFace < myMesh->getNumberFaces(); iFace++){
      EMap<EMatrix> SMap(S->getValues()->data() + iFace*lenT2, lenT, lenT);
      invS.compute(SMap);
      EMap<EVector> S0Map(S0->getValues()->data() + iFace*lenT, lenT);
      EMap<EVector> trace(fieldMap->at("Trace")->getValues()->data() + iFace*lenT, lenT);
      trace = invS.solve(S0Map);
      if(verbose){
        pb.update();
      }
    }
  }
  if(part != NULL){
    part->updateSharedInformation();
  }
  double * begining;
  int iEl = 0;
  if(verbose){
    pb.setIterIndex(&iEl);
    pb.setNumIterations(myMesh->getNumberCells());
    std::cout << "Computing Solution and Flux - Element loop (nElements = " << myMesh->getNumberCells() << ", nNodes/El = " << nNodesPEl << "):" << std::endl;
    pb.update();
  }
  for(iEl = 0; iEl < myMesh->getNumberCells(); iEl++){
    myMesh->getCell(iEl, &cell);
    myMesh->getCell2Face(iEl, &facesInCell);
    locFacesInCell = facesInCell;
    if(part != NULL){
      part->global2LocalFaceSlice(facesInCell, &locFacesInCell);
    }
    EMap<EMatrix> locU(U->getValues()->data() + iEl*lenU*lenL, lenU, lenL);
    EMap<EMatrix> locQ(Q->getValues()->data() + iEl*lenQ*lenL, lenQ, lenL);
    EMap<EVector> locU0(U0->getValues()->data() + iEl*lenU, lenU);
    EMap<EVector> locQ0(Q0->getValues()->data() + iEl*lenQ, lenQ);
    for(int i = 0; i < nFacesPEl; i++){
      if(locFacesInCell[i] != -1){
        begining = fieldMap->at("Trace")->getValues()->data()+locFacesInCell[i]*lenT;
        myMesh->getFace(locFacesInCell[i], &face);
      } else {
        begining = fieldMap->at("Trace")->getParValues()->data() + lenT*std::distance(fieldMap->at("Trace")->getParIds()->begin(), std::find(fieldMap->at("Trace")->getParIds()->begin(), fieldMap->at("Trace")->getParIds()->end(), facesInCell[i]));
        myMesh->getGhostFace(facesInCell[i], &face);
      }
      EMap<EVector> buff(begining, lenT);
      for(int j = 0; j < nNodesPFc; j++){
        int found = std::distance(face.begin(), std::find(face.begin(), face.end(), cell[nodeMap->at(i)[j]]));
        for(int ndof = 0; ndof < nDOFsPerNode; ndof++){
          locL[(i*nNodesPFc + j)*nDOFsPerNode + ndof] = buff[found*nDOFsPerNode + ndof];
        }
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
  if(part != NULL){
    part->updateSharedInformation();
  }
};//solve

void HDGSolver::setOptions(HDGSolverOpts opts){
  myOpts = opts;
  setVerbosity(myOpts.verbosity);
};//setOptions

}//hfox
