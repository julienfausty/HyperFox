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
  if(S != NULL){delete S;}
  S = new Field(myMesh, Cell, 1, std::pow(nDOFsPerNode*nNodesPFace*nFacesPEl, 2));
  if(S0 != NULL){delete S0;}
  S0 = new Field(myMesh, Cell, 1, nDOFsPerNode*nNodesPFace*nFacesPEl);
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
  calcElementalMatrices();
  applyBoundaryConditions();
  assembleSystem();
  assembled = 1;
};//assemble

void HDGSolver::calcElementalMatrices(){
  std::fill(S->getValues()->begin(), S->getValues()->end(), 0.0);
  std::fill(S0->getValues()->begin(), S0->getValues()->end(), 0.0);
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
  std::vector<double> fieldBuffer;
  std::vector<int>::iterator intIt;
  int buff = 0;
  int offset = 0;
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
    std::cout << "Local calculations - Element loop (nElements = " << myMesh->getNumberCells() << ", nNodes/El = " << nNodesPEl << "):" << std::endl;
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
          throw(ErrorHandle("HDGSolver", "calcElementalMatrices", "couldn't find cell node in face."));
        }
        face2CellMap[i*nNodesPFc + j] = std::distance(face.begin(), intIt);
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
    EMap<EMatrix> locS(S->getValues()->data() + iEl*lenL*lenL, lenL, lenL);
    EMap<EVector> locU0(U0->getValues()->data() + iEl*lenU, lenU);
    EMap<EVector> locQ0(Q0->getValues()->data() + iEl*lenQ, lenQ);
    EMap<EVector> locS0(S0->getValues()->data() + iEl*lenL, lenL);
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
    if(verbose){
      pb.update();
    }
  }
};//calcElementalMatrices

void HDGSolver::applyBoundaryConditions(){
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
  int lenT = nNodesPFc*nDOFsPerNode;
  int lenL = nFacesPEl*nNodesPFc*nDOFsPerNode;
  const AssemblyType * modAssembly;
  std::map<std::string, std::vector<double> > locFieldMap;
  std::map<std::string, std::vector<double> > faceFieldMap = prepareLocalFieldMap(Face);
  std::map<std::string, std::vector<double> > cellFieldMap = prepareLocalFieldMap(Cell);
  std::map<std::string, std::vector<double> > nodalFieldMap = prepareLocalFieldMap(Node);
  std::map<std::string, std::vector<double> >::iterator itfm;
  std::vector<int> cell(nNodesPEl);
  std::vector<int> locCell(nNodesPEl);
  std::vector<int> face(nNodesPFc);
  std::vector<int> locFace(nNodesPFc);
  std::vector<int> face2Cells(2, 0);
  std::vector<int> cell2Face(nFacesPEl, 0);
  std::vector<int> unitCell(1, 0);
  std::vector<int> locUnitCell(1, 0);
  std::vector< std::vector<double> > nodes(nNodesPFc, std::vector<double>(dim, 0.0));
  std::vector<int> cell2FaceMap(nNodesPFc);
  std::set<int>::const_iterator itFace;
  std::vector<double> fieldBuffer;
  int locCellInd;
  int locFaceInd;
  int locFaceInEl;
  int index = 0;
  ProgressBar pb;
  for(int iBoundary = 0; iBoundary < bSystem.getBoundaryList()->size(); iBoundary++){
    BoundaryModel * boundaryModel = std::get<0>(bSystem.getBoundaryList()->at(iBoundary));
    modAssembly = boundaryModel->getAssemblyType();
    const std::set<int> * boundaryFaces = std::get<1>(bSystem.getBoundaryList()->at(iBoundary));
    if(verbose){
      index = 0;
      std::cout << "Assembly - Boundary loop " << iBoundary << " (nFaces = " << boundaryFaces->size() << ", nNodes/Face = " << nNodesPFc << "):" << std::endl;
      pb.setIterIndex(&index);
      pb.setNumIterations(boundaryFaces->size());
      pb.update();
    }
    for(itFace = boundaryFaces->begin(); itFace != boundaryFaces->end(); itFace++){
      locFaceInd = *itFace;
      if(part != NULL){
        locFaceInd = part->global2LocalFace(*itFace);
      }
      myMesh->getFace(locFaceInd, &face);
      myMesh->getFace2Cell(locFaceInd, &face2Cells);
      for(int iEl = 0; iEl < face2Cells.size(); iEl++){      
        unitCell[0] = *itFace;
        if(part == NULL){
          locCellInd = face2Cells[iEl];
          myMesh->getSlicePoints(face, &nodes);
          constructLocalFields(unitCell, &faceFieldMap);
          constructLocalFields(face, &nodalFieldMap);
          unitCell[0] = face2Cells[iEl];
          locUnitCell[0] = unitCell[0];
          constructLocalFields(unitCell, &cellFieldMap);
          myMesh->getCell(face2Cells[iEl], &cell);
        } else {
          locUnitCell[0] = locFaceInd;
          constructLocalFields(unitCell, locUnitCell, &faceFieldMap);
          part->global2LocalNodeSlice(face, &locFace);
          for(int i = 0; i < locFace.size(); i++){
            if(locFace[i] != -1){
              myMesh->getPoint(locFace[i], &(nodes[i]));
            } else {
              myMesh->getGhostPoint(face[i], &(nodes[i]));
            }
          }
          constructLocalFields(face, locFace, &nodalFieldMap);
          locCellInd = part->global2LocalElement(face2Cells[iEl]);
          unitCell[0] = face2Cells[iEl];
          locUnitCell[0] = locCellInd;
          constructLocalFields(unitCell, locCell, &cellFieldMap);
          if(locCellInd != -1){
            myMesh->getCell(locCellInd, &cell);
          } else {
            myMesh->getGhostCell(face2Cells[iEl], &cell);
          }
        }
        boundaryModel->setElementNodes(&nodes);
        if(myOpts.doubleValuedTau){
          int matLen = std::pow(nDOFsPerNode, 2);
          fieldBuffer.resize(faceFieldMap["Tau"].size()/2, 0.0);
          for(int j = 0; j < nNodesPFc; j++){
            int buff = j*matLen;
            for(int k = 0; k < matLen; k++){
              fieldBuffer[buff + k] = faceFieldMap["Tau"][buff*2 + iEl*matLen + k];
            }
          }
          faceFieldMap["Tau"] = fieldBuffer;
        }
        for(itfm = nodalFieldMap.begin(); itfm != nodalFieldMap.end(); itfm++){locFieldMap[itfm->first] = itfm->second;}
        for(itfm = faceFieldMap.begin(); itfm != faceFieldMap.end(); itfm++){locFieldMap[itfm->first] = itfm->second;}
        for(int i = 0; i < nNodesPFc; i++){
          cell2FaceMap[i] = std::distance(cell.begin(), (std::find(cell.begin(), cell.end(), face[i])));
        }
        for(itfm = cellFieldMap.begin(); itfm != cellFieldMap.end(); itfm++){
          int len = itfm->second.size()/(nNodesPEl);
          locFieldMap[itfm->first].resize(nNodesPFc*len);
          for(int i = 0; i < nNodesPFc; i++){
            std::copy(itfm->second.begin() + cell2FaceMap[i]*len, itfm->second.begin() + (cell2FaceMap[i]+1)*len, locFieldMap[itfm->first].begin() + i*len);
          }
        }
        boundaryModel->setFieldMap(&locFieldMap);
        boundaryModel->compute();
        //put in S
        //FROM HERE ON LOCCELLIND MUST NOT BE -1: COMMUNICATION NECESSARY IF IT IS THE CASE (BOUNDARY FACE CANNOT ALSO BE PARTITION BOUNDARY FACES)
        EMap<EMatrix> locU(U->getValues()->data() + locCellInd*lenU*lenL, lenU, lenL);
        EMap<EMatrix> locQ(Q->getValues()->data() + locCellInd*lenQ*lenL, lenQ, lenL);
        EMap<EMatrix> locS(S->getValues()->data() + locCellInd*lenL*lenL, lenL, lenL);
        EMap<EVector> locU0(U0->getValues()->data() + locCellInd*lenU, lenU);
        EMap<EVector> locQ0(Q0->getValues()->data() + locCellInd*lenQ, lenQ);
        EMap<EVector> locS0(S0->getValues()->data() + locCellInd*lenL, lenL);
        //position face in cell
        myMesh->getCell2Face(locCellInd, &cell2Face);
        locFaceInEl = std::distance(cell2Face.begin(), std::find(cell2Face.begin(), cell2Face.end(), *itFace));
        //Zero out rows
        if(modAssembly->matrix == Set){
          locS.block(lenT*locFaceInEl, 0, lenT, lenL) = EMatrix::Zero(lenT, lenL);
        }
        if(modAssembly->rhs == Set){
          locS0.segment(lenT*locFaceInEl, lenT) = EVector::Zero(lenT);
        }
        if(boundaryModel->getBoundaryModelType() == CGType){
          locS.block(lenT*locFaceInEl, lenT*locFaceInEl, lenT, lenT) += *(boundaryModel->getLocalMatrix());
          locS0.segment(lenT*locFaceInEl, lenT) += *(boundaryModel->getLocalRHS());
        } else if(boundaryModel->getBoundaryModelType() == HDGType){
          locCellInd = face2Cells[iEl];
          if(part != NULL){
            locCellInd = part->global2LocalElement(face2Cells[iEl]);
          }
          EMatrix locUFc(lenT, lenL);
          EMatrix locQFc(dim*lenT, lenL);
          EVector locU0Fc(lenT);
          EVector locQ0Fc(lenT*dim);
          for(int iN = 0; iN < nNodesPFc; iN++){
            locUFc.block(iN*nDOFsPerNode, 0, nDOFsPerNode, lenL) = locU.block(cell2FaceMap[iN]*nDOFsPerNode, 0, nDOFsPerNode, lenL);
            locQFc.block(iN*nDOFsPerNode*dim, 0, nDOFsPerNode*dim, lenL) = locQ.block(cell2FaceMap[iN]*nDOFsPerNode*dim, 0, nDOFsPerNode*dim, lenL);
            locU0Fc.segment(iN*nDOFsPerNode, nDOFsPerNode) = locU0.segment(cell2FaceMap[iN]*nDOFsPerNode, nDOFsPerNode);
            locQ0Fc.segment(iN*nDOFsPerNode*dim, nDOFsPerNode*dim) = locQ0.segment(cell2FaceMap[iN]*nDOFsPerNode*dim, nDOFsPerNode*dim);
          }
          locS.block(lenT*locFaceInEl, 0, lenT, lenL) += boundaryModel->getLocalMatrix()->block(0, 0, lenT, lenT) *locUFc  + 
            boundaryModel->getLocalMatrix()->block(0, lenT, lenT, lenT*dim) * locQFc + 
            boundaryModel->getLocalMatrix()->block(0, lenT*(1+dim), lenT, lenT);
          locS0.segment(lenT*locFaceInEl, lenT) += *(boundaryModel->getLocalRHS()) - 
            boundaryModel->getLocalMatrix()->block(0, 0, lenT, lenT)*locU0Fc - 
            boundaryModel->getLocalMatrix()->block(0, lenT, lenT, lenT*dim) * locQ0Fc;
        } else {
          throw(ErrorHandle("HDGSolver", "applyBoundaryConditions", "the boundary type must be CG or HDG"));
        }
      }
    }
    if(verbose){
      index += 1;
      pb.update();
    }
  }
};//applyBoundaryConditions

void HDGSolver::assembleSystem(){
  if(myOpts.type != SEXPLICIT){
    if(assembled){
      linSystem->clearSystem();
    }
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
  int lenT = nNodesPFc*nDOFsPerNode;
  int lenL = nFacesPEl*nNodesPFc*nDOFsPerNode;
  const AssemblyType * modAssembly = model->getAssemblyType();
  const std::vector< std::vector<int> > * nodeMap = refEl->getFaceNodes();
  std::vector<int> cell(nNodesPEl);
  std::vector<int> face(nNodesPFc);
  std::vector<int> facesInCell(nFacesPEl, 0);
  std::vector<int> locFacesInCell(nFacesPEl, 0);
  std::vector<int> face2CellMap(nFacesPEl*nNodesPFc, 0);
  Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> locS(lenL, lenL);
  EVector locS0(lenL);
  std::vector<int> matRowCols(lenL, 0);
  std::vector<int>::iterator intIt;
  int buff;
  int iEl = 0;
  ProgressBar pb;
  if(verbose){
    pb.setIterIndex(&iEl);
    pb.setNumIterations(myMesh->getNumberCells());
    std::cout << "Assembly - Element loop (nElements = " << myMesh->getNumberCells() << ", nNodes/El = " << nNodesPEl << "):" << std::endl;
    pb.update();
  }
  for(iEl = 0; iEl < myMesh->getNumberCells(); iEl++){ 
    //compute indexes
    myMesh->getCell(iEl, &cell);
    myMesh->getCell2Face(iEl, &facesInCell);
    locFacesInCell = facesInCell;
    if(part != NULL){
      part->global2LocalFaceSlice(facesInCell, &locFacesInCell);
    }
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
          throw(ErrorHandle("HDGSolver", "assembleSystem", "couldn't find cell node in face."));
        }
        face2CellMap[i*nNodesPFc + j] = std::distance(face.begin(), intIt);
        for(int k = 0; k < nDOFsPerNode; k++){
          matRowCols[(i*nNodesPFc + j)*nDOFsPerNode + k] = (facesInCell[i] * nNodesPFc + face2CellMap[i*nNodesPFc + j])*nDOFsPerNode + k;
        }
      }
    }
    //set correct data format
    EMap<EMatrix> locSStored(S->getValues()->data() + iEl*lenL*lenL, lenL, lenL);
    locS = locSStored;
    EMap<EVector> locS0Stored(S0->getValues()->data() + iEl*lenL, lenL);
    locS0 = locS0Stored;
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
};//assembleSystem

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
