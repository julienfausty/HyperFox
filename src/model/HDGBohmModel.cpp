#include "HDGBohmModel.h"

namespace hfox{

namespace nGamma{

void HDGBohmModel::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator it;
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "b"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("b") == fieldMap.end()){
    throw(ErrorHandle("HDGBohmModel", "setFieldMap", 
          "need to give a field named b to the HDGBohmModel"));
  }
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "SoundVelocity"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("SoundVelocity") == fieldMap.end()){
    throw(ErrorHandle("HDGBohmModel", "setFieldMap", 
          "need to give a field named SoundVelocity to the HDGBohmModel"));
  }
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "D"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("D") == fieldMap.end()){
    throw(ErrorHandle("HDGBohmModel", "setFieldMap", 
          "need to give a field named D to the HDGBohmModel"));
  }
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "G"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("G") == fieldMap.end()){
    throw(ErrorHandle("HDGBohmModel", "setFieldMap", 
          "need to give a field named G to the HDGBohmModel"));
  }
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "ExteriorNormals"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("ExteriorNormals") == fieldMap.end()){
    throw(ErrorHandle("HDGBohmModel", "setFieldMap", 
          "need to give a field named ExteriorNormals to the HDGBohmModel"));
  }
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "Solution"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("Solution") == fieldMap.end()){
    throw(ErrorHandle("HDGBohmModel", "setFieldMap", 
          "need to give a field named Solution to the HDGBohmModel"));
  }
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "Trace"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("Trace") == fieldMap.end()){
    throw(ErrorHandle("HDGBohmModel", "setFieldMap", 
          "need to give a field named Trace to the HDGBohmModel"));
  }
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "Flux"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("Flux") == fieldMap.end()){
    throw(ErrorHandle("HDGBohmModel", "setFieldMap", 
          "need to give a field named Flux to the HDGBohmModel"));
  }
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "Tau"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("Tau") == fieldMap.end()){
    throw(ErrorHandle("HDGBohmModel", "setFieldMap", 
          "need to give a field named Tau to the HDGBohmModel"));
  }
  fieldSet = 1;
};//setFieldMap

void HDGBohmModel::allocate(int nDOFsPerNode){
  if(nDOFsPerNode != 2){
    throw(ErrorHandle("HDGBohmModel", "allocate", 
          "the number of degrees of freedom per node must be 2"));
  }
  assembly.matrix = Set;
  assembly.rhs = Set;
  myType = HDGType;
  initializeOperators();
  localMatrix = EMatrix::Zero(nDOFsPerNode*refEl->getNumNodes(), nDOFsPerNode*refEl->getNumNodes()*4);
  originalSystem = localMatrix;
  localRHS = EVector::Zero(nDOFsPerNode*refEl->getNumNodes());
  for(std::map<std::string, Operator*>::iterator it = operatorMap.begin(); it != operatorMap.end(); it++){
    (it->second)->allocate(1);
  }
  allocated = 1;
};//allocate

void HDGBohmModel::initializeOperators(){
};//initializeOperators

void HDGBohmModel::computeLocalMatrix(){
  if(!(fieldSet and allocated)){
    throw(ErrorHandle("HDGBohmModel", "computeLocalMatrix", 
          "the fields should be set and the object allocated before computing."));
  }
  if(!transferFunction){
    throw(ErrorHandle("HDGBohmModel", "computeLocalMatrix", 
          "the transfer function should be set before computing."));
  }
  jacobians = Operator::calcJacobians(*elementNodes, refEl);
  invJacobians = Operator::calcInvJacobians(jacobians);
  dV = Operator::calcMeasure(Operator::calcDetJacobians(jacobians), refEl);
  localMatrix = EMatrix::Zero(localMatrix.rows(), localMatrix.cols());
  originalSystem = EMatrix::Zero(originalSystem.rows(), originalSystem.cols());
  //refEl data
  int nNodes = refEl->getNumNodes();
  int nIPs = refEl->getNumIPs();
  const std::vector< std::vector<double> > * shapes = refEl->getIPShapeFunctions();
  //buffers
  double dbuff;
  std::vector<EVector> normals(nIPs, EVector::Zero(2));
  std::vector<EVector> traces(nIPs, EVector::Zero(2));
  std::vector<EMatrix> fluxes(nIPs, EMatrix::Zero(2, 2));
  std::vector<EMatrix> taus(nIPs, EMatrix::Zero(2, 2));
  std::vector<EMatrix> Ds(nIPs, EMatrix::Zero(2, 2));
  std::vector<EMatrix> Gs(nIPs, EMatrix::Zero(2, 2));
  std::vector<EVector> solution(nIPs, EVector::Zero(2));
  std::vector<double> ipSoundSpeed(nIPs, 0);
  std::vector<double> bdotN(nIPs, 0);
  double tfInput0, tfInput1;
  double dirichlet, neumann;
  std::vector<double> transferVals(nIPs, 0);
  std::vector< std::vector<double> > transferDerivVals(nIPs, std::vector<double>(2, 0));
  EMatrix normalMat = EMatrix::Zero(2, nNodes);
  EMatrix bMat = EMatrix::Zero(2, nNodes);
  EMatrix tMat = EMatrix::Zero(2, nNodes);
  EMatrix qMat = EMatrix::Zero(4, nNodes);
  EMatrix tauMat = EMatrix::Zero(4, nNodes);
  EMatrix dMat = EMatrix::Zero(4, nNodes);
  EMatrix gMat = EMatrix::Zero(4, nNodes);
  EMatrix uMat = EMatrix::Zero(2, nNodes);
  EMap<const EVector> soundSpeedVec(fieldMap["SoundVelocity"]->data(), fieldMap["SoundVelocity"]->size());
  EVector vecBuff = EVector::Zero(2);
  EMatrix matBuff;
  int dSize = fieldMap["D"]->size()/nNodes;
  int gSize = fieldMap["G"]->size()/nNodes;
  EVector idVec(4);
  idVec << 1.0, 0.0, 0.0, 1.0;
  //setup
  for(int ndof = 0; ndof < nNodes; ndof++){
    normalMat.col(ndof) = EMap<const EVector>(fieldMap["ExteriorNormals"]->data() + ndof*2, 2);
    bMat.col(ndof) = EMap<const EVector>(fieldMap["b"]->data() + ndof*2, 2);
    tMat.col(ndof) = EMap<const EVector>(fieldMap["Trace"]->data() + ndof*2, 2);
    qMat.col(ndof) = EMap<const EVector>(fieldMap["Flux"]->data() + ndof*4, 4);
    tauMat.col(ndof) = EMap<const EVector>(fieldMap["Tau"]->data() + ndof*4, 4);
    if(dSize == 1){
      dMat.col(ndof) = idVec * fieldMap["D"]->at(ndof);
    } else {
      dMat.col(ndof) = EMap<const EVector>(fieldMap["D"]->data() + ndof*4, 4);
    }
    if(gSize == 1){
      gMat.col(ndof) = idVec * fieldMap["G"]->at(ndof);
    } else {
      gMat.col(ndof) = EMap<const EVector>(fieldMap["G"]->data() + ndof*4, 4);
    }
    uMat.col(ndof) = EMap<const EVector>(fieldMap["Solution"]->data() + ndof*2, 2);
  }
  for(int ip = 0; ip < nIPs; ip++){
    EMap<const EVector> shapeVec(shapes->at(ip).data(), shapes->at(ip).size());
    normals[ip] = normalMat * shapeVec;
    ipSoundSpeed[ip] = shapeVec.dot(soundSpeedVec);
    vecBuff = bMat * shapeVec;
    bdotN[ip] = vecBuff.dot(normals[ip]);
    traces[ip] = tMat*shapeVec;
    matBuff = qMat*shapeVec;
    fluxes[ip] = EMap<EMatrix>((matBuff).data(), 2, 2);
    matBuff = tauMat*shapeVec;
    taus[ip] = EMap<EMatrix>((matBuff).data(), 2, 2);
    matBuff = dMat*shapeVec;
    Ds[ip] = EMap<EMatrix>((matBuff).data(), 2, 2);
    matBuff = gMat*shapeVec;
    Gs[ip] = EMap<EMatrix>((matBuff).data(), 2, 2);
    solution[ip] = uMat*shapeVec;
    tfInput0 = solution[ip][1]*bdotN[ip];
    tfInput1 = solution[ip][0]*ipSoundSpeed[ip]*bdotN[ip];
    transferVals[ip] = transferFunction(tfInput0, tfInput1);
    transferDerivVals[ip] = derivTransferFunction(tfInput0, tfInput1);
  }
  //computation
  //n condition
  for(int ip = 0; ip < refEl->getNumIPs(); ip++){
    vecBuff = Ds[ip]*normals[ip];
    for(int ndof = 0; ndof < nNodes; ndof++){
      for(int mdof = 0; mdof < nNodes; mdof++){
        dbuff = dV[ip]*(shapes->at(ip)[ndof])*(shapes->at(ip)[mdof]);
        originalSystem(mdof*2, ndof*2) += taus[ip](0,0)*dbuff;
        originalSystem(mdof*2, ndof*2 + 1) += taus[ip](0,1)*dbuff;
        originalSystem(mdof*2, nNodes*6 + ndof*2) -= taus[ip](0,0)*dbuff;
        originalSystem(mdof*2, nNodes*6 + ndof*2 + 1) -= taus[ip](0,1)*dbuff;
        for(int d = 0; d < 2; d++){
          originalSystem(mdof*2, 2*nNodes + (ndof*2 + d)*2) -= dbuff*vecBuff[d];
        }
      }
    }
  }
  //Gamma condition
  for(int ip = 0; ip < nIPs; ip++){
    vecBuff = Gs[ip]*normals[ip];
    for(int ndof = 0; ndof < nNodes; ndof++){
      for(int mdof = 0; mdof < nNodes; mdof++){
        dbuff = dV[ip]*(shapes->at(ip)[ndof])*(shapes->at(ip)[mdof]);
        //neumann
        originalSystem(mdof*2 + 1, ndof*2) += taus[ip](1, 0)*(1-transferVals[ip])*dbuff;
        originalSystem(mdof*2 + 1, ndof*2 + 1) += taus[ip](1, 1)*(1-transferVals[ip])*dbuff;
        originalSystem(mdof*2 + 1, nNodes*6 + ndof*2) -= taus[ip](1, 0)*(1-transferVals[ip])*dbuff;
        originalSystem(mdof*2 + 1, nNodes*6 + ndof*2 + 1) -= taus[ip](1, 1)*(1-transferVals[ip])*dbuff;
        for(int d = 0; d < 2; d++){
          originalSystem(mdof*2 + 1, 2*nNodes + (ndof*2 + d)*2 + 1) -= (1-transferVals[ip])*dbuff*vecBuff[d];
        }
        //dirichlet
        originalSystem(mdof*2 + 1, ndof*2) -= transferVals[ip]*dbuff*ipSoundSpeed[ip]*std::abs(bdotN[ip]);
        originalSystem(mdof*2 + 1, ndof*2 + 1) += transferVals[ip]*dbuff*bdotN[ip];
        //mixed
        dirichlet = dbuff*(solution[ip][1]*bdotN[ip] - solution[ip][0]*ipSoundSpeed[ip]*std::abs(bdotN[ip]));
        neumann = dbuff*(-vecBuff.dot(fluxes[ip].row(1)) + taus[ip](1,1)*(solution[ip][1] - traces[ip][1]) + taus[ip](1, 0)*(solution[ip][0] - traces[ip][0]));
        localMatrix(mdof*2 + 1, ndof*2) += (ipSoundSpeed[ip]*bdotN[ip]*transferDerivVals[ip][1])*(dirichlet - neumann);
        localMatrix(mdof*2 + 1, ndof*2 + 1) += (bdotN[ip]*transferDerivVals[ip][0])*(dirichlet - neumann);
      }
    }
  }
  //add in original system to matrix
  localMatrix += originalSystem;
};//computeLocalMatrix


void HDGBohmModel::computeLocalRHS(){
  if(!(allocated and fieldSet)){
    throw(ErrorHandle("HDGBohmModel", "computeLocalRHS", 
          "the fields should be set and the object allocated before computing."));
  }
  if(!transferFunction){
    throw(ErrorHandle("HDGBohmModel", "computeLocalRHS", 
          "the transfer function should be set before computing."));
  }
  EVector U0 = EVector::Zero(localMatrix.cols());
  int nNodes = refEl->getNumNodes();
  U0.segment(0, nNodes*2) = EMap<const EVector>(fieldMap["Solution"]->data(), fieldMap["Solution"]->size());
  U0.segment(nNodes*2, nNodes*4) = EMap<const EVector>(fieldMap["Flux"]->data(), fieldMap["Flux"]->size());
  U0.segment(nNodes*6, nNodes*2) = EMap<const EVector>(fieldMap["Trace"]->data(), fieldMap["Trace"]->size());
  localRHS = (localMatrix - originalSystem)*U0;
};//computeLocalRHS

};//nGamma

};//hfox
