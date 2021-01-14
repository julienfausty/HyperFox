#include "BohmModel.h"

namespace hfox{

namespace nGamma{

void BohmModel::setFieldMap(const std::map<std::string, std::vector<double> > * fm){
  std::map<std::string, std::vector<double> >::const_iterator it;
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "b"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("b") == fieldMap.end()){
    throw(ErrorHandle("BohmModel", "setFieldMap", 
          "need to give a field named b to the BohmModel"));
  }
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "SoundVelocity"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("SoundVelocity") == fieldMap.end()){
    throw(ErrorHandle("BohmModel", "setFieldMap", 
          "need to give a field named SoundVelocity to the BohmModel"));
  }
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "ExteriorNormals"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("ExteriorNormals") == fieldMap.end()){
    throw(ErrorHandle("BohmModel", "setFieldMap", 
          "need to give a field named ExteriorNormals to the BohmModel"));
  }
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "Trace"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("Trace") == fieldMap.end()){
    throw(ErrorHandle("BohmModel", "setFieldMap", 
          "need to give a field named Trace to the BohmModel"));
  }
  for(it = fm->begin(); it != fm->end(); it++){
    if(it->first == "Flux"){
      fieldMap[it->first] = &(it->second);
    }
  }
  if(fieldMap.find("Flux") == fieldMap.end()){
    throw(ErrorHandle("BohmModel", "setFieldMap", 
          "need to give a field named Flux to the BohmModel"));
  }
  fieldSet = 1;
};//setFieldMap

void BohmModel::allocate(int nDOFsPerNode){
  assembly.matrix = Set;
  assembly.rhs = Set;
  initializeOperators();
  localMatrix = EMatrix::Zero(nDOFsPerNode*refEl->getNumNodes(), nDOFsPerNode*refEl->getNumNodes()*4);
  originalSystem = localMatrix;
  localRHS = EVector::Zero(nDOFsPerNode*refEl->getNumNodes());
  for(std::map<std::string, Operator*>::iterator it = operatorMap.begin(); it != operatorMap.end(); it++){
    (it->second)->allocate(1);
  }
  allocated = 1;
};//allocate

void BohmModel::initializeOperators(){
};//initializeOperators

void BohmModel::computeLocalMatrix(){
  if(!(fieldSet and allocated)){
    throw(ErrorHandle("BohmModel", "computeLocalMatrix", 
          "the fields should be set and the object allocated before computing."));
  }
  if(!transferFunction){
    throw(ErrorHandle("BohmModel", "computeLocalMatrix", 
          "the transfer function should be set before computing."));
  }
  jacobians = Operator::calcJacobians(*elementNodes, refEl);
  invJacobians = Operator::calcInvJacobians(jacobians);
  dV = Operator::calcMeasure(Operator::calcDetJacobians(jacobians), refEl);
  localMatrix = EMatrix::Zero(localMatrix.rows(), localMatrix.cols());
  //refEl data
  int nNodes = refEl->getNumNodes();
  int nIPs = refEl->getNumIPs();
  const std::vector< std::vector<double> > * shapes = refEl->getIPShapeFunctions();
  //buffers
  double dbuff;
  std::vector<EVector> normals(nIPs, EVector::Zero(2));
  std::vector<EVector> traces(nIPs, EVector::Zero(2));
  std::vector<EVector> fluxes(nIPs, EMatrix::Zero(2, 2));
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
  EMap<const EVector> soundSpeedVec(fieldMap["SoundSpeed"]->data(), fieldMap["SoundSpeed"]->size());
  EVector vecBuff = EVector::Zero(2);
  EMatrix matBuff;
  //setup
  for(int ndof = 0; ndof < nNodes; ndof++){
    normalMat.col(ndof) = EMap<const EVector>(fieldMap["ExteriorNormals"]->data() + ndof*2, 2);
    bMat.col(ndof) = EMap<const EVector>(fieldMap["b"]->data() + ndof*2, 2);
    tMat.col(ndof) = EMap<const EVector>(fieldMap["Trace"]->data() + ndof*2, 2);
    qMat.col(ndof) = EMap<const EVector>(fieldMap["Flux"]->data() + ndof*4, 4);
  }
  for(int ip = 0; ip < refEl->getNumIPs(); ip++){
    EMap<const EVector> shapeVec(shapes->at(ip).data(), shapes->at(ip).size());
    normals[ip] = normalMat * shapeVec;
    ipSoundSpeed[ip] = shapeVec.dot(soundSpeedVec);
    vecBuff = bMat * shapeVec;
    bdotN[ip] = vecBuff.dot(normals[ip]);
    traces[ip] = tMat*shapeVec;
    matBuff = qMat*shapeVec;
    fluxes[ip] = EMap<EMatrix>((matBuff).data(), 2, 2);
    tfInput0 = traces[ip][1]*bdotN[ip];
    tfInput1 = traces[ip][0]*ipSoundSpeed[ip]*bdotN[ip];
    transferVals[ip] = transferFunction(tfInput0, tfInput1);
    transferDerivVals[ip] = derivTransferFunction(tfInput0, tfInput1);
  }
  //computation
  //n condition
  for(int ip = 0; ip < refEl->getNumIPs(); ip++){
    for(int ndof = 0; ndof < nNodes; ndof++){
      for(int mdof = 0; mdof < nNodes; mdof++){
        dbuff = dV[ip]*(shapes->at(ip)[ndof])*(shapes->at(ip)[mdof]);
        for(int d = 0; d < 2; d++){
          originalSystem(mdof*2, nNodes + (ndof*2 + d)*2) += dbuff*normals[ip][d];
        }
      }
    }
  }
  //Gamma condition
  for(int ip = 0; ip < nIPs; ip++){
    for(int ndof = 0; ndof < nNodes; ndof++){
      for(int mdof = 0; mdof < nNodes; mdof++){
        dbuff = dV[ip]*(shapes->at(ip)[ndof])*(shapes->at(ip)[mdof]);
        //neumann
        for(int d = 0; d < 2; d++){
          originalSystem(mdof*2 + 1, nNodes + (ndof*2 + d)*2 + 1) += (1-transferVals[ip])*dbuff*normals[ip][d];
        }
        //dirichlet
        originalSystem(mdof*2 + 1, nNodes*6 + ndof*2) -= transferVals[ip]*dbuff*ipSoundSpeed[ip]*std::abs(bdotN[ip]);
        originalSystem(mdof*2 + 1, nNodes*6 + ndof*2 + 1) += transferVals[ip]*dbuff*bdotN[ip];
        //add in original system to matrix
        localMatrix += originalSystem;
        //mixed
        dirichlet = dbuff*(traces[ip][1]*bdotN[ip] - traces[ip][0]*ipSoundSpeed[ip]*std::abs(bdotN[ip]));
        neumann = dbuff*normals[ip].dot(fluxes[ip].row(1));
        localMatrix(mdof*2 + 1, nNodes*6 + ndof*2) += (ipSoundSpeed[ip]*bdotN[ip]*transferDerivVals[ip][1])*(dirichlet - neumann);
        localMatrix(mdof*2 + 1, nNodes*6 + ndof*2 + 1) += (bdotN[ip]*transferDerivVals[ip][0])*(dirichlet - neumann);
      }
    }
  }
};//computeLocalMatrix


void BohmModel::computeLocalRHS(){
  if(!(allocated and fieldSet)){
    throw(ErrorHandle("BohmModel", "computeLocalRHS", 
          "the fields should be set and the object allocated before computing."));
  }
  if(!transferFunction){
    throw(ErrorHandle("BohmModel", "computeLocalRHS", 
          "the transfer function should be set before computing."));
  }
  EVector U0 = EVector::Zero(localMatrix.cols());
  int nNodes = refEl->getNumNodes();
  U0.segment(nNodes*2, nNodes*4) = EMap<const EVector>(fieldMap["Flux"]->data(), fieldMap["Flux"]->size());
  U0.segment(nNodes*6, nNodes*2) = EMap<const EVector>(fieldMap["Trace"]->data(), fieldMap["Trace"]->size());
  localRHS = (localMatrix - originalSystem)*U0;
};//computeLocalRHS

};//nGamma

};//hfox