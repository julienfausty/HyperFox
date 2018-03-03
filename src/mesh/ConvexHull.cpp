#include "ConvexHull.h"

namespace hfox{

void ConvexHull::setVertexes(const std::vector< std::vector<double> > * 
    pvertexes_cand){
  ErrorHandle eh("ConvexHull", "setVertexes");
  eh.checkPointList(*pvertexes_cand);
  this->pvertexes = pvertexes_cand;
  this->dimension = ((*pvertexes)[0]).size();
};

void ConvexHull::setFaces(std::vector< std::vector<int> > * 
    pfaces_cand){
  this->pfaces = pfaces_cand;
};

void ConvexHull::setAlgoType(CHullAlgo typeAlgo_cand){
  this->typeAlgo = typeAlgo_cand;
};

CHullAlgo ConvexHull::getAlgoType() const{
  return this->typeAlgo;
};

int ConvexHull::getDimension() const{
  return this->dimension;
};

void ConvexHull::computeConvexHull(){
  switch(this->typeAlgo){
    case QuickHull:
      this->quickHull();
      break;
    default:
      ErrorHandle eh("ConvexHull","computeConvexHull",
          "the algorithm type given to the ConvexHull class is not yet supported.");
      throw eh;
      break;
  }
};

void ConvexHull::quickHull(){
  std::vector< std::vector<int> > faces;
};

double ConvexHull::calcDistFace(std::vector<int> & face, std::vector<double> & point) 
  const{
    int i;
    EVector n = this->calcNormalVector(face);
    EVector pointv(dimension);
    for(i = 0; i < dimension; i++){
      pointv << point[i] - (*pvertexes)[0][i];
    }
    return n.dot(pointv); 
};

EVector ConvexHull::calcNormalVector(std::vector<int> & face) const{
  int i;
  // construct tangent matrix
  EMatrix tanMat(dimension-1, dimension);
  for(i = 1; i < dimension; i++){
    for(int j = 0; j < dimension; j++){
      tanMat(i-1, j) = (*pvertexes)[face[i]][j] - (*pvertexes)[face[0]][j];
    }
  }
  EVector normal(dimension);
  // devise 1 or -1 symbol
  int symbol;
  if(dimension % 2){
    symbol = 1;
  }else{
    symbol = -1;
  }
  // calc normal component by component using determinants
  for(i = 0; i < dimension; i++){
    EMatrix tempMat(dimension-1, dimension-1);
    tempMat = this->removeColumn(tanMat, i);
    normal(i) = symbol*tempMat.determinant();
    symbol *= -1;
  }
  return normal.normalized();
}

EMatrix ConvexHull::removeColumn(const EMatrix & mat, unsigned int colToRemove) const{
  unsigned int numRows = mat.rows();
  unsigned int numCols = mat.cols()-1;
  EMatrix tempMat = mat; 
  if( colToRemove < numCols ){
    tempMat.block(0,colToRemove,numRows,numCols-colToRemove) = 
      tempMat.rightCols(numCols-colToRemove);
  }
  tempMat.conservativeResize(numRows,numCols);
  return tempMat;
};

} //hfox
