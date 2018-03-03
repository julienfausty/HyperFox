#include "CPolytope.h"

namespace hfox{

CPolytope::CPolytope(){
};

CPolytope::CPolytope(std::vector< std::vector<double> > & points_cand){
  this->setPoints(points_cand);
  this->setFaces();
};

void CPolytope::setFaces(std::vector< std::vector<int> > & faces_cand){
  ErrorHandle eh("CPolytope", "setFaces", 
      "should never try to set faces manually.\n");
  std::cout << "Warning : " << eh.what();
  std::cout << "Computing faces using convex hull algorithm." << std::endl;
  this->setFaces();
};

void CPolytope::setFaces(){
  this->calcConvexHull();
};

void CPolytope::calcConvexHull(){
};

double CPolytope::calcDistFace(int index_face, std::vector<double> & point) 
  const{
  return 1;
};

double CPolytope::calcVolume() const{
  return 1;
};

std::vector<double> CPolytope::calcSurfaces() const{
  std::vector<double> result;
  return result;
};

std::vector< std::vector<double> > CPolytope::calcEdgeLengths() const{
  std::vector< std::vector<double> > result;
  return result;
};

double CPolytope::calcSurface() const{
  return 1;
};

double CPolytope::calcEdgeLength() const{
  return 1;
};

bool CPolytope::isPointInside(const std::vector<double> & point) const{
  return 0;
};

} //hfox
