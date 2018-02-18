#include "Polytope.h"

namespace hfox{

void Polytope::setPoints(std::vector< std::vector<double> > & points_candidate){
  ErrorHandle eh("Polytope", "setPoints");
  try{
    eh.checkPointList(points_candidate);
  }catch(ErrorHandle & e){
    std::cout << e.what() << std::endl;
    exit(0);
  }
  this->dimension = (points[0]).size();
  this->points = points_candidate;
};

void Polytope::setFaces(std::vector< std::vector<int> > & faces_candidate){
  ErrorHandle eh("Polytope", "setFaces");
};

const std::vector< std::vector<double> > * Polytope::getPoints() const{
  const std::vector< std::vector<double> > * ppoints = &points;
  return ppoints;
};

const std::vector< std::vector<int> > * Polytope::getFaces() const{
  const std::vector< std::vector<int> > * pfaces = &faces;
  return pfaces;
};

int Polytope::getNumberPoints() const{
  return (this->points).size();
};

int Polytope::getNumberFaces() const{
  return (this->faces).size();
};

int Polytope::getDimension() const{
  return (this->dimension);
};

} //hfox
