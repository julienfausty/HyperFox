#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <vector>
#include <string>
#include <iostream>
#include "CPolytope.h"
#include "ErrorHandle.h"

/*!
 * \brief a convex polytope with the dimension of space + 1 vetices.
 */

namespace hfox{

template<int d>
class Simplex : public CPolytope{
  public:
    // Constructors
    /*!
     * \brief empty constructor.
     */
    Simplex(){dimension = d;};
    /*!
     * \brief points constructor.
     */
    Simplex(std::vector< std::vector<double> > & points_candidate);
    /*!
     * \brief set points in Simplex.
     */
    void setPoints(std::vector< std::vector<double> > & points_candidate);
    /*!
     * \brief set faces in Simplex
     */
    void setFaces();
};

template<int d>
void Simplex<d>::setPoints(std::vector< std::vector<double> > & 
    points_candidate){
    CPolytope::setPoints(points_candidate);
    if(points.size() != dimension+1){
      ErrorHandle eh("Simplex", "PointsConstructor", 
          "The number of points ("+std::to_string(points.size())+
          ") is not consistent with the dimension of the Simplex ("
          +std::to_string(dimension)+")." );
      std::cout << eh.what() << std::endl;
      std::cerr << eh.what() << std::endl;
      exit(0);
    }
};

template<int d>
void Simplex<d>::setFaces(){
  this->faces.resize(d+1);
  std::vector<int> * pface;
  int i, j;
  for(i = 0; i < (d+1); i+=2){
    int delim = 0;
    for(j = 0; j < (d+1); j++){
      pface = &(this->faces[i]);
      pface->resize(d);
      if(i != j){
        (*pface)[j+delim] = i;
      }else{
        delim = -1;
      }
    }
  }
  for(i = 1; i < (d+1); i+=2){
    int delim = -1;
    for(j = 0; j < (d+1); j++){
      pface = &(this->faces[i]);
      pface->resize(d);
      if(i != j){
        (*pface)[(d+1)-(j+delim)] = i;
      }else{
        delim = 0;
      }
    }
  }
}

template<int d>
Simplex<d>::Simplex(std::vector< std::vector<double> > & points_candidate){
    dimension = d;
    this->setPoints(points_candidate);
    this->setFaces();
};

} //hfox

#endif //SIMPLEX_H
