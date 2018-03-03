#ifndef SIMPLEX_H
#define SIMPLEX_H

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include "CPolytope.h"
#include "ErrorHandle.h"
#include "DenseEigen.h"

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
    Simplex(){dimension = d;this->setFaces();};
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
    /*!
     * \brief calculate distance from simplex to point.
     * @param point coordinates of a point.
     */
    double distance(const std::vector<double> & point) const;
    /*!
     * \brief compute the normal to the simplex (only for dimMisMatch = 1).
     */
    EVector computeNormal() const;
  protected:
    /*!
     * \brief remove the column i from the matrix mat and return a new matrix.
     */
    EMatrix removeColumn(const EMatrix & mat, unsigned int colToRemove) const;
    /*!
     * integer of the dimension mismatch ($d_{points}-d$)
     */
    int dimMismatch;
};

template<int d>
void Simplex<d>::setPoints(std::vector< std::vector<double> > & 
    points_candidate){
  CPolytope::setPoints(points_candidate);
  if(points.size() != d+1){
    ErrorHandle eh("Simplex", "PointsConstructor", 
        "The number of points ("+std::to_string(points.size())+
        ") is not consistent with the dimension of the Simplex ("
        +std::to_string(dimension)+")." );
    throw eh;
  }
  this->dimMismatch = dimension - d;
  if((dimMismatch < 0){
    ErrorHandle eh("Simplex", "PointsConstructor", 
        "The dimension of the points ("+std::to_string((points[0]).size())+
        ") is not consistent with the dimension of the Simplex ("
        +std::to_string(dimension)+")." );
    throw eh;
  }
};

template<int d>
void Simplex<d>::setFaces(){
  this->faces.resize(d+1);
  int i, j, k;
  int dimmod2 = (d+1) % 2;
  for(i = 0; i < (d+1); i++){
    faces[i].resize(d);
    bool imod2 = (i+dimmod2) % 2;
    k = 0;
    for(j = 0; j < (d+1); j++){
      if(imod2){
        if(i != (d-j)){
          (faces[i])[k] = d-j;
          k++;
        }
      }
      else{
        if(i != j){
        (faces[i])[k] = j;
        k++;
        }
      }
    }
  }
}

template<int d>
Simplex<d>::Simplex(std::vector< std::vector<double> > & points_candidate){
    this->setPoints(points_candidate);
    if(d != 0){
      this->setFaces();
    }
};

template<int d>
EVector Simplex<d>::computeNormal() const{
  EVector normal(dimension);
  if(dimMisMatch != 1){
    ErrorHandle eh("Simplex". "computeNormal", "cannot compute the normal of " + 
        "a simplex that's nor embeded in a higher dimension (dimMismatch = "<< 
        this->dimMismatch << ").");
    throw eh;
    return normal;
  }else{
    int i, j;
    // construct tangent matrix
    EMatrix tanMat(dimension-1, dimension);
    for(i = 1; i < dimension; i++){
      for(j = 0; j < dimension; j++){
        tanMat(i-1, j) = points[i][j] - points[0][j];
      }
    }
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
}

template<int d>
EMatrix Simplex<d>::removeColumn(const EMatrix & mat, unsigned int colToRemove) 
  const{
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

template<int d>
double Simplex<d>::distance(const std::vector<double> & point) const{
  if(dimMismatch == 1){
    EVector normal = this->computeNormal();
    EVector pointv(dimension);
    for(i = 0; i < dimension; i++){
      pointv << point[i] - points[0][i];
    }
    return n.dot(pointv);
  }else{
    std::vector<double> distances;
    std::vector<std::vector<double> > temppoints(d);
    Simplex<d-1> hyperSim; 
    for(int i = 0; i < faces.size(); i++){
      for(int j = 0; j < (faces[i]).size(); j++){
        temppoints[j] = points[faces[i][j]];
      }
    }
  }
}

} //hfox

#endif //SIMPLEX_H
