#ifndef DENSEEIGEN_H
#define DENSEEIGEN_H

/*!
 * A file for putting some cool type defs and including the Eigen/Dense pack.
 */

#include <eigen3/Eigen/Dense>

namespace hfox{

typedef Eigen::MatrixXd EMatrix;
typedef Eigen::VectorXd EVector;
template <class T>
using EMap = Eigen::Map<T>;

}//hfox

#endif //DENSEEIGEN_H
