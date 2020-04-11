#ifndef PETSCINTERFACE_H
#define PETSCINTERFACE_H

#include "petscksp.h"
#include "LinAlgebraInterface.h"
#include "PetscOptions.h"

/*!
 * \brief Interface class the the PetsC library for linear algebra
 */

namespace hfox{

class PetscInterface : public LinAlgebraInterface{
  public:
    /*!
     * \brief default constructor
     */
    PetscInterface();
    /*!
     * \brief default constructor for the interface with options
     *
     * @param options a struct with the user controlable options
     */
    PetscInterface(PetscOptions options);
    /*!
     * \brief destructor for the class
     */
    ~PetscInterface();
    /*!
     * \brief create the PetsC objects to use (Matrix, RHS, Solver and Preconditionner)
     */
    void initialize();
    /*!
     * \brief configure the initialized objects with the PetscOptions struct
     */
    void configure();
    /*!
     * \brief allocate both the matrix and rhs
     *
     * @param ndofs number of degrees of freedom (i.e. M = ndofs x ndofs, b = ndofs)
     */
    void allocate(int ndofs);
    /*!
     * \brief set an element of the matrix
     *
     * @param i line number
     * @param j column number
     * @param val the value in the matrix
     */
    void setValMatrix(int i, int j, double & val);
    /*!
     * \brief set multiple values of the matrix
     *
     * @param ijs vector of tuples with the indexes of the matrix vals (line number, column number) 
     * @param vals a vector of values
     */
    void setValsMatrix(std::vector< std::tuple<int, int> > & ijs, std::vector<double> & vals);
    /*!
     * \brief set an element of the right hand side
     *
     * @param i line number
     * @param val the value in the right hand side
     */
    void setValRHS(int i, double & val);
    /*!
     * \brief set multiple values of the right hand side
     *
     * @param is vector of ints with the indexes of the rhs vals
     * @param vals a vector of values
     */
    void setValsRHS(std::vector< int > & is, std::vector<double> & vals);
    /*!
     * \brief command to have the package solve the assembled system
     *
     * @param solution vector to insert the solution into
     */
    void solve(std::vector<double> * solution);
    /*!
     * \brief clear the values set in the system but keep the structure if possible
     */
    void clearSystem();
    /*!
     * \brief get a const pointer to the options struct
     */
    const PetscOptions * getOptions() const{return &myOptions;};
    /*!
     * \brief get a const pointer to the lhs matrix
     */
    const Mat * getMatrix() const{return &M;};
    /*!
     * \brief get a const pointer to the rhs vector
     */
    const Vec * getRHS() const{return &b;};
    /*!
     * \brief get a const pointer to the KSP solver instance
     */
    const KSP * getKSP() const{return &kspSolver;};
    /*!
     * \brief get a const pointer to the preconditionner instance
     */
    const PC * getPreConditionner() const{return &preCond;};
  protected:
    /*!
     * \brief the options struct determining the behavior of the solvers
     */
    PetscOptions myOptions;
    /*!
     * \brief the left hand side (lhs) matrix where it all happens
     */
    Mat M;
    /*!
     * \brief the right hand side (rhs) vector
     */
    Vec b;
    /*!
     * \brief the Krylov Subspace linear system solver
     */
    KSP kspSolver;
    /*!
     * \brief the preconditionner context
     */
    PC preCond;

};//PetscInterface

}//hfox

#endif//PETSCINTERFACE_H
