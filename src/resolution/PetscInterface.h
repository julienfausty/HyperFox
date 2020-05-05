#ifndef PETSCINTERFACE_H
#define PETSCINTERFACE_H

#include <vector>
#include <iostream>
#include <string>
#include "petscksp.h"
#include "petscerror.h"
#include "LinAlgebraInterface.h"
#include "ErrorHandle.h"
#include "PetscOpts.h"

namespace hfox{

/*!
 * \brief Interface class the the PetsC library for linear algebra
 */

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
    PetscInterface(PetscOpts options);
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
     * @param diagSparsePattern a pointer to a vector of the number of non zero diagonal entries of the matrix per row
     * @param offSparsePattern a pointer to a vector of the number of non zero off-diagonal entries of the matrix per row
     *
     * diagSparsePattern and offSparsePattern must be passed together (ie both NULL or both valid).
     */
    void allocate(int ndofs, const std::vector<int> * diagSparsePattern = NULL, 
        const std::vector<int> * offSparsePattern = NULL);
    /*!
     * \brief add an element to the matrix
     *
     * @param i line number
     * @param j column number
     * @param val the value in the matrix
     */
    void addValMatrix(int i, int j, const double & val);
    /*!
     * \brief add multiple values to the matrix
     *
     * @param ijs vector of tuples with the indexes of the matrix vals (line number, column number) 
     * @param vals a vector of values
     */
    void addValsMatrix(std::vector<int> & is, std::vector<int> & js, const double * vals);
    /*!
     * \brief add an element to the right hand side
     *
     * @param i line number
     * @param val the value in the right hand side
     */
    void addValRHS(int i, const double & val);
    /*!
     * \brief add multiple values to the right hand side
     *
     * @param is vector of ints with the indexes of the rhs vals
     * @param vals a vector of values
     */
    void addValsRHS(std::vector< int > & is, const double * vals);

    /*!
     * \brief set an element of the matrix
     *
     * @param i line number
     * @param j column number
     * @param val the value in the matrix
     */
    void setValMatrix(int i, int j, const double & val);
    /*!
     * \brief set multiple values of the matrix
     *
     * @param is vector of indexes of the matrix vals (row number) 
     * @param js vector of indexes of the matrix vals (col number) 
     * @param vals a vector of values
     */
    void setValsMatrix(std::vector<int> & is, std::vector<int> & js, const double * vals);
    /*!
     * \brief set an element of the right hand side
     *
     * @param i line number
     * @param val the value in the right hand side
     */
    void setValRHS(int i, const double & val);
    /*!
     * \brief set multiple values of the right hand side
     *
     * @param is vector of ints with the indexes of the rhs vals
     * @param vals a vector of values
     */
    void setValsRHS(std::vector< int > & is, const double * vals);
    /*!
     * \brief zero out the set of rows
     *
     * @param is vector of the ints of the rows to zero out
     */
    virtual void zeroOutRows(std::vector<int> & is);
    /*!
     * \brief command to have the package solve the assembled system
     *
     * @param solution vector to insert the solution into
     */
    void solve(std::vector<double> * solution);
    /*!
     * \brief assemble the system
     */
    void assemble();
   /*!
     * \brief assemble an intermediate system
     */
    void assembleFlush();
    /*!
     * \brief clear the values set in the system but keep the structure if possible
     */
    void clearSystem();
    /*!
     * \brief destroy the structures of the system setting the system to an uninitialized state
     */
    void destroySystem();
    /*!
     * \brief get a const pointer to the options struct
     */
    const PetscOpts * getOptions() const{return &myOptions;};
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
     * \brief helper method for setting up the class
     */
    void setUp();
    /*!
     * \brief the options struct determining the behavior of the solvers
     */
    PetscOpts myOptions;
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
    /*!
     * \brief boolean that tracks whether Petsc was initialized elsewhere or not
     */
    bool masterPetscClass = 0;
    /*!
     * \brief buffer PETSC error code to not deal with the setup and teardown all the time
     */
    PetscErrorCode pErr;

};//PetscInterface

}//hfox

#endif//PETSCINTERFACE_H
