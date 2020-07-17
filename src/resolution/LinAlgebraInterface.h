#ifndef LINALGEBRAINTERFACE_H
#define LINALGEBRAINTERFACE_H

#include <map>
#include <string>
#include <tuple>
#include <vector>

namespace hfox{

/*!
 * \brief The interface for plugging linear algebra libraries.
 *
 * A virtual class that sets up the interface for any linear algebra package one would
 * like to use for solving the sparse matrix problem.
 *
 * AX = b
 *
 * Warning: the convention here will be a column major data ordering such that when using 
 * Eigen one must transpose the data before inputting
 */

class LinAlgebraInterface{

  public:
    /*!
     * \brief initialize the objects that the package needs
     */
    virtual void initialize()=0;
    /*!
     * \brief configure the initialized objects with the user defined parameters (these parameters should be passed in during construction of the object)
     */
    virtual void configure()=0;
    /*!
     * \brief allocate for the size of the linear system
     *
     * @param local ndofs number of degrees of freedom (i.e. A = ndofs x ndofs, b = ndofs) locally on partition
     * @param diagSparsePattern a pointer to a vector of the number of non zero diagonal entries of the matrix per row
     * @param offSparsePattern a pointer to a vector of the number of non zero off-diagonal entries of the matrix per row
     *
     * diagSparsePattern and offSparsePattern must be passed together (ie both NULL or both valid).
     */
    virtual void allocate(int ndofs, 
        const std::vector<int> * diagSparsePattern = NULL, 
        const std::vector<int> * offSparsePattern = NULL)=0;
    /*!
     * \brief add an element to the matrix
     *
     * @param i line number
     * @param j column number
     * @param val the value in the matrix
     */
    virtual void addValMatrix(int i, int j, const double & val)=0;
    /*!
     * \brief add multiple values to the matrix
     *
     * @param ijs vector of tuples with the indexes of the matrix vals (line number, column number) 
     * @param vals a vector of values
     */
    virtual void addValsMatrix(std::vector<int> & is, std::vector<int> & js, const double * vals)=0;
    /*!
     * \brief add an element to the right hand side
     *
     * @param i line number
     * @param val the value in the right hand side
     */
    virtual void addValRHS(int i, const double & val)=0;
    /*!
     * \brief add multiple values to the right hand side
     *
     * @param is vector of ints with the indexes of the rhs vals
     * @param vals a vector of values
     */
    virtual void addValsRHS(std::vector< int > & is, const double * vals)=0;
    /*!
     * \brief set an element of the matrix
     *
     * @param i line number
     * @param j column number
     * @param val the value in the matrix
     */
    virtual void setValMatrix(int i, int j, const double & val)=0;
    /*!
     * \brief set multiple values of the matrix
     *
     * @param ijs vector of tuples with the indexes of the matrix vals (line number, column number) 
     * @param vals a vector of values
     */
    virtual void setValsMatrix(std::vector<int> & is, std::vector<int> & js, const double * vals)=0;
    /*!
     * \brief set an element of the right hand side
     *
     * @param i line number
     * @param val the value in the right hand side
     */
    virtual void setValRHS(int i, const double & val)=0;
    /*!
     * \brief set multiple values of the right hand side
     *
     * @param is vector of ints with the indexes of the rhs vals
     * @param vals a vector of values
     */
    virtual void setValsRHS(std::vector< int > & is, const double * vals)=0;
    /*!
     * \brief zero out the set of rows
     *
     * @param is vector of the ints of the rows to zero out
     */
    virtual void zeroOutRows(std::vector<int> & is) = 0;
    /*!
     * \brief assemble the system
     */
    virtual void assemble()=0;
    /*!
     * \brief assemble an intermediate system
     */
    virtual void assembleFlush()=0;
    /*!
     * \brief command to have the package solve the assembled system
     *
     * @param solution vector to insert the solution into
     */
    virtual void solve(std::vector<double> * solution)=0;
    /*!
     * \brief get the solution ownership range of the partition
     * @param range vector of indexes to fill with the range
     */
    virtual void getSolutionOwnership(std::vector<int> * range)=0;
    /*!
     * \brief clear the values set in the system but keep the structure if possible
     */
    virtual void clearSystem()=0;
    /*!
     * \brief destroy the structures of the system setting the system to an uninitialized state
     */
    virtual void destroySystem()=0;
    /*!
     * \brief get the number of degrees of freedom in the system
     */
    virtual int getNumDofs() const{return nDOFs;};
  protected:
    /*!
     * \brief number of degrees of freedom and thus the size of the linear system
     */
    int nDOFs;
    /*!
     * \brief boolean value tracking if the object has been initialized
     */
    bool initialized = 0;
    /*!
     * \brief boolean value tracking if the object has been configured
     */
    bool configured = 0;
    /*!
     * \brief boolean value tracking if the object has been allocated
     */
    bool allocated = 0;
    /*!
     * \brief boolean value tracking if the object has been assembled (at least once)
     */
    bool assembled = 0;
};//LinAlgebraInterface

}//hfox

#endif//LINALGEBRAINTERFACE_H
