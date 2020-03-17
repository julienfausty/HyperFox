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
 */

class LinAlgebraInterface{

  public:
    /*!
     * \brief initialize the objects that the package needs
     */
    virtual void initialize()=0;
    /*!
     * \brief configure the initialized objects with the user defined parameters
     *
     * @param configuration a map with key value pairs set to strings that should completely configure the
     * behavior of the package.
     */
    virtual void configure(std::map<std::string, std::string> configuration)=0;
    /*!
     * \brief allocate for the size of the linear system
     *
     * @param ndofs number of degrees of freedom (i.e. A = ndofs x ndofs, b = ndofs)
     */
    virtual void allocate(int ndofs)=0;
    /*!
     * \brief set an element of the matrix
     *
     * @param i line number
     * @param j column number
     * @param val the value in the matrix
     */
    virtual void setValMatrix(int i, int j, double val)=0;
    /*!
     * \brief set multiple values of the matrix
     *
     * @param ijs vector of tuples with the indexes of the matrix vals (line number, column number) 
     * @param vals a vector of values
     */
    virtual void setValsMatrix(std::vector< std::tuple<int, int> > ijs, std::vector<double> vals)=0;
    /*!
     * \brief set an element of the right hand side
     *
     * @param i line number
     * @param val the value in the right hand side
     */
    virtual void setValRHS(int i, double val)=0;
    /*!
     * \brief set multiple values of the right hand side
     *
     * @param is vector of ints with the indexes of the rhs vals
     * @param vals a vector of values
     */
    virtual void setValsRHS(std::vector< int > is, std::vector<double> vals)=0;
    /*!
     * \brief command to have the package solve the assembled system and return the result as a vector
     */
    virtual std::vector<double> solve()=0;


  protected:
    /*!
     * \brief number of degrees of freedom and thus the size of the linear system
     */
    int nDOFs;
};//LinAlgebraInterface

}//hfox

#endif//LINALGEBRAINTERFACE_H
