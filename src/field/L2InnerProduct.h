#ifndef L2INNERPRODUCT_H
#define L2INNERPRODUCT_H

#include <vector>
#include <functional>
#include "DenseEigen.h"
#include "Field.h"
#include "Integrator.h"

namespace hfox{

/*!
 * \brief This integrator computes an L2 inner product either between two fields or functions or a combination of both.
 */

class L2InnerProduct : public Integrator{

  public:

    using Integrator::Integrator;
    
    /*!
     * \brief useful typedef for defining the functions that come into play here
     */
    typedef std::function<void(std::vector<double>&, std::vector<double>*)> FieldFunc;

    /*!
     * \brief set the left and right sides of the product in FieldField mode 
     *
     * @param lfield the left hand field of the product
     * @param rfield the right hand field of the product
     */
    void setProduct(Field* lfield, Field* rfield);

     /*!
     * \brief set the left and right sides of the product in FieldFunction mode 
     *
     * @param lfield the left hand field of the product
     * @param rfield the right hand field of the product
     */
    void setProduct(FieldFunc lfield, Field* rfield); 

     /*!
     * \brief set the left and right sides of the product in FieldFunction mode 
     *
     * @param lfield the left hand field of the product
     * @param rfield the right hand field of the product
     */
    void setProduct(Field* lfield, FieldFunc rfield); 

     /*!
     * \brief set the left and right sides of the product in FieldFunction mode 
     *
     * @param lfield the left hand field of the product
     * @param rfield the right hand field of the product
     */
    void setProduct(FieldFunc lfield, FieldFunc rfield); 

  protected:

    /*!
     * \brief evalute the integrand at the integration points of the entity indexed by iEnt
     *
     * @param iEnt index of the entity we need the integrand on
     * @param ipVals an array to fill with the values
     */
    void evaluateIntegrand(int iEnt, std::vector<double> * ipVals);

    /*!
     * \brief helper function to check individual properties of a Field
     */
    void checkField(Field * pField);

    /*!
     * \brief an enum for the different modes the product can be in
     */
    enum ProductMode{FieldField, FieldFunction, FunctionFunction};

    /*!
     * \brief the mode of the product
     */
    ProductMode pmode;

    /*!
     * \brief the vector dimension of the fields
     */
    int vectorDim = 0;

    /*!
     * \brief a list of Fields (max 2)
     */
    std::vector<Field*> fields;

    /*!
     * \brief a list of functions that act like fields
     */
    std::vector<FieldFunc> fieldFuncs;

};//L2InnerProduct

};//hfox

#endif//L2INNERPRODUCT_H
