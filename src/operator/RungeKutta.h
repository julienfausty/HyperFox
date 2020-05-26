#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include "TimeScheme.h"
#include "RKType.h"
#include "Field.h"

namespace hfox{

/*!
 * \brief A time scheme object implementing Runge-Kutta type time integration
 *
 * For now: can take ERK and DIRK types of schemes but does not implement fully implicit type RK
 */

class RungeKutta : public TimeScheme{

  public:
    /*!
     * \brief Constructor
     *
     * @param re the reference element
     * @param type the Runge-Kutta method type
     * @param field a list of auxiliary field names
     */
    RungeKutta(const ReferenceElement * re, RKType type=CrankNicolson, std::vector<std::string> fields=std::vector<std::string>());
    /*!
     * \brief sets the field map for the class
     *
     * @param fm the field map
     */
    void setFieldMap(const std::map<std::string, std::vector<double> > * fm);
    /*!
     * \brief the method that applies the time scheme to the local stiffness matrix and rhs inplace
     *
     * @param stiffness pointer to the stiffness matrix that will become the total matrix
     * @param rhs pointer to the right hand side
     */
    void apply(EMatrix * stiffness, EVector * s);
    /*!
     * \brief sets the butcher table for the class
     *
     * @param butchTable a matrix defining a Butcher table : $\begin{array}{c|c} c & A\\ & b\end{array}$
     */
    void setButcherTable(EMatrix butchTable);
    /*!
     * \brief sets the butcher table for the class using a preset
     *
     * @param type an RKType associated to a specific Butcher table
     */
    void setButcherTable(RKType type);
    /*!
     * \brief set the auxiliary field names
     *
     * @param fields a list of auxiliary field names
     */
    void setAuxiliaryFields(std::vector<std::string> fields);
    /*!
     * \brief compute the current stage from the solution of the system
     */
    void computeStage(std::map<std::string, Field*> * fm);
    /*!
     * \brief compute the solution from the stages of the system
     */
    void computeSolution(std::map<std::string, Field*> * fm);
    /*!
     * \brief a method for getting the current stage
     */
    int getStage() const{return stageCounter;};
    /*!
     * \brief get number of stages
     */
    int getNumStages() const;
  protected:
    /*!
     * \brief the Butcher table of the method
     */
    EMatrix bTable;
    /*!
     * \brief names of the auxiliary fields that need to be considered
     */
    std::vector<std::string> auxiliaryFields;
    /*!
     * \brief a stage counter for the multi-stage process
     */
    int stageCounter;
    /*!
     * \brief a database of Butcher tables preset for different methods
     */
    static std::map<RKType, EMatrix> butcherDB;
    /*!
     * \brief set up the butcher DB
     */
    static void setUpDB();

};//RungeKutta

}//hfox

#endif//RUNGEKUTTA_H
