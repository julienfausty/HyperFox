#ifndef RUNGEKUTTA_H
#define RUNGEKUTTA_H

#include "TimeScheme.h"
#include "RKType.h"

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
     */
    RungeKutta(const ReferenceElement * re, RKType=CrankNicolson);
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
    void apply(EMatrix * stiffness, EVector * rhs);
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
     * \brief a method for getting the current stage
     */
    int getStage() const;
  protected:
    /*!
     * \brief the Butcher table of the method
     */
    EMatrix bTable;
    /*!
     * \brief boolean tracking if the Butcher table is currently set
     */
    bool bset = 0;
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

#endif//RUNEGKUTTA_H
