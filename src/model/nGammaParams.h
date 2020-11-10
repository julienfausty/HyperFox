#ifndef NGAMMAPARAMS_H
#define NGAMMAPARAMS_H

namespace hfox{

namespace nGamma{
/*!
 * \brief a simple struct for parameterizing nGamma models
 */

enum DriftType{
  /*!
   * \brief no drift terms
   */
  None
};//DriftType

enum TermAssembly{
  Implicit,
  Explicit
};//TermAssembly

struct nGammaParams{
  /*!
   * \brief the type of drift term to use
   */
  DriftType drift = None;
  /*!
   * \brief how to assemble the drift term
   */
  TermAssembly driftAssembly = Implicit;
  /*!
   * \brief how to assemble the linear velocity terms
   */
  TermAssembly linearVelocityAssembly = Implicit;
  /*!
   * \brief how to assemble the non-linear velocity terms
   */
  TermAssembly nonLinearVelocityAssembly = Implicit;
  /*!
   * \brief the speed of sound
   */
  double cs = 1.0;
};//nGammaParams

};//nGamma

};//hfox


#endif//NGAMMAPARAMS_H
