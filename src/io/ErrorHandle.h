#ifndef ERRORHANDLE_H
#define ERRORHANDLE_H

#include <exception>
#include <string>
#include <vector>

namespace hfox{

/*!
 * \brief class for checking and handling errors.
 *
 * Class that inherits from std::exception and should be used to handle most 
 * types of errors in HyperFox.
 */

class ErrorHandle: public std::exception{
  public:
    // Constructor
    /*!
     * \brief empty constructor.
     */
    ErrorHandle() : err_msg("something went wrong.\n"), userclass("SomeClass"), 
    userfunction("someFunction"){};
    /*!
     * \brief personalized constructor.
     */
    ErrorHandle(std::string msg_cand) : err_msg(msg_cand), 
      userclass("SomeClass"), userfunction("someFunction"){};
    /*!
     * \brief constructor where the user specifies all.
     */
    ErrorHandle(std::string userclass_cand, std::string userfunction_cand, 
      std::string msg_cand) : err_msg(msg_cand), userclass(userclass_cand), 
      userfunction(userfunction_cand){};
    /*!
     * \brief constructor where the user specifies only the class and function.
     */
    ErrorHandle(std::string userclass_cand, std::string userfunction_cand) : 
      err_msg("something went wrong.\n"), userclass(userclass_cand), 
      userfunction(userfunction_cand){};
    //Destructor
    virtual ~ErrorHandle(){};
    /*!
     * \brief overload what function from std::exception.
     */
    virtual const char * what() const noexcept;
    //Checkers
    /*!
     * \brief check if a points list is dimensionally consistent.
     */
    void checkPointList(const std::vector< std::vector<double> > & 
        points);
    /*!
     * \brief check if an list of index lists is in right range.
     */
    void checkIndexList(int max_index, const std::vector< std::vector<int> > & 
        faces);
  protected:
    std::string err_msg;
    std::string userclass;
    std::string userfunction;
};

} //hfox

#endif //ERRORHANDLE_H
