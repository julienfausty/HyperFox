#ifndef IO_H
#define IO_H

#include <string>
#include <boost/filesystem.hpp>

namespace hfox{

/*!
 * \brief The Input/Output class.
 *
 * A virtual class that sets up the interface for any input or output operations one may perform.
 * Every instance of the class should have a unique format (e.g. vtu, hdf5, txt, ...).
 */

class Io{

  public:
    /*!
     * \brief Every Io instance should have an input method that takes a file name and loads it.
     *
     * @param filename the name of the file on wishes to load.
     */
    virtual void load(std::string filename)=0;
    /*!
     * \brief Every Io instance should have an ouput method that takes a file name and writes it.
     *
     * @param filename the name of the file on wishes to write to.
     */
    virtual void write(std::string filename)=0;
  protected:
    /*!
     * \brief A hgelper function used to get the extension of a file name.
     *
     * @param filename the name of the file on wishes to write t.
     */
    std::string getExtension(std::string filename);

}; //Io

} //hfox

#endif
