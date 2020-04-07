#ifndef IO_H
#define IO_H

#include <string>
#include <boost/filesystem.hpp>
#include <map>
#include "ErrorHandle.h"
#include "Mesh.h"
#include "Field.h"

namespace hfox{

/*!
 * \brief The Input/Output interface.
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
    /*!
     * \brief a method to set the mesh pointer (the mesh should already have a reference element)
     *
     * @param mesh a pointer to a Mesh object
     */
    virtual void setMesh(Mesh * mesh);
    /*!
     * \brief a method to set a field to be written
     *
     * @param name name one would like to give to the field
     * @param field pointer to the field object
     */
    virtual void setField(std::string name, Field * field);
  protected:
    /*!
     * \brief A hgelper function used to get the extension of a file name.
     *
     * @param filename the name of the file on wishes to write t.
     */
    static std::string getExtension(std::string filename);
    /*!
     * \brief pointer to mesh object (when applicable)
     */
    Mesh * myMesh;
    /*!
     * \brief map from names of field to pointers to those fields (when applicable)
     */
    std::map<std::string, Field*> fieldMap;


}; //Io

} //hfox

#endif
