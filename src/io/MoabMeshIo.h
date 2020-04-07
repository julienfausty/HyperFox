#ifndef MOABMESHIO_H
#define MOABMESHIO_H

#include <iostream>
#include <moab/Core.hpp>
#include "Io.h"
#include "Mesh.h"

namespace hfox{

/*!
 * \brief The Input/Output class that plugs in to moab's capabilities (load_file, write_file).
 *
 * A class to be used to import meshes directly through moab's interface. The file formats that are 
 * supported are those that moab supports with the conventions that it uses.
 */

class MoabMeshIo : public Io{

  public:
    /*!
     * \brief an empty constructor
     */
    MoabMeshIo();
    /*!
     * \brief a constructor passing a mesh
     *
     * @param mesh a pointer to a Mesh object (the mesh should already have a reference element)
     */
    MoabMeshIo(Mesh * mesh);
    /*!
     * \brief method to load a mesh in a file into the mesh object
     *
     * @param filename name of file to read from
     */
    void load(std::string filename);
    /*!
     * \brief method to write a mesh into a file
     *
     * @param filename name of file to write to
     */
    void write(std::string filename);
  protected:
    /*!
     * \brief a method for entering a moab meshset into a Mesh
     *
     * @param mbIface the interface to the moab instance
     * @param meshset the mesh set to input into mesh
     */
    void setMeshSet(moab::Interface * mbIFace, moab::EntityHandle & meshset);
    /*! 
     * \brief method for determining MOAB element type
     */
    moab::EntityType determineMOABType() const;

};//MoabMeshIo

};//hfox

#endif
