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
     * \brief a method to set the mesh pointer (the mesh should already have a reference element)
     *
     * @param mesh a pointer to a Mesh object
     */
    void setMesh(Mesh * mesh);
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
     * \brief a pointer to the associated mesh
     */
    Mesh * myMesh;
    /*!
     * \brief a pointer to a moab interface
     */
    moab::Interface * mbIFace;

};//MoabMeshIo

};//hfox

#endif
