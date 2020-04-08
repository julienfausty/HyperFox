#ifndef HDF5IO_H
#define HDF5IO_H

#include <string>
#include "Io.h"
#include "ErrorHandle.h"
#include "hdf5.h"

/*!
 * \brief The Input/Output class that manages the native HDF5 format
 *
 * The HDF5 format described here is simple. The hiearchy is as follows:
 *
 * - Mesh
 *   - Nodes: vals
 *   - Cells: vals
 * - FieldData
 *   - NameOfField: fieldVals (attributes: (ftype : field type integer 0->nodes, 2->faces, 3->cells))
 *
 * Default behavior is to write everything held in the class
 */

namespace hfox{

class HDF5Io : public Io{

  public:
    /*!
     * \brief empty constructor
     */
    HDF5Io();
    /*!
     * \brief mesh constructor
     */
    HDF5Io(Mesh * mesh);
    /*!
     * \brief loading HDF5 files
     */
    void load(std::string filename);
    /*!
     * \brief writing HDF5 files
     */
    void write(std::string filename);
  protected:
    /*!
     * \brief for loading a mesh from the file
     *
     * @param meshGrp HDF5 handle to opened mesh group
     */
    void loadMesh(hid_t meshGrp);
    /*!
     * \brief for loading fields from the file
     *
     * @param fieldGrp HDF5 handle to opened field group
     */
    void loadFields(hid_t fieldGrp);
    /*!
     * \brief for writing a mesh to the file
     *
     * @param meshGrp HDF5 handle to created mesh group
     */
    void writeMesh(hid_t meshGrp);
    /*!
     * \brief for writing fields to the file
     *
     * @param fieldGrp HDF5 handle to created field group
     */
    void writeFields(hid_t fieldGrp);

};//HDF5Io

}//hfox

#endif
