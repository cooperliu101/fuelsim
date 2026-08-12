#ifndef FUELSIM_EXODUS_MESH_IO_HPP
#define FUELSIM_EXODUS_MESH_IO_HPP

#include "fuelsim/mesh.hpp"

#include <string>

namespace fuelsim {

class ExodusMeshIo final {
  public:
    static UnstructuredQuad4Mesh read_quad4(const std::string& path);
    static void write_quad4(const std::string& path, const UnstructuredQuad4Mesh& mesh);
    static UnstructuredHex8Mesh read_hex8(const std::string& path);
    static void write_hex8(const std::string& path, const UnstructuredHex8Mesh& mesh);
};

} // namespace fuelsim

#endif
