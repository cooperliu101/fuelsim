#include "exodus_file.hpp"

#include <exodusII.h>

#include <stdexcept>

namespace fuelsim::exodus_detail {

void check_exodus(int status, const std::string& operation) {
    if (status < 0)
        throw std::runtime_error(operation + ": " + ex_strerror(status));
}

ExodusFile::ExodusFile(int id) : _id(id) {}

ExodusFile::~ExodusFile() {
    if (_id >= 0)
        ex_close(_id);
}

int ExodusFile::id() const noexcept {
    return _id;
}

void ExodusFile::close(const std::string& operation) {
    const int id = _id;
    _id = -1;
    check_exodus(ex_close(id), operation);
}

} // namespace fuelsim::exodus_detail
