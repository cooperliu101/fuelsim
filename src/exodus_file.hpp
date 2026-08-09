#ifndef FUELSIM_EXODUS_FILE_HPP
#define FUELSIM_EXODUS_FILE_HPP

#include <exodusII.h>

#include <stdexcept>
#include <string>

namespace fuelsim::exodus_detail {

inline void check_exodus(int status, const std::string& operation) {
    if (status < 0)
        throw std::runtime_error(operation + ": " + ex_strerror(status));
}

class ExodusFile final {
  public:
    explicit ExodusFile(int id) : _id(id) {}

    ExodusFile(const ExodusFile&) = delete;
    ExodusFile& operator=(const ExodusFile&) = delete;

    ~ExodusFile() {
        if (_id >= 0)
            ex_close(_id);
    }

    int id() const noexcept {
        return _id;
    }

    void close(const std::string& operation = "Could not close Exodus file") {
        const int id = _id;
        _id = -1;
        check_exodus(ex_close(id), operation);
    }

  private:
    int _id;
};

} // namespace fuelsim::exodus_detail

#endif
