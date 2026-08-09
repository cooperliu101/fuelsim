#ifndef FUELSIM_EXODUS_FILE_HPP
#define FUELSIM_EXODUS_FILE_HPP

#include <string>

namespace fuelsim::exodus_detail {

void check_exodus(int status, const std::string& operation);

class ExodusFile final {
  public:
    explicit ExodusFile(int id);

    ExodusFile(const ExodusFile&) = delete;
    ExodusFile& operator=(const ExodusFile&) = delete;

    ~ExodusFile();

    int id() const noexcept;

    void close(const std::string& operation = "Could not close Exodus file");

  private:
    int _id;
};

} // namespace fuelsim::exodus_detail

#endif
