#ifndef FUELSIM_TIME_TABLE_HPP
#define FUELSIM_TIME_TABLE_HPP

#include <string>
#include <vector>

namespace fuelsim {

class PiecewiseLinearTimeTable final {
  public:
    PiecewiseLinearTimeTable(std::string name, std::vector<double> times, std::vector<double> values);

    const std::string& name() const noexcept;
    const std::vector<double>& times() const noexcept;
    const std::vector<double>& values() const noexcept;
    double value(double time) const;

  private:
    std::string _name;
    std::vector<double> _times;
    std::vector<double> _values;
};

} // namespace fuelsim

#endif
