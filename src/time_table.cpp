#include "fuelsim/time_table.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace fuelsim {

PiecewiseLinearTimeTable::PiecewiseLinearTimeTable(std::string name,
                                                   std::vector<double> times,
                                                   std::vector<double> values)
    : _name(std::move(name)), _times(std::move(times)),
      _values(std::move(values)) {
    if (_name.empty())
        throw std::invalid_argument("Time-table name must not be empty");
    if (_times.size() < 2 || _times.size() != _values.size())
        throw std::invalid_argument(
            "Time table requires at least two time/value pairs: " + _name);
    for (std::size_t index = 0; index < _times.size(); ++index) {
        if (!std::isfinite(_times[index]) || _times[index] < 0.0 ||
            !std::isfinite(_values[index]))
            throw std::invalid_argument(
                "Time-table entries must be finite with nonnegative times: " +
                _name);
        if (index > 0 && !(_times[index] > _times[index - 1]))
            throw std::invalid_argument(
                "Time-table times must be strictly increasing: " + _name);
    }
}

const std::string& PiecewiseLinearTimeTable::name() const noexcept {
    return _name;
}

const std::vector<double>& PiecewiseLinearTimeTable::times() const noexcept {
    return _times;
}

const std::vector<double>& PiecewiseLinearTimeTable::values() const noexcept {
    return _values;
}

double PiecewiseLinearTimeTable::value(double time) const {
    if (!std::isfinite(time) || time < 0.0)
        throw std::invalid_argument(
            "Time-table evaluation time must be finite and nonnegative");
    if (time <= _times.front())
        return _values.front();
    if (time >= _times.back())
        return _values.back();
    const auto upper = std::upper_bound(_times.begin(), _times.end(), time);
    const std::size_t right = static_cast<std::size_t>(upper - _times.begin());
    const std::size_t left = right - 1;
    const double fraction =
        (time - _times[left]) / (_times[right] - _times[left]);
    return (1.0 - fraction) * _values[left] + fraction * _values[right];
}

} // namespace fuelsim
