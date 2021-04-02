// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_TIMING_COMP_TIMES_H_
#define SRC_TIMING_COMP_TIMES_H_

#include <array>
#include <cassert>

// Helper class for accumulating times during tests.
template <std::size_t kTimeLength>
class Times {
 public:
  Times(void) { std::fill(times.begin(), times.end(), 0.0); }

  Times& operator+=(const Times& a_rhs) {
    for (std::size_t n = 0; n < kTimeLength; ++n) {
      times[n] += a_rhs[n];
    }
    return *this;
  }

  double& operator[](const std::size_t a_index) {
    assert(a_index < kTimeLength);
    return times[a_index];
  }

  const double& operator[](const std::size_t a_index) const {
    assert(a_index < kTimeLength);
    return times[a_index];
  }

  double* data(void) { return times.data(); }

  const double* data(void) const { return times.data(); }

 private:
  std::array<double, kTimeLength> times;
};

#endif  // SRC_TIMING_COMP_TIMES_H_
