// Copyright (C) 2020 Robert Chiodi <robert.chiodi@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#ifndef SRC_TIMING_COMP_FILES_H_
#define SRC_TIMING_COMP_FILES_H_

#include <iostream>
#include <string>

// Simple RAII wrapper for files written to during testing.
struct Files {
  Files(void) = delete;

  Files(std::string irl_name, std::string r3d_name, std::string voftools_name) {
    irl = fopen(irl_name.c_str(), "w");
    r3d = fopen(r3d_name.c_str(), "w");
    voftools = fopen(voftools_name.c_str(), "w");
  }

  void writeToFiles(const std::string& a_string) {
    fprintf(irl, "%s", a_string.c_str());
    fprintf(r3d, "%s", a_string.c_str());
    fprintf(voftools, "%s", a_string.c_str());
  }

  ~Files(void) {
    fclose(irl);
    fclose(r3d);
    fclose(voftools);
  }

  FILE* irl;
  FILE* r3d;
  FILE* voftools;
};

#endif  // SRC_TIMING_COMP_FILES_H_
