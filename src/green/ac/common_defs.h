/*
* Copyright (c) 2023 University of Michigan
*
* Permission is hereby granted, free of charge, to any person obtaining a copy of this
* software and associated documentation files (the “Software”), to deal in the Software
* without restriction, including without limitation the rights to use, copy, modify,
* merge, publish, distribute, sublicense, and/or sell copies of the Software, and to
* permit persons to whom the Software is furnished to do so, subject to the following
* conditions:
*
* The above copyright notice and this permission notice shall be included in all copies or
* substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED,
* INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
* PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
* FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
* OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*/

#ifndef GREEN_AC_COMMON_DEFS_H
#define GREEN_AC_COMMON_DEFS_H

#include <stdexcept>
#include <cstdio>

namespace green::ac {
  enum AC_KIND {
    Nevanlinna
  };

  inline int compare_version_strings(const std::string& v1, const std::string& v2) {
    int major_V1 = 0, minor_V1 = 0, patch_V1 = 0;
    int major_V2 = 0, minor_V2 = 0, patch_V2 = 0;

    char suffix_V1[32] = "";
    char suffix_V2[32] = "";

    int parsed_1 = std::sscanf(v1.c_str(), "%d.%d.%d%30s", &major_V1, &minor_V1, &patch_V1, suffix_V1);
    int parsed_2 = std::sscanf(v2.c_str(), "%d.%d.%d%30s", &major_V2, &minor_V2, &patch_V2, suffix_V2);

    if (parsed_1 < 3) {
      throw std::runtime_error("First version string (v1) failed to parse: '" + v1 + "'. Expected format: major.minor.patch[suffix]");
    }
    if (parsed_2 < 3) {
      throw std::runtime_error("Second version string (v2) failed to parse: '" + v2 + "'. Expected format: major.minor.patch[suffix]");
    }

    if (major_V1 != major_V2) {
      return major_V1 > major_V2 ? 1 : -1;
    }
    if (minor_V1 != minor_V2) {
      return minor_V1 > minor_V2 ? 1 : -1;
    }
    if (patch_V1 != patch_V2) {
      return patch_V1 > patch_V2 ? 1 : -1;
    }

    return 0;
  }
}

#endif  // GREEN_AC_COMMON_DEFS_H
