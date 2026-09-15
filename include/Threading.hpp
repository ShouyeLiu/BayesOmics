// SPDX-License-Identifier: GPL-3.0-or-later
#ifndef BAYESOMICS_THREADING_HPP
#define BAYESOMICS_THREADING_HPP
#include <string>
namespace Threading {
struct Selection {
    int affinity=1, available=1, requested=0, selected=1;
    bool automatic=true;
};
// Available means permitted CPUs / scheduler allocation, not transient idle CPU.
Selection select(const std::string &request);
}
#endif
