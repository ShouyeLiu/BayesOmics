// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <cstdint>
#include <string>

namespace MemoryBudget {
struct Status { uint64_t limit=0, resident=0, available=0; };
void configure(const std::string &megabytes, const std::string &genotypes, const std::string &chains);
Status status();
bool enabled();
bool streamGenotypes(uint64_t rows, uint64_t columns);
bool streamChains(uint64_t bytes);
void require(uint64_t bytes, const std::string &purpose);
void requireTextRows(const std::string &file,uint64_t perRow,const std::string &purpose);
uint64_t bytes(uint64_t rows,uint64_t columns,uint64_t element=8);
std::string describe();
}
