// SPDX-License-Identifier: GPL-3.0-or-later
#include "Threading.hpp"
#include <algorithm>
#include <cerrno>
#include <climits>
#include <cstdlib>
#include <stdexcept>
#include <thread>
#ifdef __linux__
#include <sched.h>
#endif

namespace {
int positive(const std::string &text,bool firstLevel=false) {
    if(text.empty() || text[0]<'0' || text[0]>'9')return 0;
    errno=0;char *end=nullptr;long value=std::strtol(text.c_str(),&end,10);
    if(errno || value<1 || value>INT_MAX || (*end && !(firstLevel && *end==',')))return 0;
    return static_cast<int>(value);
}
int environmentLimit(const char *name,bool firstLevel=false) {
    const char *value=std::getenv(name);return value?positive(value,firstLevel):0;
}
}
Threading::Selection Threading::select(const std::string &request) {
    Selection out;out.affinity=std::max(1u,std::thread::hardware_concurrency());
#ifdef __linux__
    cpu_set_t allowed;CPU_ZERO(&allowed);
    if(sched_getaffinity(0,sizeof(allowed),&allowed)==0)out.affinity=std::max(1,CPU_COUNT(&allowed));
#endif
    out.available=out.affinity;
    for(const char *name:{"SLURM_CPUS_PER_TASK","SLURM_CPUS_ON_NODE","OMP_THREAD_LIMIT"}) {
        int limit=environmentLimit(name);if(limit)out.available=std::min(out.available,limit);
    }
    out.automatic=request=="auto";
    if(out.automatic) {
        out.selected=out.available;
        // An explicit OpenMP environment preference remains meaningful in auto mode.
        int preference=environmentLimit("OMP_NUM_THREADS",true);
        if(preference)out.selected=std::min(out.selected,preference);
    } else {
        out.requested=positive(request);
        if(!out.requested)throw std::invalid_argument("--thread must be a positive integer or auto");
        out.selected=std::min(out.requested,out.available);
    }
    return out;
}
