// SPDX-License-Identifier: GPL-3.0-or-later
#include "MemoryBudget.hpp"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <limits>
#include <stdexcept>
#include <sys/resource.h>
#include <unistd.h>

namespace MemoryBudget {
namespace {
uint64_t cap=0;
std::string genotypeMode="auto",chainMode="auto";
constexpr uint64_t MiB=1024*1024;
uint64_t resident() {std::ifstream f("/proc/self/statm");uint64_t total=0,rss=0;f>>total>>rss;return rss*uint64_t(sysconf(_SC_PAGESIZE));}
uint64_t number(const std::string &path) {std::ifstream f(path);std::string s;f>>s;if(s.empty()||s=="max")return 0;try{return std::stoull(s);}catch(...){return 0;}}
uint64_t detected() {
    uint64_t available=0;std::ifstream f("/proc/meminfo");std::string line;
    while(std::getline(f,line)){std::istringstream in(line);std::string key;uint64_t kb;in>>key>>kb;if(key=="MemAvailable:")available=kb*1024;}
    uint64_t result=available?available+resident():uint64_t(sysconf(_SC_PHYS_PAGES))*sysconf(_SC_PAGESIZE);
    // Walk all visible ancestors: a leaf may be unlimited under a capped job.
    std::ifstream cg("/proc/self/cgroup");
    while(std::getline(cg,line)) {
        auto a=line.find(':'),b=line.find(':',a+1);if(b==std::string::npos)continue;
        const std::string controllers=line.substr(a+1,b-a-1);bool v2=controllers.empty();
        if(!v2 && (","+controllers+",").find(",memory,")==std::string::npos)continue;
        std::string base=v2?"/sys/fs/cgroup":"/sys/fs/cgroup/memory", path=line.substr(b+1);
        if(path.find("..")!=std::string::npos)path="/";
        for(;;) {
            auto limit=number(base+path+(v2?"/memory.max":"/memory.limit_in_bytes"));
            auto used=number(base+path+(v2?"/memory.current":"/memory.usage_in_bytes"));
            if(limit && limit<(uint64_t(1)<<60))result=std::min(result,resident()+(limit>used?limit-used:0));
            if(path.empty()||path=="/")break;auto slash=path.find_last_of('/');path=slash==0?"/":path.substr(0,slash);
        }
    }
    rlimit lim{};if(getrlimit(RLIMIT_AS,&lim)==0 && lim.rlim_cur!=RLIM_INFINITY) {
        std::ifstream vm("/proc/self/statm");uint64_t pages=0;vm>>pages;uint64_t used=pages*sysconf(_SC_PAGESIZE);
        result=std::min(result,resident()+(lim.rlim_cur>used?lim.rlim_cur-used:0));
    }
    // Slurm's declared allocation is an additional ceiling even without cgroups.
    const char *node=getenv("SLURM_MEM_PER_NODE"),*cpu=getenv("SLURM_MEM_PER_CPU"),*count=getenv("SLURM_CPUS_PER_TASK");
    try {uint64_t allocation=0;if(node)allocation=std::stoull(node)*MiB;else if(cpu && count)allocation=std::stoull(cpu)*std::stoull(count)*MiB;
        if(allocation)result=std::min(result,allocation); // Slurm --mem=0 means all node memory
    }catch(...){}
    return result;
}
}
uint64_t bytes(uint64_t rows,uint64_t columns,uint64_t element) {
    if(columns && rows>std::numeric_limits<uint64_t>::max()/columns)throw std::overflow_error("Memory dimension overflow");
    uint64_t n=rows*columns;if(element && n>std::numeric_limits<uint64_t>::max()/element)throw std::overflow_error("Memory size overflow");return n*element;
}
void configure(const std::string &mb,const std::string &g,const std::string &m) {
    for(const auto &mode:{g,m})if(mode!="auto" && mode!="dense" && mode!="stream")throw std::invalid_argument("Storage must be auto, dense or stream");
    genotypeMode=g;chainMode=m;uint64_t detectedLimit=detected();
    if(mb=="auto")cap=detectedLimit/10*8;
    else {size_t end=0;double value=0;try{value=std::stod(mb,&end);}catch(...){throw std::invalid_argument("--memory requires auto or a positive MiB value");}
        if(end!=mb.size()||!(value>0)||value>double(std::numeric_limits<uint64_t>::max()/MiB))throw std::invalid_argument("Invalid --memory MiB value");
        cap=std::min(uint64_t(value*MiB),detectedLimit/10*9);
    }
    if(!cap)throw std::runtime_error("No usable memory budget detected");
}
bool enabled(){return cap!=0;}
Status status(){
    Status s;s.resident=resident();
    // The requested budget is a ceiling, not a reservation. Other processes or
    // a tighter cgroup/rlimit can reduce headroom after configure() returned.
    // Refresh before allocation/storage decisions rather than trusting the
    // startup snapshot for the lifetime of a long MCMC run.
    s.limit=cap?std::min(cap,detected()):0;
    s.available=s.limit>s.resident?s.limit-s.resident:0;
    return s;
}
bool streamGenotypes(uint64_t n,uint64_t p){
    if(!enabled())return false;
    auto needed=bytes(n,p);if(genotypeMode=="stream")return true;
    // Leave half of remaining memory for model state, molecular data and output.
    if(genotypeMode=="dense"){require(needed,"dense genotype matrix");return false;}
    return needed>status().available/2;
}
bool streamChains(uint64_t n){if(!enabled())return false;if(chainMode=="stream")return true;if(chainMode=="dense"){require(n,"retained MCMC matrices");return false;}return n>status().available/3;}
void require(uint64_t n,const std::string &purpose){if(!enabled())return;auto s=status();if(n>s.available/10*9)throw std::runtime_error("Memory budget insufficient for "+purpose+": needs "+std::to_string((n+MiB-1)/MiB)+" MiB additional; "+std::to_string(s.available/MiB)+" MiB remains. Increase --memory/allocation or reduce the problem/block size.");}
void requireTextRows(const std::string &file,uint64_t perRow,const std::string &purpose){
    if(!enabled())return;std::ifstream in(file);if(!in)return;uint64_t rows=0;std::string line;while(std::getline(in,line))if(!line.empty())++rows;require(bytes(rows,perRow,1),purpose);
}
std::string describe(){auto s=status();return "Memory budget: "+std::to_string(s.limit/MiB)+" MiB; current RSS "+std::to_string(s.resident/MiB)+" MiB (host/cgroup/rlimit/allocation ceilings; not an OS reservation)";}
}
