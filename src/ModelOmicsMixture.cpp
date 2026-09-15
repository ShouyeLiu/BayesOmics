// SPDX-License-Identifier: GPL-3.0-or-later
#include "Data.hpp"
#include "OmicsMixture.hpp"
#include "Stat.hpp"
OmicsMixture::Draws::Draws(const Data &data) {
    offset.resize(data.numIncdSnps+1,0);
    for(unsigned j=0;j<data.numIncdSnps;++j) {
        auto genes=data.gwasSnpID2geneIDMap.find(data.snpEffectNames[j]);
        offset[j+1]=offset[j]+(genes==data.gwasSnpID2geneIDMap.end()?1:2*genes->second.size());
    }
    uniform.resize(offset.back());normal.resize(offset.back());
    for(Eigen::Index i=0;i<uniform.size();++i){uniform[i]=Stat::ranf();normal[i]=Stat::snorm();}
}
