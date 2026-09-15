// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include "Stat.hpp"
#include "MemoryBudget.hpp"
#include <map>
#include <numeric>
#include <string>
#include <vector>
#include <stdexcept>

// Scheduling and random-number storage only. Conditional distributions, effect
// updates and sufficient statistics stay in each independent EIEO method.
namespace IndependentBlocks {
inline std::vector<std::vector<unsigned>> groups(
    const std::map<int,std::vector<int>> &blocks,
    const std::map<int,std::string> &ids,
    const std::map<std::string,int> &geneIndex,
    const std::map<std::string,std::vector<std::string>> &snpGenes) {
    std::vector<unsigned> parent(blocks.size());
    std::iota(parent.begin(),parent.end(),0);
    auto root=[&](unsigned b){while(parent[b]!=b){parent[b]=parent[parent[b]];b=parent[b];}return b;};
    std::vector<int> geneOwner(geneIndex.size(),-1),snpOwner(ids.size(),-1);
    for(unsigned b=0;b<blocks.size();++b) for(int j:blocks.at(b)) {
        if(j<0 || unsigned(j)>=ids.size() || snpOwner[j]!=-1)
            throw std::invalid_argument("LD blocks must partition SNPs without duplicate membership");
        snpOwner[j]=b;
        auto mapped=snpGenes.find(ids.at(j));
        if(mapped==snpGenes.end())continue;
        for(const auto &gene:mapped->second) {
            const unsigned g=geneIndex.at(gene);
            if(geneOwner.at(g)<0)geneOwner[g]=b;
            else {auto a=root(b),c=root(geneOwner[g]);if(a!=c)parent[std::max(a,c)]=std::min(a,c);}
        }
    }
    for(int owner:snpOwner)if(owner<0)throw std::invalid_argument("SNP missing from LD block partition");
    std::map<unsigned,std::vector<unsigned>> merged;
    for(unsigned b=0;b<blocks.size();++b)merged[root(b)].push_back(b);
    std::vector<std::vector<unsigned>> result;
    for(auto &group:merged)result.push_back(std::move(group.second));
    return result;
}

struct Draws {
    std::vector<std::size_t> offset;
    std::vector<double> uniform,normal;
    Draws(const std::map<int,std::string> &ids,
          const std::map<std::string,std::vector<std::string>> &snpGenes) {
        MemoryBudget::require(MemoryBudget::bytes(ids.size()+1,1,sizeof(std::size_t)),"SNP random-offset workspace");
        offset.resize(ids.size()+1,0);
        for(unsigned j=0;j<ids.size();++j) {
            auto mapped=snpGenes.find(ids.at(j));
            offset[j+1]=offset[j]+(mapped==snpGenes.end()?1:2*mapped->second.size());
        }
        MemoryBudget::require(MemoryBudget::bytes(offset.back(),2),"Independent-block random variates");
        uniform.resize(offset.back());normal.resize(offset.back());
        // Fixed SNP/trait order; no global RNG or Armadillo RNG inside workers.
        // As in GCTB, reserve normals even for components which prove inactive.
        for(std::size_t i=0;i<uniform.size();++i){uniform[i]=Stat::ranf();normal[i]=Stat::snorm();}
    }
};
}
