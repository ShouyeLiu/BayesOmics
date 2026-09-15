// SPDX-License-Identifier: GPL-3.0-or-later
#include "GenotypeView.hpp"
#include "MemoryBudget.hpp"
#include <fcntl.h>
#include <unistd.h>
#include <sys/stat.h>
#include <stdexcept>
#include <cerrno>

namespace {void readAt(int fd,void *ptr,size_t size,uint64_t offset) {
    auto out=static_cast<unsigned char*>(ptr);
    while(size){auto n=pread(fd,out,size,offset);if(n<0 && errno==EINTR)continue;if(n<=0)throw std::runtime_error("Truncated/unreadable SNP-major BED file");out+=n;size-=n;offset+=n;}
}}
BedColumns::BedColumns(const std::string &path,unsigned rawIndividuals,unsigned rawSnps,std::vector<unsigned> inds,std::vector<unsigned> markers):
stride((uint64_t(rawIndividuals)+3)/4),individuals(std::move(inds)),snps(std::move(markers)) {
    MemoryBudget::require(MemoryBudget::bytes(rows(),4)+MemoryBudget::bytes(cols(),3)+stride,"one BED column and scaling metadata");
    fd=open(path.c_str(),O_RDONLY);if(fd<0)throw std::runtime_error("Cannot open BED: "+path);
    try {
        unsigned char h[3];readAt(fd,h,3,0);if(h[0]!=0x6c||h[1]!=0x1b||h[2]!=1)throw std::runtime_error("BED must be SNP-major");
        struct stat info{};if(fstat(fd,&info)!=0 || uint64_t(info.st_size)<3+MemoryBudget::bytes(rawSnps,stride,1))throw std::runtime_error("BED is truncated relative to BIM/FAM dimensions");
        for(auto i:individuals)if(i>=rawIndividuals)throw std::out_of_range("BED individual mapping");
        for(auto j:snps)if(j>=rawSnps)throw std::out_of_range("BED SNP mapping");
        packed.resize(stride);decoded.resize(rows());means.setZero(cols());scales.setOnes(cols());
    }catch(...){close(fd);fd=-1;throw;}
}
BedColumns::~BedColumns(){if(fd>=0)close(fd);}
Eigen::VectorXd BedColumns::raw(Eigen::Index j)const {
    if(j<0||j>=cols())throw std::out_of_range("BED column");
    readAt(fd,packed.data(),packed.size(),3+stride*snps[j]);Eigen::VectorXd result(rows());const int decode[4]={2,-9,1,0};
    for(Eigen::Index i=0;i<rows();++i){auto raw=individuals[i];result[i]=decode[(packed[raw/4]>>((raw%4)*2))&3];}return result;
}
Eigen::Ref<const Eigen::VectorXd> BedColumns::col(Eigen::Index j)const {
    if(cached!=j){decoded=raw(j);for(Eigen::Index i=0;i<rows();++i)if(decoded[i]==-9)decoded[i]=means[j];
        decoded.array()-=means[j];decoded.array()/=scales[j];if(weights.size())decoded.array()*=weights.array();cached=j;}
    return decoded;
}
void BedColumns::retain(const std::vector<int> &columns){
    std::vector<unsigned> kept;Eigen::VectorXd m(columns.size()),s(columns.size());
    for(size_t i=0;i<columns.size();++i){kept.push_back(snps.at(columns[i]));m[i]=means[columns[i]];s[i]=scales[columns[i]];}
    snps=std::move(kept);means=std::move(m);scales=std::move(s);invalidate();
}
const Eigen::MatrixXd &GenotypeView::dense()const {if(bed)throw std::runtime_error("This operation requires dense genotypes; use --genotype-storage dense with sufficient --memory");return *matrix;}
