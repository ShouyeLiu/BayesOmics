// SPDX-License-Identifier: GPL-3.0-or-later
#pragma once
#include <Eigen/Core>
#include <memory>
#include <string>
#include <vector>

// BED stays packed on disk. One decoded column is cached; no N x P scratch file.
// col() is valid until the next different column on the SAME store. Callers must
// acquire columns before launching observation workers, not from those workers.
class BedColumns {
    int fd=-1;
    uint64_t stride=0;
    std::vector<unsigned> individuals, snps;
    mutable std::vector<unsigned char> packed;
    mutable Eigen::VectorXd decoded;
    mutable Eigen::Index cached=-1;
public:
    Eigen::VectorXd means, scales, weights;
    BedColumns(const std::string &path,unsigned rawIndividuals,unsigned rawSnps,
               std::vector<unsigned> selectedIndividuals,std::vector<unsigned> selectedSnps);
    ~BedColumns();
    BedColumns(const BedColumns&)=delete;
    Eigen::Index rows()const{return individuals.size();}
    Eigen::Index cols()const{return snps.size();}
    Eigen::VectorXd raw(Eigen::Index j)const;
    Eigen::Ref<const Eigen::VectorXd> col(Eigen::Index j)const;
    void invalidate(){cached=-1;}
    void retain(const std::vector<int> &columns);
};
class GenotypeView {
    const Eigen::MatrixXd *matrix=nullptr;
    std::shared_ptr<BedColumns> bed;
public:
    GenotypeView(const Eigen::MatrixXd &m):matrix(&m){}
    GenotypeView(std::shared_ptr<BedColumns> b):bed(std::move(b)){}
    Eigen::Index rows()const{return bed?bed->rows():matrix->rows();}
    Eigen::Index cols()const{return bed?bed->cols():matrix->cols();}
    bool streamed()const{return bool(bed);}
    Eigen::Ref<const Eigen::VectorXd> col(Eigen::Index j)const{return bed?bed->col(j):Eigen::Ref<const Eigen::VectorXd>(matrix->col(j));}
    const Eigen::MatrixXd &dense()const;
};
