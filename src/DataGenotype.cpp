// SPDX-License-Identifier: GPL-3.0-or-later
#include "Data.hpp"
#include "MemoryBudget.hpp"

void Data::readBedColumns(const bool noscale,const string &path,bool molecular,bool materialize) {
    std::vector<unsigned> inds,markers;
    if(molecular){for(unsigned i=0;i<numInds;++i)if(indInfoVec[i]->hasEQTL && indInfoGeneMap.at(indInfoVec[i]->catID)->kept)inds.push_back(i);}
    else {inds.resize(numKeptInds);for(unsigned i=0;i<numInds;++i)if(indInfoVec[i]->kept)inds.at(indInfoVec[i]->index)=i;}
    if(inds.size()!=(molecular?numKeptIndsGene:numKeptInds))throw std::runtime_error("BED phenotype row mapping mismatch");
    for(unsigned j=0;j<numSnps;++j)if(snpInfoVec[j]->included && (!molecular || snpInfoVec[j]->iseQTL))markers.push_back(j);
    auto source=std::make_shared<BedColumns>(path,numInds,numSnps,std::move(inds),markers);
    if(!molecular)source->weights=RinverseSqrt;
    VectorXd variances(markers.size()), norms(markers.size());std::vector<int> keep;
    for(unsigned j=0;j<markers.size();++j) {
        auto *snp=snpInfoVec[markers[j]];VectorXd x=source->raw(j);double sum=0;unsigned missing=0;
        for(Eigen::Index i=0;i<x.size();++i){if(x[i]==-9)++missing;else sum+=x[i];}
        if(missing==x.size())throw std::runtime_error("All genotypes missing at SNP "+snp->rsID);
        const double mean=sum/double(x.size()-missing);source->means[j]=mean;
        for(Eigen::Index i=0;i<x.size();++i)if(x[i]==-9)x[i]=mean;
        const double af=.5*mean, twopq=2*af*(1-af);
        if(molecular) {
            auto *eqtl=eqtlInfoMap.at(snp->rsID);eqtl->af=af;eqtl->twopq=twopq;
            if(twopq==0 || x.maxCoeff()-x.minCoeff()<1e-6){eqtl->included=false;continue;}
        } else {snp->af=af;snp->twopq=twopq;}
        x.array()-=mean;
        double scale=noscale?1.0:std::sqrt(twopq);
        if(!noscale && !molecular && matchedGWAS)scale=std::sqrt(x.squaredNorm()/(x.size()-1.0));
        if(!noscale && !suppliedGenotypeScale.empty())scale=suppliedGenotypeScale.at(snp->rsID);
        if(!(scale>0) || !std::isfinite(scale))throw std::runtime_error("Invalid genotype scale (monomorphic SNP): "+snp->rsID);
        source->scales[j]=scale;if(!molecular)snp->scaleFactor=scale;
        if(!noscale)x.array()/=scale;
        variances[j]=noscale?twopq:Gadget::calcVariance(x);
        // Both GWAS and molecular empirical variances use the sample (N-1)
        // denominator of the supplied R normalization/covariance reference.
        if(!noscale && (matchedGWAS || !suppliedGenotypeScale.empty()))variances[j]*=x.size()/(x.size()-1.0);
        if(!molecular)x.array()*=RinverseSqrt.array();
        norms[j]=x.squaredNorm();keep.push_back(j);
    }
    if(molecular) {
        source->retain(keep);
        if(materialize){ZGene.resize(source->rows(),source->cols());for(Eigen::Index j=0;j<source->cols();++j)ZGene.col(j)=source->col(j);geneBedColumns.reset();}
        else {geneBedColumns=source;ZGene.resize(0,0);}
        snp2pqEqtl=variances(keep);ZPZdiagGene=norms(keep);
        incdEqtlInfoVec=makeIncdEqtlInfoVec(eqtlInfoVec);numIncdEqtls=incdEqtlInfoVec.size();
        if(source->cols()!=cisSnpIDVec.size())throw std::runtime_error("Molecular BED/QC column mapping mismatch");
        ZGeneDat.emplace_back(cisSnpIDVec,geneGenotype());
        std::vector<int> non;for(unsigned j=0;j<numIncdSnps;++j)if(!eqtlInfoMap.count(snpEffectNames[j]) || !eqtlInfoMap.at(snpEffectNames[j])->included)non.push_back(j);
        snp2pqNonEqtl=snp2pq(non);
    } else {
        if(materialize){Z.resize(source->rows(),source->cols());for(Eigen::Index j=0;j<source->cols();++j)Z.col(j)=source->col(j);bedColumns.reset();}
        else {bedColumns=source;Z.resize(0,0);}
        snp2pq=variances;ZPZdiag=norms;
    }
    LOGGER << (molecular?"Molecular":"GWAS") << (materialize?" genotypes: dense matrix; ":" genotypes: BED-backed stream; ") << source->rows() << " individuals x " << source->cols()
           << " SNPs; BED decoder column cache (" << MemoryBudget::bytes(source->rows(),1)/1024 << " KiB), no expanded genotype file." << endl;
}
