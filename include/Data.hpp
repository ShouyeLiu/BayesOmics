//
//  data.hpp
//  gctb
//
//  Created by Jian Zeng on 14/06/2016.
//  Copyright © 2016 Jian Zeng. All rights reserved.
//

#ifndef data_hpp
#define data_hpp
#include <complex>
#define lapack_complex_double std::complex<double>

#include <iostream>
// #include <cstdio>
#include <cstdlib>
#include <fstream>
#include <set>
#include <bitset>
#include <iomanip>     
#include <unordered_map>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
#include <boost/format.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/filesystem.hpp>
#include <boost/foreach.hpp> 
#include <boost/regex.hpp>
#include <omp.h>
// #include <julia.h>

#include "Macro.hpp"
#include "Gadgets.hpp"
#include "Stat.hpp"
#include "Logger.hpp"

using namespace std;
using namespace Eigen;
// using namespace arma;
using namespace boost::multiprecision;
using boost::multiprecision::cpp_dec_float_100;

typedef SparseMatrix<double, Eigen::ColMajor, long long> SpMat;

struct VectorDat
{
public:
    const vector<string> names;
    std::map<string, int> name2index;
    unsigned size;
    Eigen::VectorXd values;

    VectorDat(const vector<string> &names, const Eigen::VectorXd &values)
        : names(names), size(int(names.size())), values(values)
    {
        for (unsigned j = 0; j < size; j++)
            name2index.insert(pair<string, int>(names[j], j));
    }
    double sval(string nameIdx) const { return values(name2index.at(nameIdx)); } // return single value of vector by indexing single name
};

struct MatrixDat
{
public:
    const vector<string> colnames;
    std::map<string, int> colname2index;
    vector<string> rownames;
    std::map<string, int> rowname2index;
    unsigned ncol;
    unsigned nrow;
    Eigen::MatrixXd values;

    MatrixDat(const vector<string> &colnames, const Eigen::MatrixXd &values)
        : colnames(colnames), ncol(int(colnames.size())), values(values)
    {
        nrow = values.rows();
        rownames.resize(nrow);
        for (unsigned j = 0; j < ncol; j++)
            colname2index.insert(pair<string, int>(colnames[j], j));
        for (unsigned j = 0; j < nrow; j++)
        {
            rownames[j] = "row" + to_string(j);
            rowname2index.insert(pair<string, int>(rownames[j], j));
        }
    }
    MatrixDat(vector<string> &rownames, const vector<string> &colnames, const Eigen::MatrixXd &values)
        : colnames(colnames), ncol(int(colnames.size())), rownames(rownames), nrow(int(rownames.size())), values(values)
    {
        for (unsigned j = 0; j < ncol; j++)
            colname2index.insert(pair<string, int>(colnames[j], j));
        for (unsigned j = 0; j < nrow; j++)
            rowname2index.insert(pair<string, int>(rownames[j], j));
    }
    Eigen::VectorXd col(string nameIdx) const { return values.col(colname2index.at(nameIdx)); }
    Eigen::VectorXd row(string nameIdx) const { return values.row(rowname2index.at(nameIdx)); }
};

class AnnoInfo;

class EqtlInfo
{
public:
    const int chrom;
    const string rsID;
    string alterId; // used for alternative names; 
    double genPos;
    const int physPos;
    const string a1; // ALT (effect allele). like plink
    const string a2; // REF (other)
    double af; // the frequency of the effect allele
    double twopq;
    int index;
    bool included; // flag for inclusion in panel
    bool flipped;  // flag for eqtl flip when eqtl->a1 == snp->a2
    // eQTL summary statistics
    double eqtl_b;
    double eqtl_se;
    double eqtl_n;
    double eqtl_chisq;
    double eqtl_af;
    int eqtl2gene_count; // how many genes have this eqtl;
    vector<double> afVec; // pleiotropic eqtl have different af for different gene;
    vector<double> samVec; 
    // ld info
    int ld_n; // sample size for eqtl ld frequency
    // double eqtl_af;
    double eqtl_pvalue;
    bool isInGene; // this will be used in generate gene eigen matrix; default value is false;
    EqtlInfo(const int idx, const string &id, const string &allele1, const string &allele2,
             const int chr, const double gpos, const int ppos, const double frequency)
        : rsID(id), index(idx), a1(allele1), a2(allele2), chrom(chr), genPos(gpos), physPos(ppos), af(frequency)
    {
        included = true;
        isInGene = false;
        flipped = false;
        eqtl_n = -999;
        twopq = -999;
        eqtl_b = -999;
        eqtl_se = -999;
        eqtl_n = -999;
        eqtl_af = -999;
        eqtl_chisq = 0.0;
        eqtl_pvalue = 1.0;
        eqtl2gene_count = 0;
    };
};

class GeneInfo
{
public:
    const int chrom;
    string probeID; // can be the ID of an exon or a transcript for RNA-seq data
    double geneticDis; // genetic distance (can be any arbitary value)
    double phyPos;
    const string ensemblID;  // unique identifier
    string geneOri; // gene orientation
    int start;
    int end;
    double midPhyPos;
    int geneLength;
    int index;
    bool kept;
    int windows;
    long sampleSize;
    double sampleSizeSD;

    // eQTL information related
    bool hasEqtl; // this will make sure gene have at least one eqtls; default value is false;
    vector<string> cisSnpNameVec; // store information from based line as a bechmark
    map<string,int> cisSnpSampleSizeMap;
    map<string,int> cisSnpID2IdxMapInGene;
    set<string> cisSnpIDSetInGene;
    VectorXd  eQTLMarginEffect;   // eQTL marginal effect
    VectorXd eQTLMarginEffectSE; // eQTL standard error of marginal effect
    VectorXd gwasMEffWithinGene;
    vector<double> eQTLChisqVec;
    vector<double> eQTLPvalVec;

    vector<int> gene2GwasSnpVec; // store snps that belong to this gene
    vector<int> gene2CisSnpVec; 
    int numSnpInGene;

    vector<string> cisSnpOrderedByGWASLD; // store snplist index reordered based on gwas LD order
    // impute eqtl or not
    // vector<string> cisSnpInBesdNotInEigenLDVec;
    // here we need to add a new variable since we read gene LD first and then read besd format.
    // vector<string> cisSnpNameVecBESD; // store information from based line as a bechmark
    // map<string,int> cisSnpID2IdxMapInGeneBESD;
    vector<int> typSnpIdxBesd;  // typed snpidx in besd snplist
    vector<int> typSnpIdxGeneLD; // typed snpidx in gene eigen ld snplist
    vector<int> impSnpIdxGeneLD; // impe snpidx in gene eigen ld snplist
    bool impCisSnpEffBool;

    // eigenValue and eivenVector
    MatrixXd eigenVector;
    MatrixXd eigenValue;


    GeneInfo(const int idx, const string &gene, const int chr)
        : index(idx), ensemblID(gene), chrom(chr)
    {
        probeID = "";
        geneticDis = -999;
        phyPos = -999;
        midPhyPos = - 999;
        geneOri = "NA";
        start = -999;
        end = -999;
        geneLength = -999;
        kept = true;
        hasEqtl = false;
        impCisSnpEffBool = false;
        windows = -999;
        sampleSize = -999;
        sampleSizeSD = -999;
        numSnpInGene =0;
        impCisSnpEffBool = false;

    }
};
 
class SnpInfo {
public:
    const string rsID;
    const string a1; // ALT (effect)
    const string a2; // REF (other)
    const int chrom;
    double genPos;
    const int physPos;
    
    
    int index;
    int window;
    int windStart;  // for window surrounding the SNP
    int windSize;   // for window surrounding the SNP
    int windEnd;
    int windStartOri;  // original value from .info file
    int windSizeOri;
    int windEndOri;
    double af;       // allele frequency of a1
    double twopq;
    bool included;  // flag for inclusion in panel
    bool isQTL;     // for simulation
    bool recoded;   // swap A1 and A2: use A2 as the reference allele and A1 as the coded allele
    bool skeleton;  // skeleton snp for sbayes
    bool flipped;   // A1 A2 alleles are flipped in between gwas and LD ref samples
    bool isInLD;
    long sampleSize;
    
    string block;
    
    VectorXd genotypes; // temporary storage of genotypes of individuals used for building sparse Z'Z
    
    vector<AnnoInfo*> annoVec;
    vector<unsigned> annoIdx;   // the index of SNP in the annotation
    map<int, AnnoInfo*> annoMap;
    
    VectorXd annoValues;
    
    double effect;   // estimated effect
    
    // GWAS summary statistics
    double gwas_b;
    double gwas_se;
    double gwas_n;
    double gwas_af;
    double gwas_chisq;
    double gwas_pvalue;

    double ldSamplVar;    // sum of sampling variance of LD with other SNPs for summary-bayes method
    double ldSum;         // sum of LD with other SNPs
    double ldsc;          // LD score: sum of r^2
    double scaleFactor; // like snp2pq, but nore stable

    int ld_n;  // LD reference sample size 
    int numNonZeroLD;   // may be different from windSize in shrunk ldm
    unsigned numAnnos;
    // eqtl summary statistics
    vector<string> geneNameVec;
    vector<int> gwasSnp2geneVec;     // store genes that harbor this snp
    bool iseQTL;                     // for eQTL, this is used in SBayesE.
    bool isInBlock; // judge whether the snp is located in ld block or not

    SnpInfo(const int idx, const string &id, const string &allele1, const string &allele2,
            const int chr, const double gpos, const int ppos)
    : rsID(id), index(idx), a1(allele1), a2(allele2), chrom(chr), genPos(gpos), physPos(ppos) {
        window = 0;
        windStart = -1;
        windSize  = 0;
        windEnd   = -1;
        windStartOri = -1;
        windSizeOri = 0;
        windEndOri = -1;
        af = -999;
        twopq = -1;
        included = true;
        isQTL = false;
        iseQTL = false;
        isInBlock = false;
        recoded = false;
        skeleton = false;
        flipped = false;
        isInLD = false;
        sampleSize = 0;
        effect = 0;
        gwas_b  = -999;
        gwas_se = -999;
        gwas_n  = -999;
        gwas_af = -999;
        gwas_chisq = 0.0;
        gwas_pvalue = 1.0;
        ldSamplVar = 0.0;
        ldSum = 0.0;
        ldsc = 0.0;
        numNonZeroLD = 0;
        numAnnos = 0;
        ld_n = -999;
        scaleFactor = 1;
        block = "NA";
    };
    
    void resetWindow(void) {windStart = -1; windSize = 0;};
    bool isProximal(const SnpInfo &snp2, const double genWindow) const;
    bool isProximal(const SnpInfo &snp2, const unsigned physWindow) const;
};

class LDBlockInfo {
public:
    const string ID;
    const int chrom;
    int index;
    // block
    int startPos;
    int endPos;
    double start_cm;
    double stop_cm;
    int pdist;
    double gdist;
    bool kept;

    // the following variables aim to read svd-ld matrix from SBayesRC-Eigen
    int startSnpIdx;
    int endSnpIdx;
    
    //int idxStart;
    //int idxEnd;
    int preBlock;
    int postBlock;
    //
    vector<string> snpNameVec;
    vector<SnpInfo*> snpInfoVec;
    vector<int> block2GwasSnpVec; // store snps that belong to this block;
    vector<bool> block2EqtlBoolVec; //  given snp is eqtl or not.
    int numSnpInBlock;
    
    VectorXd eigenvalues;
    float sumPosEigVal; // sum of all positive eigenvalues

    LDBlockInfo(const int idx, const string id, const int chr) : index(idx), ID(id), chrom(chr)
    {
        // block info
        startPos = -999;
        endPos = -999;
        start_cm = -999;
        stop_cm = -999;
        pdist = -999;
        gdist = -999;
        // svd'ed ld info
        startSnpIdx = -999;
        endSnpIdx = -999;
        //idxStart = -999;
        //idxEnd = -999;
        preBlock = -999;
        postBlock = -999;
        numSnpInBlock = 0;
        kept = true;
        sumPosEigVal = 0;
    }
};

class locus_bp {
public:
    string locusName;
    int chr;
    double bp;

    locus_bp(string locusNameBuf, int chrBuf, double bpBuf)
    {
        locusName = locusNameBuf;
        chr = chrBuf;
        bp = bpBuf;
    }

    bool operator()(const locus_bp &other)
    {
        return (chr == other.chr && bp <= other.bp);
    }
};

class ChromInfo {
public:
    const int id;
    const unsigned size;
    const int startSnpIdx;
    const int endSnpIdx;
    const unsigned startSnPhyPos;
    const unsigned endSnpPhyPos;
    string startSnpID;
    string endSnpID;

    
    ChromInfo(const int id, const unsigned size, const int startSnp, const int endSnp, const unsigned startSnPhyPos,const unsigned endSnpPhyPos):
     id(id), size(size), startSnpIdx(startSnp), endSnpIdx(endSnp),startSnPhyPos(startSnPhyPos),endSnpPhyPos(endSnpPhyPos) {
        startSnpID = "";
        endSnpID = "";
     }
};

class AnnoInfo {  // annotation info for SNPs
public:
    int idx;
    const string label;
    unsigned size;
    double fraction;   // fraction of all SNPs in this annotation
    unsigned chrom;   // for continuous annotation
    unsigned startBP; // for continuous annotation
    unsigned endBP;   // for continuous annotation
    
    vector<SnpInfo*> memberSnpVec;
    map<int, SnpInfo*> memberSnpMap;
    VectorXd snp2pq;
    
    AnnoInfo(const int idx, const string &lab): idx(idx), label(lab){
        size = 0;
        chrom = 0;
        startBP = 0;
        endBP = 0;
    }
    
    void getSnpInfo(void);
    void print(void);
};

class IndInfo {
public:
    const string famID;
    const string indID;
    const string catID;    // catenated family and individual ID
    const string fatherID;
    const string motherID;
    const int famFileOrder; // original fam file order
    const int sex;  // 1: male, 2: female
    
    int index;
    bool kept;
    
    double phenotype;
    double rinverse;
    
    VectorXd covariates;  // covariates for fixed effects
    VectorXd randomCovariates;

    std::map<string,double> genePheMap; // gene expression (protein ) phenotype;
    bool hasEQTL;                     // for eQTL, this is used in SBayesE.

    IndInfo(const int idx, const string &fid, const string &pid, const string &dad, const string &mom, const int sex)
    : famID(fid), indID(pid), catID(fid+":"+pid), fatherID(dad), motherID(mom), index(idx), famFileOrder(idx), sex(sex) {
        phenotype = -9;
        rinverse = 1;
        kept = true;
        hasEQTL = false;
    }
};

class Data {
public:
    MatrixXd X;              // coefficient matrix for fixed effects
    MatrixXd W;              // coefficient matrix for random effects
    MatrixXd Z;              // coefficient matrix for SNP effects
    VectorXd D;              // 2pqn
    VectorXd y;              // phenotypes
    
    //SpMat ZPZ; // sparse Z'Z because LE is assumed for distant SNPs
    vector<VectorXd> ZPZ;
    MatrixXd ZPZmat;
    vector<SparseVector <double>  > ZPZsp;
    SpMat ZPZspmat;
    SpMat ZPZinv;
    
    MatrixXd annoMat;        // annotation coefficient matrix
    MatrixXd APA;            // annotation X'X matrix
    VectorXd annoMean;       // column mean of annotation coefficient matrix
    VectorXd annoSD;         // column SD of annotation coefficient matrix

    MatrixXd XPX;            // X'X the MME lhs
    MatrixXd WPW;
    MatrixXd ZPX;            // Z'X the covariance matrix of SNPs and fixed effects
    VectorXd XPXdiag;        // X'X diagonal
    VectorXd WPWdiag;
    VectorXd ZPZdiag;        // Z'Z diagonal
    VectorXd XPy;            // X'y the MME rhs for fixed effects
    VectorXd ZPy;            // Z'y the MME rhs for snp effects
    
    VectorXd snp2pq;         // 2pq of SNPs
    VectorXd scalingGWASFactorVec; // we use this to scaling snp effect assuming 1 phenotypic variance
    VectorXd se;             // se from GWAS summary data
    VectorXd tss;            // total ss (ypy) for every SNP
    VectorXd b;              // beta from GWAS summary data
    VectorXd n;              // sample size for each SNP in GWAS
    VectorXd Dratio;         // GWAS ZPZdiag over reference ZPZdiag for each SNP
    VectorXd DratioSqrt;     // square root of GWAS ZPZdiag over reference ZPZdiag for each SNP
    VectorXd chisq;          // GWAS chi square statistics = D*b^2
    VectorXd varySnp;        // per-SNP phenotypic variance
    
    VectorXi windStart;      // leading snp position for each window
    VectorXi windSize;       // number of snps in each window

    // for Eigen dec
    VectorXi blockStarts;    // each LD block startings index in SNP included scale
    VectorXi blockSizes;     // each LD block size;
    VectorXd nGWASblock;     // median GWAS sample size for each block in GWAS
    VectorXd numSnpsBlock;   // number of SNPs for each block
    VectorXd numEigenvalBlock;  // number of eigenvalues kept for each block
    
    VectorXd LDsamplVar;     // sum of sampling variance of LD for each SNP with all other SNPs; this is for summary-bayes methods
    VectorXd LDscore;        // sum of r^2 over SNPs in significant LD
    
    VectorXd RinverseSqrt;   // sqrt of the weights for the residuals in the individual-level model
    VectorXd Rsqrt;
    
    double ypy;               // y'y the total sum of squares adjusted for the mean
    double varGenotypic;
    double varResidual;
    double varPhenotypic;
    double varRandom;         // variance explained by random covariate effects
    
    bool reindexed;
    bool sparseLDM;
    bool shrunkLDM;
    bool readLDscore;
    bool makeWindows;
    bool weightedRes;
    
    bool lowRankModel;
    
    vector<SnpInfo*> snpInfoVec;
    vector<IndInfo*> indInfoVec;

    vector<AnnoInfo*> annoInfoVec;
    vector<string> annoNames;
    vector<string> snpAnnoPairNames;
    
    map<string, SnpInfo*> snpInfoMap;
    map<string, IndInfo*> indInfoMap;

    vector<SnpInfo*> incdSnpInfoVec;
    vector<IndInfo*> keptIndInfoVec;
    
    vector<string> fixedEffectNames;
    vector<string> randomEffectNames;
    vector<string> snpEffectNames;
    
    set<int> chromosomes;
    vector<ChromInfo*> chromInfoVec;
    map<int, ChromInfo *> chromInfoMap;

    vector<bool> fullSnpFlag;    
    vector<unsigned> numSnpMldVec;
    vector<unsigned> numSnpAnnoVec;
    VectorXd numAnnoPerSnpVec;
    vector<SpMat> annowiseZPZsp;
    vector<VectorXd> annowiseZPZdiag;
    vector<vector<unsigned> > windowSnpIdxVec;
    
    //////// ld block begin ///////
     vector<LDBlockInfo *> ldBlockInfoVec;
     vector<LDBlockInfo *> keptLdBlockInfoVec;
     map<string, LDBlockInfo *> ldBlockInfoMap;
     vector<string> ldblockNames;
     vector<VectorXd> eigenValLdBlock; // store lambda  (per LD block matrix = U * diag(lambda)* V')  per gene LD
     vector<MatrixXd> eigenVecLdBlock; // store U   (per  LD block matrix = U * diag(lambda)* V')  per gene LD
     vector<VectorXd> wcorrBlocks;
     vector<MatrixXd> Qblocks;
     vector<MatrixDat> QblocksDat;
     ///////// ld block end  ////////
    map<string,MatrixXd> geneEigenVecMap;
    map<string,VectorXd> geneEigenValMap;
    map<int, vector<int>> ldblock2gwasSnpMap;
    vector<VectorXd> gwasEffectInBlock;  // gwas marginal effect;
    vector<VectorXd> pseudoGwasEffectTrn;
    vector<VectorXd> pseudoGwasEffectVal;
    double pseudoGwasNtrn;
    VectorXd b_val;


    //////////// BayesOmics parameter begin /////////
    vector<IndInfo *> indInfoGeneVec;  // store phenotype from gene, different from indInfoVec;
    vector<IndInfo *> keptindInfoGeneVec;
    vector<EqtlInfo *> eqtlInfoVec;
    map<string, EqtlInfo *> eqtlInfoMap;
    map<string, EqtlInfo *> eqtlInfoAlterMap; // used for eqtl alternative names;
    vector<EqtlInfo *> incdEqtlInfoVec;
    vector<GeneInfo *> geneInfoVec;
    map<string, GeneInfo *> geneInfoMap;
    vector<GeneInfo *> keptGeneInfoVec;
    
    // double piEffNonEqtlVal;
    // double piEffEqtlVal;
    vector<VectorDat> neQTLVec;                // sample size for each SNP in eQTL
    vector<VectorXd> snp2pqeQTL;               // 2pq for each eQTL
    vector<VectorXd> scalingeQTLFactorVecVec;  // // we use this to scaling eQTL effect assuming 1 phenotypic variance
    VectorXd varGenotypiceQTL;                 // genetic variance of eQTL
    VectorXd varPhenotypiceQTL;                // variance of gene expression for each gene
    VectorXd varResidualeQTL;
    VectorXd varRandomeQTL;
    VectorXd ypyeQTL;

    vector<string> cisSnpIDVec; // Initiate in function makeIncdEqtlInfoVec();
    vector<string> geneEffectNames;
    vector<string> gwasAndGeneEffectNames; // used in SBayesE mixed model to incorporate SBayesC and SBayesE

    // Variables used in summary part.
    vector<VectorXd> eQTLEffAcrossGenes;
    VectorXd eQTLEffMeanAcrossGenes;
    VectorXd eQTLEffSEMeanAcrossGenes;
    vector<double> sampleSizeAcrossGene; // used to summarise per snp sample size across genes;
    VectorXd numSnpsInEigenAcrossGene;   // number of SNPs for each gene 
    VectorXd numEigenvalAcrossGene;  // number of eigenvalues kept for each gene

    map<string, IndInfo *> indInfoGeneMap;
    map<string, vector<int> > genePheIdxMap;
    map<int, vector<int>> gene2cisSnpMap; // this will show the location of eQTLs in cis-region, used in BayesE.
    map<string , int> cisSnpID2IdxMap;
    map<int, string> gwasSnpIdx2snpIDMap;
    map<string, int> geneID2IdxMap;  // geneID -> geneIdx
    map<string, vector<string>> gwasSnpID2geneIDMap; // use for switching SBayesE or SBayesC
    map<string, vector<int> > gwasSnpID2geneIdxMap;
    map<int, vector<int>> gene2gwasSnpMap; // this will show the location of eQTL in the whole snplist.
    map<int, vector<string>> gene2cisSnpIDMap;
    map<int, vector<string>> gwas2SnpIDMap; // use in SBayesCO class


    int numKeptIndsGene;
    MatrixXd ZGene; // coefficient matrix for SNP effects from cis-region
    VectorXd ZPZdiagGene; // Z'Z diagonal
    VectorXd snp2pqEqtl;     // 2pq of SNPs
    VectorXd snp2pqNonEqtl;         // 2pq of SNPs
    VectorXd wbcorr;
    
    vector <VectorXd > genePheVec;
    vector<VectorXd> wAcorr; 
    vector<VectorXd> wbcorrGene;
    vector<MatrixDat> ZGeneDat;
    vector<MatrixDat> ZDat;

    unsigned numIndGenes;
    unsigned numIncdEqtls;
    unsigned numKeptGenes;
    unsigned numGenes;
    unsigned numeQTLs;
    unsigned numKeptIndseQTL;
    unsigned numNonEqtl;
    unsigned numEqtl;
    unsigned numEqtlOverlap;
    unsigned numEqtlFlip;
    //////////// BayesOmics parameter end ///////////

    //////////// SBayesOmics begin //////////////
    vector<MatrixXd> rvalGene;
    vector<VectorXd> eigenValGene; // store lambda  (per gene LD matrix = U * diag(lambda)* V')  per gene LD
    vector<MatrixXd> eigenVecGene; // store U   (per gene LD matrix = U * diag(lambda)* V')  per gene LD
    map<string,float> geneEigenSumPosEigValMap;
    map<string,float> geneEigenCutoffMap;
    // vector<MatrixDat> Qgene;
    vector<MatrixDat> QgeneDat;


    // // save ind names for debug
    vector<string> gwasIndNames;
    vector<string> eqtlIndNameVec;
    // //////////// SBayesOmics end /////////////////
    unsigned numFixedEffects;
    unsigned numRandomEffects;
    unsigned numSnps;
    unsigned numInds;
    unsigned numIncdSnps;
    unsigned numKeptInds;
    unsigned numChroms;
    unsigned numSkeletonSnps;
    unsigned numAnnos;
    unsigned numWindows;
    unsigned numLDBlocks;
    unsigned numKeptLDBlocks;
    unsigned numGWASFlip;
    
    string label;
    string title;
    
    int niter;
    int burnIn;

    Data(){
        numFixedEffects = 0;
        numRandomEffects = 0;
        numSnps = 0;
        numInds = 0;
        numIncdSnps = 0;
        numKeptInds = 0;
        numChroms = 0;
        numSkeletonSnps = 0;
        numAnnos = 0;
        numWindows = 0;
        numLDBlocks = 0;
        numKeptLDBlocks = 0;
        numGWASFlip = 0;

        numIndGenes = 0;
        numIncdEqtls = 0;
        numKeptGenes = 0;
        numGenes = 0;
        numeQTLs = 0;
        numKeptIndseQTL = 0;
        numNonEqtl = 0;
        numEqtl = 0;
        numEqtlOverlap = 0;
        numEqtlFlip = 0;

        niter = 0;
        burnIn = 0;

        reindexed = false;
        sparseLDM = false;
        readLDscore = false;
        makeWindows = false;
        weightedRes = false;
        lowRankModel = false;
    }
    
    void readFamFile(const string &famFile);
    void readBimFile(const string &bimFile);
    void readBedFile(const bool noscale, const string &bedFile);
    void readPhenotypeFile(const string &phenFile, const unsigned mphen = 1);
    void readCovariateFile(const string &covarFile);
    void readRandomCovariateFile(const string &covarFile);
    void readGwasSummaryFile(const string &gwasFile, const double afDiff, const double mafmin, const double mafmax, const double pValueThreshold, const bool imputeN, const bool removeOutlierN);
    void readLDmatrixInfoFileOld(const string &ldmatrixFile);
    void readLDmatrixInfoFile(const string &ldmatrixFile);
    void readLDmatrixBinFile(const string &ldmatrixFile);
    void readLDmatrixTxtFile(const string &ldmatrixFile);
    void readGeneticMapFile(const string &freqFile);
    void readfreqFile(const string &geneticMapFile);
    void keepMatchedInd(const string &keepIndFile, const unsigned keepIndMax,bool ldBool= false);
    void includeSnp(const string &includeSnpFile);
    void excludeSnp(const string &excludeSnpFile);
    void includeChr(const unsigned chr);
    void includeChreQTL(const unsigned chr);
    void includeChrGene(const unsigned chr);
    void includeBlock(const unsigned block);
    void includeBlockID(const string includeBlockID);
    void excludeMHC(void);
    void excludeAmbiguousSNP(void);
    void excludeSNPwithMaf(const double mafmin, const double mafmax);
    void excludeRegion(const string &excludeRegionFile);
    void includeSkeletonSnp(const string &skeletonSnpFile);

    void includeMatchedSnp(void);
    vector<SnpInfo*> makeIncdSnpInfoVec(const vector<SnpInfo*> &snpInfoVec);
    vector<IndInfo*> makeKeptIndInfoVec(const vector<IndInfo*> &indInfoVec);
    void getWindowInfo(const vector<SnpInfo*> &incdSnpInfoVec, const unsigned windowWidth, VectorXi &windStart, VectorXi &windSize);
    void getNonoverlapWindowInfo(const unsigned windowWidth);
    void buildSparseMME(const string &bedFile, const unsigned windowWidth);
//    void makeLDmatrix(const string &bedFile, const unsigned windowWidth, const string &filename);
    string partLDMatrix(const string &partParam, const string &outfilename, const string &LDmatType);
    void makeLDmatrix(const string &bedFile, const string &LDmatType, const double chisqThreshold, const double LDthreshold, const unsigned windowWidth,
                      const string &snpRange, const string &filename, const bool writeLdmTxt);
    void makeshrunkLDmatrix(const string &bedFile, const string &LDmatType, const string &snpRange, const string &filename, const bool writeLdmTxt, const double effpopNE, const double cutOff, const double genMapN);
    void resizeWindow(const vector<SnpInfo*> &incdSnpInfoVec, const VectorXi &windStartOri, const VectorXi &windSizeOri,
                      VectorXi &windStartNew, VectorXi &windSizeNew);
    void computeAlleleFreq(const MatrixXd &Z, vector<SnpInfo*> &incdSnpInfoVec, VectorXd &snp2pq);
    void reindexSnp(vector<SnpInfo*> snpInfoVec);
    void initVariances(const double heritability, const double propVarRandom);
    
    void outputSnpResults(const VectorXd &posteriorMean, const VectorXd &posteriorSqrMean, const VectorXd &pip, const bool noscale, const string &filename) const;
    void outputSnpResults(const VectorXd &posteriorMean, const VectorXd &posteriorSqrMean, const VectorXd &lastSample, const VectorXd &pip, const bool noscale, const string &filename) const;
    void outputGeneEffectResults(const VectorXd &posteriorMean, const VectorXd &posteriorSqrMean, const VectorXd &lastSample,const VectorXd &pip,const string &mcmcType,const string &filename) const;
    // void outputGeneParameter(const VectorXd &geneEffMean, const VectorXd &geneEffSqrMean,
    //                               const VectorXd &cisMean, const VectorXd &cisSqrMean,
    //                               const VectorXd &gwasCisMean, const VectorXd &gwasCisSqrMean,
    //                               const VectorXd &genicGwasEnrichMean, const VectorXd &genicGwasEnrichSqrMean,
    //                               const VectorXd &genicCisEnrichMean, const VectorXd &genicCisEnrichSqrMean,
    //                               const string &filename) const;
    void outputFixedEffects(const MatrixXd &fixedEffects, const string &filename) const;
    void outputRandomEffects(const MatrixXd &randomEffects, const string &filename) const;
    void outputWindowResults(const VectorXd &posteriorMean, const string &filename) const;
    void summarizeSnpResults(const SpMat &snpEffects, const string &filename) const;
    void buildSparseMME(const bool sampleOverlap, const bool noscale);
    void readMultiLDmatInfoFile(const string &mldmatFile);
    void readMultiLDmatBinFile(const string &mldmatFile);
    void outputSnpEffectSamples(const SpMat &snpEffects, const unsigned burnin, const unsigned outputFreq, const string &snpResFile, const string &filename) const;
    void resizeLDmatrix(const string &LDmatType, const double chisqThreshold, const unsigned windowWidth, const double LDthreshold, const double effpopNE, const double cutOff, const double genMapN);
    void outputLDmatrix(const string &LDmatType, const string &filename, const bool writeLdmTxt) const;
    void displayAverageWindowSize(const VectorXi &windSize);
    
    void inputSnpResults(const string &snpResFile);
    void inputSnpInfoAndResults(const string &snpResFile, const string &bayesType);
    void readLDmatrixBinFileAndShrink(const string &ldmatrixFile);
    void readMultiLDmatBinFileAndShrink(const string &mldmatFile, const double genMapN);
    void directPruneLDmatrix(const string &ldmatrixFile, const string &outLDmatType, const double chisqThreshold, const string &title, const bool writeLdmTxt);
    void jackknifeLDmatrix(const string &ldmatrixFile, const string &outLDmatType, const string &title, const bool writeLdmTxt);
    void addLDmatrixInfo(const string &ldmatrixFile);
    
    void readLDscoreFile(const string &ldscFile);
    void imputePerSnpSampleSize(vector<SnpInfo*> &snpInfoVec, unsigned &numIncdSnps, double sd);
    void getZPZspmat(void);
    void getZPZmat(void);
    void binSnpByLDrsq(const double rsqThreshold, const string &title);
    void readWindowFile(const string &windowFile);
    void binSnpByWindowID(void);
    void filterSnpByLDrsq(const double rsqThreshold);
    void readResidualDiagFile(const string &resDiagFile);
    
    void mergeLdmInfo(const string &outLDmatType, const string &dirname);

    /////////// eigen decomposition for LD blocks
    void readLDBlockInfoFile(const string &ldBlockInfoFile);
    void getEigenDataFromFullLDM(const string &filename, const float eigenCutoff = 0.995);
    void eigenDecomposition( const MatrixXf &X, const float &prop, VectorXf &eigenValAdjusted, MatrixXf &eigenVecAdjusted, float &sumPosEigVal);
    MatrixXd generateLDmatrixPerBlock(const string &bedFile, const vector<string> &snplists); // generate full LDM for block using genotype directly
    MatrixXd generateLDmatrixPerBlock(const vector<int> &snplistIdx); // generate full LDM for block using genotype Z matrix
    

    void makeBlockLDmatrix(const string &bedFile, const string &LDmatType, const unsigned block, const string &filename, const bool writeLdmTxt, int ldBlockRegionWind = 0);
    void readBlockLdmBinaryAndDoEigenDecomposition(const string &dirname, const unsigned block, const float eigenCutoff, const bool writeLdmTxt);
    // get eigen matrix from LD matrix
    void getEigenDataForLDBlock(const string &bedFile, const string &ldBlockInfoFile, int ldBlockRegionWind, const string &filename, const float eigenCutoff);
    // get eigen matrix from genotype
    void calcBlockLDEigenDecom(const string &bedFile, const string &ldBlockInfoFile, int ldBlockRegionWind, const string &filename, const double eigenCutoff = 0.9995,bool outInfoOnly = false);
    void outputBlockLDmatrixInfo(const LDBlockInfo &block, const string &outSnpfile, const string &outldmfile) const;
    void impG(const unsigned block, double diag_mod = 0.1);
    ///////////// read LD matrix eigen-decomposition data for LD blocks
    void readEigenMatrix(const string &eigenMatrixFile, const double eigenCutoff = 0.9995,const bool qcBool = false);
    vector<LDBlockInfo *> makeKeptLDBlockInfoVec(const vector<LDBlockInfo *> &ldBlockInfoVec);
    void readEigenMatrixBinaryFile(const string &eigenMatrixFile, const double eigenCutoff);
    void readEigenMatrixBinaryFileAndMakeWandQ(const string &dirname, const double eigenCutoff, const vector<VectorXd> &GWASeffects, const double nGWAS, const bool noscale, const bool makePseudoSummary);
    void UseLDBlockEigenMakeWAndQblocks(const vector<VectorXd> &GWASeffects, const double nGWAS, const bool noscale,const bool makePseudoSummary = true);
    ///////////// merge eigen matrices
    void mergeMultiEigenLDMatrices(const string & infoFile, const string &filename, const string LDmatType);
    //////////// Step 2.2 Build multiple maps
    void buildMMEigen(const string &dirname, const bool sampleOverlap, const double eigenCutoff, const bool noscale, const bool gwasInfoOnly); // for eigen decomposition
    void includeMatchedBlocks(const bool haveEqtlInfo = false);
    void readEqtlSummaryFileFromLD2BESD(const string &besdFile, const string &eqtlSummaryQueryFile,const unsigned includeChr, const double afDiff, const double mafmin, const double mafmax, const double pValueThreshold, const bool imputeN, const double eigenCutoff);
    //////////// Step 2.3 build model matrix
//    void constructWandQ(const double eigenCutoff, const bool noscale);
 
    void imputeSummaryData(void);
    void truncateEigenMatrix(const float sumPosEigVal, const double eigenCutoff, const VectorXd &oriEigenVal, const MatrixXd &oriEigenVec, VectorXd &newEigenVal, MatrixXd &newEigenVec);
    void constructPseudoSummaryData(void);
    void constructWandQ(const vector<VectorXd> &GWASeffects, const double nGWAS, const bool noscale);
    void scaleGwasEffects(void);
    void mapSnpsToBlocks(int ldBlockRegionWind = 0);
    void mergeBlockGwasSummary(const string &gwasSummaryFile, const string &title);

    ///////////// Individual-level BayesOmics function begin ///////////
    vector<EqtlInfo *> makeIncdEqtlInfoVec(const vector<EqtlInfo *> &eqtlInfoVec);
    void makeIncdEqtlInfoVecBasedOnGeneInfoVec();
    vector<GeneInfo *> makeKeptGeneInfoVec(const vector<GeneInfo *> &geneInfoVec);
    vector<IndInfo*> makeKeptIndInfoGeneVec(const vector<IndInfo*> &indInfoGeneVec);
    void includeMatchedEqtl();
    void initeQTLVariances(const double heritability, const double propVarRandom);
    void readPhenoFromGeneFile(const string geneID, const string & singleGenePheFile,const string subGenePath);
    void readGeneInfoFile(const string geneFile,const string subGenePath);
    void keepMatchedIndGene(const string &keepIndFile, const unsigned keepIndMax);
    void mapGwasSnpToGeneCisRegion(const string &bedFile, const bool noscale, const double cisRegionWind);
    void readBedFileForGene (const string &bedFile, const bool noscale);
    void ConstructGwasEqtlGeneMaps();
    void ConstructGenePheAndwAcorr();
    void buildBayesOmicsMME(const string bedFile, const bool noscale, const bool haveGene);

    ///////////// data management begin ///////////
    void readSNPRSIDMaps(const string snpRsidMapFile, const string grchType,const string title);
    void readProMaps(const string proteinMapFile,const double pValueThreshold, const string title);
    void readGWASFromProtein(const string geneID, vector<string> & cisSnpNameVec, vector<double>  &eQTLMarginEffect, vector<double> &eQTLMarginEffectSE, map<string,int> &cisSnpID2IdxMapInGene,  set<string> & cisSnpIDSetInGene, const string proPath,const double pValueThreshold, bool &hasEqtl);
    void readFlistFile(const string geneInfoFile,const string subGenePath);
    void readSingleEsdFile(const string geneID,const string singleGenePheFile, int eqtlIdx,const string subGenePath);
    void includexQTLSnp(const string &includeSnpFile);

    void saveBesdFormat(const string title,const bool makeBesdSmrBool);
    void saveAnnoPlainMatFormat(const string title, const bool isAnnoBinary,const bool hasSnpInfo);
    void saveQueryMoQTLInfo(const string title);
    void saveSumstatsFormat(const string title);
    void saveQueryCojoMaMoQTLInfo(const string title, const int slurmArrayLimit);
    void filterMoQTLBasedOnCisWindow(const double cisRegionWind);
    void filterMoQTLBasedOnPvalue(const double pValueThreshold);
    // void filterMoQTLBasedPerSNPSampleSize(const double CoeffSD);
    void mergeMultiBesdData(const string besdListFile,const string title);
    bool readMultiEpiFile(const string &epiFile, vector<GeneInfo *> &geneInfoVecLocal,map<string, GeneInfo *> &geneInfoMapLocal,const bool hasSnpInfo, const bool hasGeneLdmInfo);
    bool readMultiEsiFile(const string &esiFile,vector<EqtlInfo *> &eqtlInfoVecLocal,map<string, EqtlInfo *> &eqtlInfoMapLocal,bool &smrBesdBool,const bool hasSnpInfo, const bool hasGeneLdmInfo);
    bool readMultiBesdFile(const string &besdFile,vector<GeneInfo *> &geneInfoVecLocal,map<string, GeneInfo *> &geneInfoMapLocal,
                            vector<EqtlInfo *> &eqtlInfoVecLocal,map<string, EqtlInfo *> &eqtlInfoMapLocal, bool &smrBesdBool,const bool hasSnpInfo,const bool hasGeneLdmInfo );

    void reSampleXqtlEffect(const int reSamType,const string &besdFile,const string &geneEigenMatrixFile,const int aimN, const double eigenCutoff,const double diag_mod  = 0.1);
    void reSampleEffectByImputation( VectorXd &betaImp, VectorXd &seImp,const int aimN, const VectorXd currN,const VectorXd snp2pq, const VectorXd beta, const VectorXd se, MatrixXd localLD);
    void reSampleEffectByBivariateLDSC( VectorXd &betaImp, VectorXd &seImp,const int aimN, const VectorXd currN,const VectorXd snp2pq, const VectorXd beta, const VectorXd se, MatrixXd localLD);
    ///////////// Step 2.3 read LD matrix of gene regions info
    void readEigenMatLDBlockInfoFile(const string &infoFile);
    void readEigenMatLDBlockSnpInfoFile(const string &snpInfoFile,const string LDMatType);
    // void readEigenMatLDBlockBinFile(const string &geneEigenMatrixFile, const double eigenCutoff,const double diag_mod = 0.1);
    // read single gene eigen file
    void readEigenMatrixGene(const string &infoFile, const double eigenCutoff,const string LDMatType,const unsigned includeChr,const double diag_mod = 0.1);
    void readEigenMatGeneInfoFile(const string &infoFile);
    void readEigenMatGeneSnpInfoFile(const string &snpInfoFile,const unsigned includeChr);
    void readEigenMatBinFile(const string &geneEigenMatrixFile, const double eigenCutoff,const string LDMatType,const double diag_mod = 0.1);
    void readSBayesRCBlockLdmInfoFile(const string &infoFile);
    void readSBayesRCBlockLdmSnpInfoFile(const string &snpInfoFile);

    int readEigenMatGeneInfoFile(vector<GeneInfo *> &geneInfoVecLD, map<string, GeneInfo *> &geneInfoMapLD,const string &infoFile);
    void readEigenMatGeneSnpInfoFile(vector<EqtlInfo *> &eqtlInfoVecLD, map<string, EqtlInfo *> &eqtlInfoMapLD,const string &snpInfoFile);
    void readEigenMatGeneBinFile(vector<GeneInfo*> &geneInfoVecLD,map<string, GeneInfo *> &geneInfoMapLD,vector<EqtlInfo*> &eqtlInfoVecLD,map<string, EqtlInfo *> eqtlInfoMapLD, const string &geneEigenMatrixFile, const double eigenCutoff,const double diag_mod = 0.1);
    void alignxQTLGeneAndBESDGeneInfo(vector<GeneInfo *> &geneInfoVecBESD,map<string, GeneInfo *> &geneInfoMapBESD,vector<EqtlInfo *> &eqtlInfoVecBESD,map<string, EqtlInfo *> &eqtlInfoMapBESD, double diag_mod = 0.1);
    void alignxQTLGeneAndBESDGeneInfoWithoutImp(vector<GeneInfo *> &geneInfoVecBESD,map<string, GeneInfo *> &geneInfoMapBESD,vector<EqtlInfo *> &eqtlInfoVecBESD,map<string, EqtlInfo *> &eqtlInfoMapBESD, double diag_mod = 0.1);
    void imputexQTLEffect(GeneInfo * &gene,GeneInfo * &geneBESD,map<string, EqtlInfo *> &eqtlInfoMapBESD,vector<int> gene2CisSnpVecBESDLocal, double diag_mod);
    ///////////// Summary-level BayesOmics function ///////////
    /////////////////////////////////////////////////////////////////////////////////////
    ///////////             make low-rank matrix for gene region           //////////////
    /////////////////////////////////////////////////////////////////////////////////////
    ///////////// Step 2.1 read LD matrix of ld blocks info
    /// Read SMR BESD format
    void readEsiFile(const string &esiFile, bool &smrBesdBool,const bool haveSnpInfo);
    void readEpiFile(const string &epiFile);
    void readBesdFile(const string &besdFile,const bool &hasSnpInfo,const bool &makeGeneLDEigenBool, bool &smrBesdBool);
    void readQeuryGZFormat(const string &eqtlSummaryQueryFile,vector<GeneInfo *> &geneInfoVecLocal,map<string, GeneInfo *> &geneInfoMapLocal,
                            vector<EqtlInfo *> &eqtlInfoVecLocal,map<string, EqtlInfo *> &eqtlInfoMapLocal,
                            bool &smrBesdBool,const bool hasSnpInfo,const bool hasGeneLdmInfo);
    void readGeneSampleSizeFile(const string &geneSamSizeFile);
    void includeGene(const string &includeGeneFile);
    void excludeGene(const string &excludeGeneFile);
    void includeSpecificGene(const string &includeSpecificGeneID);
    // Do eigen decomposition
    // get eigen matrix from genotype directly
    void calcGeneLDEigenDecomBasedBesd(const string &eqtlSummaryFile, const string &eqtlSummaryQueryFile, const string &geneListFile,const string specificGeneID, const string &filename, const double cisRegionWind =2 ,const double eigenCutoff = 0.9995,bool debugBool = false);
    void outputEigenMatIndInfoFromLDMat(const string LDMatType, const string filename, const bool &mergeBool);
    void calcAndOutputEigenMatBinFromLDMat(const string bedFile, const vector<string> snplists,const string filename,const string LDMatType,const double eigenCutoff);

    void MergeMultiBinEigen(const string filename,const string LDMatType);

    void mergeMultiEigenMat(const string title, const string besdListFile,const string LDMatType, const double &eigenCutoff);
    bool readMultiEigenMatInfoFile(const string title,vector<GeneInfo *> &geneInfoVecLD, map<string, GeneInfo *> &geneInfoMapLD);
    bool readMultiEigenMatSnpInfoFile(const string title,vector<EqtlInfo *> &eqtlInfoVecLD, map<string, EqtlInfo *> &eqtlInfoMapLD);
    bool readMultiEigenMatBinFile(const string title, vector<GeneInfo*> &geneInfoVecLD,map<string, GeneInfo *> &geneInfoMapLD,vector<EqtlInfo*> &eqtlInfoVecLD,map<string, EqtlInfo *> eqtlInfoMapLD, float eigenCutoff);
    long estimateSamSize(VectorXd &beta,VectorXd &se);
    double adjeQTLSE(double beta,double p);

    /////////////////////////////////////////////
    ///////////      debug mode     /////////////
    /////////////////////////////////////////////
    void extractParameterFromIndModel(const string bedFile);
    void extractParameterFromSumModel(const string bedFile);
    void cleanUpUselessParameters();
};

#endif /* data_hpp */
