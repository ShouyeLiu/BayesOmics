#include "MemoryBudget.hpp"
// SPDX-License-Identifier: GPL-3.0-or-later
//
// This file is part of BayesOmics, a statistical genetics software package
// developed by Shouye Liu.
//
// Copyright (C) 2025 Shouye Liu <syliu.xue@foxmail.com>
//
// BayesOmics is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// BayesOmics is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with BayesOmics. If not, see <https://www.gnu.org/licenses/>.


#include "Data.hpp"


vector<EqtlInfo*> Data::makeIncdEqtlInfoVec(const vector<EqtlInfo*> &eqtlInfoVec){
    vector<EqtlInfo*> includedEqtls;
    cisSnpIDVec.clear();
    cisSnpID2IdxMap.clear();
    includedEqtls.clear();
    EqtlInfo * eqtl = NULL;
    map<string, SnpInfo*>::iterator iterSnp;
    for (unsigned i=0,j=0; i< numeQTLs; ++i) {
        eqtl = eqtlInfoVec[i];
        if(eqtl->included) { // 
            eqtl->index = j;  // reindex snps
            includedEqtls.push_back(eqtl);
           // snpEffectNames.push_back(snp->rsID);
            cisSnpIDVec.push_back(eqtl->rsID);
            cisSnpID2IdxMap.insert(pair<string, int>( eqtl->rsID ,j) );
            j++;
        }
    }
    if(includedEqtls.size() == 0){
        LOGGER.e(0,"Zero eQTLs are kept.");
    }
    return includedEqtls;
}

void Data::makeIncdEqtlInfoVecBasedOnGeneInfoVec(){
    EqtlInfo * eqtl = NULL;
    GeneInfo * gene = NULL;
    map<string, EqtlInfo*>::iterator iterEqtl;
    string eqtlID;
    // set all flag as false
    for (unsigned i=0; i<numIncdEqtls; ++i) {
        eqtl = incdEqtlInfoVec[i];
        eqtl->included = false;
        eqtl->isInGene = false;
    }
    //// now we filter xQTLs within the gene need to be removed
    for (unsigned i=0; i< numKeptGenes; ++i) {
        gene = keptGeneInfoVec[i];
        for(unsigned j = 0; j < gene->cisSnpNameVec.size(); j++ ){
            eqtlID = gene->cisSnpNameVec[j];
            if(gene->cisSnpID2IdxMapInGene.find(eqtlID) != gene->cisSnpID2IdxMapInGene.end()){
                // this eqtl is here
                eqtl = eqtlInfoMap.find(eqtlID)->second;
                eqtl->included = true;
                eqtl->isInGene = true;
            }
        } // loop xqtls for each gene
    } // loop all available genes
    ////////////////////////////
    // now we need to remove xQTLs within non-available genes
    incdEqtlInfoVec = makeIncdEqtlInfoVec(eqtlInfoVec);
    numIncdEqtls = incdEqtlInfoVec.size();
}

vector<GeneInfo*> Data::makeKeptGeneInfoVec(const vector<GeneInfo*> &geneInfoVec){
    vector<GeneInfo*> keptGenes; 
    geneEffectNames.clear();
    GeneInfo *gene = NULL;
    for (unsigned i=0, j=0; i< numGenes; ++i) {
        gene = geneInfoVec[i];
        // currently, we will remove genes whose xQTLs less thann one
        // if(gene->cisSnpNameVec.size() < 2) gene->kept = false;
        if(gene->kept) {
            gene->index = j++;  // reindex inds
            keptGenes.push_back(gene);
            geneEffectNames.push_back(gene->ensemblID);
        }
    }
    if(keptGenes.size() == 0){
        LOGGER.e(0,"Zero gene are kept.");
    }
    return keptGenes;
}



vector<IndInfo*> Data::makeKeptIndInfoGeneVec(const vector<IndInfo*> &indInfoGeneVec){
    vector<IndInfo*> keptInds;
    map<string, IndInfo*>::iterator it;
    // map<string, IndInfo*>::iterator iterIndGene,iterInd;
    keptInds.reserve(numIndGenes);
    // eqtlIndNameVec.clear();
    IndInfo *ind = NULL,*indGene = NULL;
    for (unsigned i=0, j=0; i<numIndGenes; ++i) {
        ind = indInfoGeneVec[i];
        if(ind->kept) {
            ind->index = j++;  // reindex inds
            keptInds.push_back(ind);
            // eqtlIndNameVec.push_back(ind->famID); // ind
        } else{
            it = indInfoMap.find(ind->catID);
            if(it != indInfoMap.end()){
                it->second->hasEQTL = false;
            }
            // ind->hasEQTL = false;
        }
    }
    if(keptInds.size() == 0){
        LOGGER.e(0,"Zero individuals in genes are kept.");
    }
    return keptInds;
}

void Data::includeMatchedEqtl(){
    incdEqtlInfoVec = makeIncdEqtlInfoVec(eqtlInfoVec);
    numIncdEqtls = (unsigned) incdEqtlInfoVec.size();
    // chromosome distribuion of SNPs from cis-region
    map<int, vector<EqtlInfo*> > chrmap;
    map<int, vector<EqtlInfo*> >::iterator it;
    for (unsigned i=0; i<numIncdEqtls; ++i) {
        EqtlInfo *eqtl = incdEqtlInfoVec[i];
        if (chrmap.find(eqtl->chrom) == chrmap.end()) {
            chrmap[eqtl->chrom] = *new vector<EqtlInfo*>;
        }
        chrmap[eqtl->chrom].push_back(eqtl);
    }
    // numChroms = (unsigned) chrmap.size();
    if(chromInfoVec.size() == 0){
        for (it=chrmap.begin(); it!=chrmap.end(); ++it) {
            int id = it->first;
            vector<EqtlInfo*> &vec = it->second;
            ChromInfo *chr = new ChromInfo(id, (unsigned)vec.size(), vec[0]->index, vec.back()->index,vec[0]->physPos,vec.back()->physPos);
            chromInfoVec.push_back(chr);
            chr->startSnpID = vec[0]->rsID;
            chr->endSnpID = vec.back()->rsID; 
            chromInfoMap.insert(pair<int, ChromInfo*>(id, chr));
        }
    }
    if(numIncdEqtls == 0) { LOGGER.e(0, to_string(numIncdEqtls) + " SNPs from cis-region on " + std::to_string(chrmap.size()) + " chromosomes are included.");}
    // LOGGER << numIncdEqtls << " SNPs from cis-region on " << (unsigned) chrmap.size() << " chromosomes are included." << endl;
    LOGGER << numIncdEqtls << " SNPs from cis-region on " << (unsigned) chrmap.size() << " chromosomes (ID: ";
    for(unsigned k = 0; k < chromInfoVec.size(); k ++){
        LOGGER << chromInfoVec[k]->id << " ";
    }
    LOGGER << ") are included." << endl;
}

////////////////////////////////////////////////////////////
//////////// Individual level BayesOmics
////////////////////////////////////////////////////////////
/// Add five functions to build BayesOmics model.
// Step 1. readGeneInfoFile to read geneInfo
// these function will be used in GCTB::inputIndInfo();
// Build GeneInfo class

void Data::initeQTLVariances(const double heritability, const double propVarRandom){
    varGenotypiceQTL.resize(numKeptGenes);
    varResidualeQTL.resize(numKeptGenes);
    varRandomeQTL.resize(numKeptGenes);
    for(unsigned j = 0; j < numKeptGenes; ++j){
        // double varPhenotypic = ypyeQTL[j]/numKeptIndseQTL;
        double varPhenotypic = varPhenotypiceQTL[j];
        varGenotypiceQTL[j] = varPhenotypic * heritability;
        varResidualeQTL[j]  = varPhenotypic - varGenotypiceQTL[j];
        varRandomeQTL[j]    = varPhenotypic * propVarRandom;
    }

}

void Data::readGeneInfoFile(const string geneInfoFile,const string subGenePath){
    if (geneInfoFile.empty()) return;
    ifstream in(geneInfoFile.c_str());
    if (!in) LOGGER.e(0, " can not open the gene info file [" + geneInfoFile + "] to read.");
    LOGGER << "Reading gene info file from [" + geneInfoFile + "]." << endl;

    geneInfoVec.clear();
    geneInfoMap.clear();
    indInfoGeneVec.clear();
    indInfoGeneMap.clear();

    unsigned chr;
    int start, end;
    string probe, geneID, geneOri,singleGenePheFile;
    unsigned idx = 0;

    Gadget::Tokenizer colData;
    string header;
    string sep(" \t\r");
    getline(in, header);
    // Step 1. read plist file, where information per gene was stored ( chr, gene start, gene end, ensgid ID, genePath)
    // genePath is absolute path;
    // So when genotypes from gene cis-regions are needed, we should map genotypes to gene based on gene cis-windows. 
    while (in >> chr >> start >> end >> geneID >> singleGenePheFile) {
        GeneInfo *gene = new GeneInfo(idx++, geneID, chr);
        gene->start = start;
        gene->end  = end;
        gene->midPhyPos = (double) (start + end)/ (double) 2.0;  // here we ued middle position as cis-region
        geneInfoVec.push_back(gene);
        if (geneInfoMap.insert(pair<string, GeneInfo*>(geneID, gene)).second == false) {
            LOGGER.e(0, " Duplicate gene ID found: \"" + geneID + "\".");
        }
        // Step 2. read phenotype information from each gene. Here based on genePath from plist file, we located gene file first,
        // and then extract relevant genes. Please note, individuals from each gene should be a subset of individuals from plink fam file,
        // and don't need to be consistent with individuals from complex trait, considering the possibility that gene sample size may be 
        // not equal to sample size of complex trait.
        readPhenoFromGeneFile(gene->ensemblID, singleGenePheFile,subGenePath);
    }
    // gene info summary
    numGenes = (unsigned) geneInfoVec.size();
    LOGGER << numGenes << " genes to be included from [" + geneInfoFile + "]." << endl;

    // individuals with gene info
    numIndGenes = (unsigned) indInfoGeneVec.size();
    LOGGER << numIndGenes << " individuals have at least one genes." << endl;
}

///////////////////////////////////////////////////////////////////////////////////
//// this function aims to per protein values and will be used in readGeneInfoFile
////  build IndInfo for each gene and select individual with gene info
void Data::readPhenoFromGeneFile(const string geneID, const string & singleGenePheFile, const string subGenePath){
    string filename;
    if (singleGenePheFile.empty()) return;
    if(!subGenePath.empty()){
        filename = subGenePath + "/" + singleGenePheFile;
    } else {
        filename = singleGenePheFile;
    }
    ifstream in(filename.c_str());
    if (!in) LOGGER.e(0, " can not open the file [" + filename + "] to read.");
    LOGGER << "Reading molecular phenotype of "+ geneID +  " file from [" + filename + "]." << endl;
    
    map<string, IndInfo*>::iterator it,iterIndGene;
    map<string, IndInfo*> indInfoGeneLocalMap;
    IndInfo *ind = NULL;
    Gadget::Tokenizer colData;
    string inputStr;
    string sep(" \t\r");
    string id;
    unsigned line=0;
    while (getline(in,inputStr)) {
        colData.getTokens(inputStr, sep);
        id = colData[0] + ":" + colData[1];
        it = indInfoMap.find(id);  // can be found in fam file
        if (it != indInfoMap.end() && colData[2] != "NA") {
            if(line%4096==0)MemoryBudget::require(4096ULL*1152,"molecular phenotype records and ID mappings");
            ind = it->second;
            ind->hasEQTL = true;
            // judge if this individual is in indGene or not;
            iterIndGene = indInfoGeneMap.find(id);
            if(iterIndGene != indInfoGeneMap.end()){
                // the individual already in indInfoGeneVec class
                iterIndGene->second-> genePheMap.insert(pair<string,double>(geneID, boost::lexical_cast<double>(colData[2].c_str())  ));
                // cout << "geneInd: " << iterIndGene->second->indID << " vlaue: " << boost::lexical_cast<double>(colData[2].c_str())  << endl;
            } else {
                // the individual is new here, then build new indGene info
                IndInfo *indGene = new IndInfo( ind->index,
                    ind->famID, 
                    ind->indID,
                    ind->fatherID,
                    ind->motherID, 
                    ind->sex);
                    indGene->hasEQTL = true;
                indGene->genePheMap.insert(pair<string,double>(geneID, boost::lexical_cast<double>(colData[2].c_str())));
                indInfoGeneVec.push_back(indGene);
                indInfoGeneMap.insert(pair<string, IndInfo*>(indGene->catID, indGene));
            }
            if (indInfoGeneLocalMap.insert(pair<string, IndInfo*>(ind->catID, ind)).second == false) {
                LOGGER.e(0, " Duplicate individual ID found: \"" + ind->catID + " in file " + filename);
            }
            ++line;
        } else {
            continue;
        }
    }
    in.close();
}

void Data::keepMatchedIndGene(const string &keepIndFile, const unsigned keepIndMax){  // keepIndFile is optional
    map<string, IndInfo*>::iterator it;
    IndInfo *ind = NULL;
    set<string> keep;

    if (!keepIndFile.empty()) {
        ifstream in(keepIndFile.c_str());
        if (!in) LOGGER.e(0,"Can not open the file [" + keepIndFile + "] to read.");
        LOGGER.i(0, "Reading kept individuals with molecular phenotypes from [" + keepIndFile + "].");
        string fid, pid;
        keep.clear();
        while (in >> fid >> pid) {
            keep.insert(fid + ":" + pid);
        }
        in.close();

        numKeptIndsGene = 0;

        for (unsigned i=0; i<indInfoGeneVec.size(); ++i) {
            ind = indInfoGeneVec[i];
            if (keep.find(ind->catID) == keep.end()) ind->kept = false;
        }
    }
    
    
    keptindInfoGeneVec = makeKeptIndInfoGeneVec(indInfoGeneVec);
    numKeptIndsGene =  (unsigned) keptindInfoGeneVec.size();

    if (numKeptIndsGene == 0) LOGGER.e(0,"No individual with gene is retained for analysis, please check individual difference between phenotype and keep file" + keepIndFile);
    LOGGER << numKeptIndsGene << " matched individuals with at least one gene are kept in the molecular phenotype." << endl;
}


///// Step 2. function mapGwasSnpToGeneCisRegion aims to construct various maps
///// this function will be used in GCTB::inputSnpInfo();
/// Build EqtlInfo class 
void Data::mapGwasSnpToGeneCisRegion(const string &bedFile, const bool noscale, const double cisRegionWind){
    ///////////////////////////////////////////////////////
    /////// Step 1. construct snp name map
    ///////////////////////////////////////////////////////
    int i,j;
    vector<locus_bp> snpVec;
    SnpInfo *snp = NULL;
    map<int, string>  chrEndSnp;
    map<string, int> snpNameMap;
    for (i = 1; i < numIncdSnps; i++) {
        snp = incdSnpInfoVec[i];
        if(incdSnpInfoVec[i]->chrom != incdSnpInfoVec[i-1]->chrom){
            chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[i - 1]->chrom,incdSnpInfoVec[i - 1]->rsID ));
        }
    }
    chrEndSnp.insert(pair<int, string>(incdSnpInfoVec[numIncdSnps - 1]->chrom,incdSnpInfoVec[numIncdSnps - 1]->rsID ));
    for (i = 0; i < numIncdSnps; i++) {
        snp = incdSnpInfoVec[i];
        snpVec.push_back(locus_bp(snp->rsID, snp->chrom, snp->physPos ));
        snpNameMap.insert(pair<string,int>(snp->rsID, i));
    }
    if (!geneSnpMapFile.empty()) {
        ifstream mapIn(geneSnpMapFile);
        if (!mapIn) throw string("Cannot open gene-SNP map: ")+geneSnpMapFile;
        string header, geneID, snpID;
        getline(mapIn,header);
        istringstream columns(header); string a,b,extra,line;
        if (!(columns >> a >> b) || a!="Gene" || b!="SNP" || (columns >> extra))
            throw string("Gene-SNP mapping header must be: Gene SNP");
        set<pair<string,string>> seen;
        for (auto *g : geneInfoVec) { g->cisSnpNameVec.clear(); g->hasEqtl=false; }
        while (getline(mapIn,line)) {
            if (line.find_first_not_of(" \t\r")==string::npos) continue;
            istringstream row(line);
            if (!(row >> geneID >> snpID) || (row >> extra))
                throw string("Malformed gene-SNP mapping row: ")+line;
            if (!geneInfoMap.count(geneID) || !snpInfoMap.count(snpID) || !snpInfoMap.at(snpID)->included)
                throw string("Unknown or excluded gene/SNP in explicit map: ")+geneID+" "+snpID;
            if (!seen.emplace(geneID,snpID).second) throw string("Duplicate gene-SNP mapping: ")+geneID+" "+snpID;
            geneInfoMap.at(geneID)->cisSnpNameVec.push_back(snpID);
            geneInfoMap.at(geneID)->hasEqtl=true;
            snpInfoMap.at(snpID)->iseQTL=true;
        }
        if (!mapIn.eof()) throw string("Malformed gene-SNP mapping: ")+geneSnpMapFile;
        for (auto *g : geneInfoVec) g->kept=g->hasEqtl;
        LOGGER << "Read " << seen.size() << " explicit gene-SNP pairs." << endl;
    } else {
    /////////////////////////////////////////
    // Step 2. Map snps to genes
    /////////////////////////////////////////
    LOGGER << "Mapping the physical positions of genes to SNP data (gene boundaries: " << cisRegionWind / 1000 << "Kb away from the middle position of the gene) ..." << endl;
    vector<string> gene2snp_1(numGenes), gene2snp_2(numGenes);
    vector<locus_bp>::iterator iter;
    map<int, string>::iterator chrIter;
    GeneInfo *gene;

    #pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numGenes; i++) {
        // find lowest snp_name in the gene
        gene = geneInfoVec[i];
        // iter = find_if(snpVec.begin(), snpVec.end(), locus_bp( geneInfoVec[i]->geneID ,geneInfoVec[i]->chrom,geneInfoVec[i]->start - cisRegionWind));
        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp( gene->ensemblID ,gene->chrom, gene->midPhyPos - cisRegionWind));
        if (iter != snpVec.end()) gene2snp_1[i] = iter->locusName;
        else gene2snp_1[i] = "NA";
    }
    #pragma omp parallel for private(iter, chrIter)
    for (i = 0; i < numGenes; i++) { 
        gene = geneInfoVec[i];               
        if (gene2snp_1[i] == "NA") {
            gene2snp_2[i] = "NA";
            continue;
        }
        iter = find_if(snpVec.begin(), snpVec.end(), locus_bp(gene->ensemblID, gene->chrom, gene->midPhyPos + cisRegionWind));
        if (iter != snpVec.end()){
            if (iter->bp ==  gene->midPhyPos + cisRegionWind){
                gene2snp_2[i] = iter->locusName;
            }else {
                if(iter!=snpVec.begin()){
                    iter--;
                    gene2snp_2[i] = iter->locusName;
                }
                else gene2snp_2[i] = "NA";
            }
        }
        else {
            chrIter = chrEndSnp.find(gene->chrom);
            if (chrIter == chrEndSnp.end()) gene2snp_2[i] = "NA";
            else gene2snp_2[i] = chrIter->second;
        }
    }
    int mapped = 0;
    for (i = 0; i < numGenes; i++) {
        gene = geneInfoVec[i];
        if (gene2snp_1[i] != "NA" && gene2snp_2[i] != "NA") 
        {
            mapped++;
            // gene->kept = true;
        }else {
            gene->kept = false;
        }
    }
    if (mapped < 1) LOGGER.e(0, " No gene can be mapped to the SNP data. Please check the input data regarding chromosome and bp.");
    else LOGGER << mapped << " genes have been mapped to SNP data." << endl;
    //////////////////////////////////////////////////////////
    map<string, SnpInfo*>::iterator iterSnp;
    map<string, int>::iterator iter1, iter2;
    // geneEffectNames.clear();
    vector<int> snpNumInGene(numGenes);
    for (i = 0; i < numGenes; i++) {
        gene = geneInfoVec[i]; 
        gene->windows = cisRegionWind;
        iter1 = snpNameMap.find(gene2snp_1[i]);
        iter2 = snpNameMap.find(gene2snp_2[i]);
        if (iter1 == snpNameMap.end() || iter2 == snpNameMap.end() || iter1->second >= iter2->second) gene->kept = false;
        snpNumInGene[i] = iter2->second - iter1->second + 1;
        if(!gene->kept) {
            LOGGER.w(0,"Gene [" + gene->ensemblID + "] is removed due to no moQTLs there.");
            continue;
        }
        gene->hasEqtl = true;
        for (j = iter1->second; j <=  iter2->second; j++) {
            // collect eqtls into gene
            gene->cisSnpNameVec.push_back(incdSnpInfoVec[j]->rsID);
            // here we label the gwas snps as an eqtl, which
            // will be used when constructing EqtlInfo class later.
            iterSnp = snpInfoMap.find(incdSnpInfoVec[j]->rsID);
            iterSnp->second->iseQTL = true;
            // cout << incdSnpInfoVec[j]->rsID <<  " ";
        } 
    }
    }
    ///////////////////////////////////////////////////////
    /////// Step 5. Construct EqtlInfo class
    ///////////////////////////////////////////////////////
    eqtlInfoVec.clear();
    eqtlInfoMap.clear();
    unsigned idx = 0;
    for (i = 0; i < numIncdSnps ; i++) {
        snp = incdSnpInfoVec[i];
        if(!snp->iseQTL) {continue; }
        EqtlInfo *eqtl = new EqtlInfo(idx++,
        snp->rsID,
        snp->a1,
        snp->a2,
        snp->chrom,
        snp->genPos,
        snp->physPos,
        snp->af);
        eqtlInfoVec.push_back(eqtl);
        if (eqtlInfoMap.insert(pair<string, EqtlInfo*>(snp->rsID, eqtl)).second == false) {
            LOGGER.e(0, " Duplicate eQTL SNP ID found: \"" + snp->rsID + "\".");
        }
    }
    numeQTLs = eqtlInfoVec.size();
    ///////////////////////////////////////////////////////
    /////// Step 4. re-construct geneInfo class
    ///////////////////////////////////////////////////////
    // this step is to remove gene without eqtls
    keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    numKeptGenes = (unsigned) keptGeneInfoVec.size();
    // remove eqtls
    incdEqtlInfoVec = makeIncdEqtlInfoVec(eqtlInfoVec);
    numIncdEqtls = (unsigned) incdEqtlInfoVec.size();

    keptindInfoGeneVec = makeKeptIndInfoGeneVec(indInfoGeneVec);
    numKeptIndsGene =  (unsigned) keptindInfoGeneVec.size();
    readBedFileForGene(bedFile,noscale);
    includeMatchedEqtl();
    keptGeneInfoVec  = makeKeptGeneInfoVec(geneInfoVec);
    numKeptGenes = (unsigned) keptGeneInfoVec.size();
    LOGGER << "" << numIncdEqtls << " SNPs are included after mapping snps to "
            << numKeptGenes << " genes." << endl;
    // construct various maps
    ConstructGwasEqtlGeneMaps();
}

///// Step 2. function  aims to transfer various parameter to BayesOmics function
///// this function will be used in GCTB::buildModel();
///// Please note only one genotype matrix from cis-region is constructed and gene-specific
///// cis-region matrix will be constructed via genePheIdxMap due to unequal smaple size between 
///// complex trait and molecular trait
void Data::readBedFileForGene(const string &bedFile,const bool noscale) {
    if(!numeQTLs || !numKeptIndsGene)throw std::runtime_error("No molecular SNP/individual retained for BED analysis");
    const bool streaming=boundedIndividual && MemoryBudget::streamGenotypes(numKeptIndsGene,numeQTLs);
    if(!streaming)MemoryBudget::require(MemoryBudget::bytes(numKeptIndsGene,numeQTLs),"dense molecular genotypes");
    readBedColumns(noscale,bedFile,true,!streaming);
}

// this function will be shared in both individual and summary level model
void Data::ConstructGwasEqtlGeneMaps(){
    //////////////////////////////////////////////////////////
    /////// Step 3. Construct various maps
    //////////////////////////////////////////////////////////
    //// Our final aim is to construct the following variables
    ///////////////////////////////////////////////////////////
    gwasAndGeneEffectNames.clear();
    gwasAndGeneEffectNames.push_back("nonEqtl");

    gwasSnpID2geneIDMap.clear();
    gwasSnpID2geneIdxMap.clear();
    gene2gwasSnpMap.clear();
    gene2cisSnpMap.clear();
    gwas2SnpIDMap.clear();
    gene2cisSnpIDMap.clear();
    geneEffectNames.clear();
    geneID2IdxMap.clear();

    numNonEqtl = numIncdSnps;
    numEqtlOverlap = 0;
    numEqtl = 0;

    set<string> eQTLUniqSet; // used to calculate num of eqtls in genic region.
    eQTLUniqSet.clear();
    map<string, SnpInfo*>::iterator iterSnp;
    map<string, int>::iterator iter1, iter2;

    vector <double>  snp2pqeQTLPerGene;
    vector<int> gene2cisSnpIdxPerGene;
    vector<string> cisSnpNameVecPerGene;
    int numEqtlPerGene;
    // for intergenic region
    vector<string> snpInterGenic;
    for(unsigned j = 0; j < numIncdSnps;j++){
        SnpInfo *snp = incdSnpInfoVec[j];
        if(snp->iseQTL) continue;
        snpInterGenic.push_back(snp->rsID);
    }
    gwas2SnpIDMap.insert(pair<int,vector<string>> (0,snpInterGenic));

    snp2pqeQTL.clear();
    auto startTime = std::chrono::steady_clock::now();
    for(unsigned j = 0, geneIdx = 0; j < numKeptGenes; j++){
        Gadget::showProgressBar(j,numKeptGenes , startTime,"Map gene info to GWAS snps");
        snp2pqeQTLPerGene.clear();
        GeneInfo * gene = keptGeneInfoVec[j];
        cisSnpNameVecPerGene = gene->cisSnpNameVec;
        numEqtlPerGene = cisSnpNameVecPerGene.size();
        gene->cisSnpNameVec.clear();
        gene->gene2CisSnpVec.clear();
        gene->gene2GwasSnpVec.clear();
        for (unsigned i = 0; i < numEqtlPerGene; ++i){
            EqtlInfo * eqtl = eqtlInfoMap[cisSnpNameVecPerGene[i]];
            if(!eqtl->included) continue;
            snp2pqeQTLPerGene.push_back(eqtl->twopq);
            // link gwas index and eqtl index
            iterSnp = snpInfoMap.find(eqtl->rsID);
            // iterSnp->second->iseQTL = true;
            gene->cisSnpNameVec.push_back(eqtl->rsID);
            // find cis index location in each genes
            gene->gene2CisSnpVec.push_back(cisSnpID2IdxMap.at(eqtl->rsID));
            // find relative cis-index in gwas snplist
            gene->gene2GwasSnpVec.push_back(iterSnp->second->index);
            eQTLUniqSet.insert(eqtl->rsID);
            numEqtlOverlap ++;
            eqtl->eqtl2gene_count++; // for each eqtl, how many genes;
            /// map snp to gene 
            if(gwasSnpID2geneIDMap.find(iterSnp->first) != gwasSnpID2geneIDMap.end()){
            //// already in gene 
                gwasSnpID2geneIDMap[iterSnp->first].push_back(gene->ensemblID);
                gwasSnpID2geneIdxMap[iterSnp->first].push_back(geneIdx);
            } else{
            // create new 
                gwasSnpID2geneIDMap.insert({iterSnp->first,{gene->ensemblID} });
                std::vector<int> geneVec(1, geneIdx);
                gwasSnpID2geneIdxMap.insert({iterSnp->first,geneVec} );
            }
        }
        if(gene->cisSnpNameVec.size() == 0) {
            gene->kept = false;
            LOGGER.w(0,"Gene [" + gene->ensemblID + "] is removed due to no moQTLs there.");
            continue;
        }
        snp2pqeQTL.push_back(Eigen::Map<Eigen::VectorXd>(snp2pqeQTLPerGene.data(), snp2pqeQTLPerGene.size()));
        gene2cisSnpMap.insert(pair<int, vector<int>>(geneIdx,gene->gene2CisSnpVec));
        gene2cisSnpIDMap.insert(pair<int,vector<string>>(geneIdx,gene->cisSnpNameVec));
        gene2gwasSnpMap.insert(pair<int,vector<int>> (geneIdx,gene->gene2GwasSnpVec));
        gwasAndGeneEffectNames.push_back(gene->ensemblID);
        gwas2SnpIDMap.insert(pair<int,vector<string>> (geneIdx +1, gene->cisSnpNameVec));
        ///// map gene to snp
        geneEffectNames.push_back(gene->ensemblID);
        geneID2IdxMap.insert(pair<string,int>(gene->ensemblID,geneIdx));
        geneIdx ++;
    }
    numEqtl = eQTLUniqSet.size();
    numNonEqtl = numIncdSnps - numEqtl;
}

// this function will be used after readGeneInfoFile function in GCTB::inputIndInfo() function to solve match issue.
void Data::ConstructGenePheAndwAcorr(){  
    map<string, IndInfo*>::iterator iterIndGene,iterInd;
    std::map<string,double>::iterator iterGenePhe;
    IndInfo *indGene, *ind;
    GeneInfo * gene;

    /////// now we loop genes to store the data
    unsigned indIdx = 0;
    VectorXd nTmp;  
    genePheVec.clear();
    wAcorr.resize(numKeptGenes);
    vector <double>  genePhePerGene;
    vector<int> genePheIdxPerGene;
    ypyeQTL.resize(numKeptGenes);
    varPhenotypiceQTL.resize(numKeptGenes);

    vector<string> eqtlIndNamePerGene;  
    // LOGGER << "num gene: " << numKeptGenes << endl;
    for (unsigned i =0; i < numKeptGenes; i++){
        genePhePerGene.clear();
        genePheIdxPerGene.clear(); // map phen to fam file
        eqtlIndNamePerGene.clear();
        gene = keptGeneInfoVec[i];
        nTmp.resize(gene->cisSnpNameVec.size());
        // here we loop individuals just want to keep consistent order between gwas y and gene phenotype based on the order of fam file
        for (unsigned j = 0,indIdx = 0; j < numInds; ++ j ) {
            ind = indInfoVec[j];
            if(!ind->hasEQTL) continue;
            // find ind with gene info
            iterIndGene = indInfoGeneMap.find(ind->catID);
            indGene = iterIndGene->second;
            if(!indGene->kept) continue;
            
            // extrect gene expression value from indGene info based on given gene
            iterGenePhe = indGene->genePheMap.find(gene->ensemblID);
            if(iterGenePhe != indGene->genePheMap.end()){
                genePhePerGene.push_back(iterGenePhe->second);
                // when constructing index, we need to map to raw genotype matrix X, so that we 
                // just need to share data
                genePheIdxPerGene.push_back(indIdx); // will be used to select individuals to 
                eqtlIndNamePerGene.push_back(indGene->famID + "_" + indGene->indID );
            }
            // since the dimension of genotype from cis-region may be different from that of whole genome,
            // we need genePheIdxPerGene to store indIdx for each gene
            indIdx ++;
        }
        /// now store gene value into vectorDat
        VectorXd geneValEigen = VectorXd::Map(genePhePerGene.data(), genePhePerGene.size());
        nTmp.setConstant(genePhePerGene.size());
        VectorDat vectorDat = VectorDat(gene->cisSnpNameVec,nTmp);
        // ypyeQTL[i] = (geneValEigen.array()-geneValEigen.mean()).square().sum();
        // varPhenotypiceQTL[i] = ypyeQTL[i]/((double) geneValEigen.size() -1);
        varPhenotypiceQTL[i] = (geneValEigen.array()-geneValEigen.mean()).square().sum()/(geneValEigen.size()-1.0);
        // cout << "gene: " << gene->ensemblID << " i: " << i << " varPHe: " << varPhenotypiceQTL[i]  << endl;
        genePheVec.push_back(geneValEigen);
        neQTLVec.push_back(vectorDat);
        wAcorr[i] = geneValEigen.array() - Gadget::calcMean(geneValEigen);
        /// Store index
        genePheIdxMap.insert(pair<string,vector<int>> (gene->ensemblID,genePheIdxPerGene));
    }
}

void Data::buildBayesOmicsMME(const string bedFile, const bool noscale, const bool haveGene){
    ////////////////////////////////////////////////////////////
    // Construct gwas information
    ////////////////////////////////////////////////////////////
    n.resize(numIncdSnps);
    n.setConstant(genotype().rows());
    ZDat.emplace_back(snpEffectNames,genotype());
    wbcorr.array() = y.array() - Gadget::calcMean(y);
    // check map
    for (unsigned j =0; j < gene2gwasSnpMap.size(); j ++){
        for (unsigned i = 0; i < gene2gwasSnpMap[j].size(); i ++) {
            if (snpEffectNames[gene2gwasSnpMap[j][i]] != cisSnpIDVec[gene2cisSnpMap[j][i]]){
                cout << snpEffectNames[gene2gwasSnpMap[j][i]] << "::" << cisSnpIDVec[gene2cisSnpMap[j][i]] << endl;
                LOGGER.e(0, " gwas snp doesn't match eqtl snp.");
            }
        }
    }
    ////////////////////////////////////////////////////////////
    // Construct gene part information
    ////////////////////////////////////////////////////////////
    if(haveGene){
        LOGGER <<  "BayesOmics-Bivariate mode is chosen." << endl;
    }else {
        numNonEqtl = numIncdSnps;
        numEqtlOverlap = 0;
        numEqtl = 0;
        numKeptGenes = 0;
        gwasAndGeneEffectNames.clear();
        gwasAndGeneEffectNames.push_back("nonEqtl");
        gwasSnpID2geneIDMap.clear();
        gene2gwasSnpMap.clear();
        geneID2IdxMap.clear();
        snp2pqNonEqtl = snp2pq;
        LOGGER <<  "BayesOmics-Single mode is chosen." << endl;
    }
    ////////////////////////////////////////////////////////////
    ///////////////  Data summary
    ////////////////////////////////////////////////////////////
    // LOGGER << "Data summary: " << endl;
    // LOGGER << "Total GWAS SNPs: " << numIncdSnps << endl 
    //     << "Total genes: " << numKeptGenes << endl 
    //     << "Total cis-SNPs in gene cis-region: "        << numEqtl << endl
    //     // << "Total cis-SNPs in gene overlap resion (including overlapping cis-SNPs): "   <<  numEqtlOverlap  << endl
    //     << "Total GWAS SNPs in non-gene cis-region: " << numNonEqtl  << endl;
    // LOGGER << "\nData summary:" << endl;
    // LOGGER << boost::format("%60s %8s %8s\n") %"" %"mean" %"sd";
    // LOGGER << boost::format("%60s %8.3f %8.3f\n") %"GWAS Phenotypic variance" %Gadget::calcMean(varySnp) %sqrt(Gadget::calcVariance(varySnp));
    // LOGGER << boost::format("%60s %8.3f %8.3f\n") %"GWAS SNP heterozygosity" %Gadget::calcMean(snp2pq) %sqrt(Gadget::calcVariance(snp2pq));
    // LOGGER << boost::format("%60s %8.0f %8.0f\n") %"GWAS SNP sample size" %Gadget::calcMean(n) %sqrt(Gadget::calcVariance(n));
    // // xQTL info
    // LOGGER << boost::format("%60s %8.3f %8.3f\n") %"xQTL Phenotypic variance" %Gadget::calcMean(varPhenotypiceQTL) %sqrt(Gadget::calcVariance(varPhenotypiceQTL));
    // LOGGER << boost::format("%60s %8.3f %8.3f\n") %"xQTL SNP heterozygosity across genes " %Gadget::calcMean(snp2pqeQTL) %sqrt(Gadget::calcVariance(snp2pqeQTL));
    // LOGGER << boost::format("%60s %8.0f %8.0f\n") %"xQTL SNP average sample size across genes" %Gadget::calcMean(sampleSizeAcrossGene) %sqrt(Gadget::calcVariance(sampleSizeAcrossGene));
    // LOGGER << boost::format("%60s %8.3f %8.3f\n") %"xQTL SNP size in non-mQTLs region" % numNonEqtl % 0;
    // LOGGER << boost::format("%60s %8.3f %8.3f\n") %"xQTL SNP size in mQTLs region" % numEqtl % 0;
    // LOGGER << boost::format("%60s %8.3f %8.3f\n") %"xQTL average SNPs per molecular" % (double)(numEqtl/(double)numKeptGenes) % 0;
}


