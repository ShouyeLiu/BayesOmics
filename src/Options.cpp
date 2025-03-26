

// #include <cstdlib>
// #include <cmath>
#include <vector>
#include <string>
#include <sstream>
#include <set>
#include <algorithm>
#include <utility>
#include <numeric>

#include "Options.hpp"


  using std::vector;
  using std::string;
  using std::cerr;
  using std::endl;
//   namespace ublas = boost::numeric::ublas;
// Options::~Options() = default;
void Options::customCheck()
{
	// std::LOGGER<< "Generally, the number of CPU and GPU in the system will be checked.";
}
void Options::printCommandLine(int argc, const char * const argv[]){
  LOGGER << "Options:" << endl;
    // LOGGER << argv[0];
  for (int i = 1; i < argc; i++) {
    if (strlen(argv[i]) >= 2 && argv[i][0] == '-' && argv[i][1] == '-')
      LOGGER << "\n";
      bool hasSpace = false;
      for (int j = 0; j < (int) strlen(argv[i]); j++)
        if (isspace(argv[i][j]))
	    hasSpace = true;
      if (hasSpace) {
        if (argv[i][0] == '-') {
	      bool foundEquals = false;
	      for (int j = 0; j < (int) strlen(argv[i]); j++) {
	        LOGGER <<  argv[i][j];
	        if (argv[i][j] == '=' && !foundEquals) {
	          LOGGER << "\"";
	          foundEquals = true;
	        }
	      }
	      LOGGER << "\" ";
        }
      else
	      LOGGER << " " << argv[i];
      }
    else
      LOGGER << " " << argv[i];
  }
  LOGGER << endl << endl;
}

void Options::readProgramOptions(int argc, const char *const argv[])
// #ifndef UNIT_TEST   // #ifndef	means: if not defined
// 	try
// #endif
{
	boostProgramOptionsRoutine(argc, argv);
  customCheck();
}

void Options::setThread(void){
    omp_set_num_threads(numThreads);
    if (numThreads == 1) return;
    // if (multiThreadEigen) {
    //     Eigen::initParallel();
    //     Eigen::setNbThreads(numThread);
    //     cout << "Eigen library is using " << Eigen::nbThreads( ) << " threads." << endl;
    // }
#pragma omp parallel
    printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
}

  // populates members; error-checks
  bool Options::boostProgramOptionsRoutine(int argc, const char *const argv[])
  {
    vector <string> bimFileTemplates, bedFileTemplates;
    vector <string> removeFileTemplates, excludeFileTemplates, modelSnpsFileTemplates;
    vector <string> covarColTemplates, qCovarColTemplates;
    vector <string> dosageFileTemplates, bgenFileTemplates;
    string dosage2FileList, impute2FileList;
    string remlGuessStr;
    string sampleFile1, bgenSampleFileList;

    namespace po = boost::program_options;
    po::options_description typical("Basic options:");
    typical.add_options()
      ("help,h", "print help message with  options")
      // ("helpFull", "print help message with full option list")
      ("version,v", "show program version information")
      ("seed", po::value<int>(&seed), "random seed")
      ("thread", po::value<int>(&numThreads),"number of computational threads")
       ("out", po::value<string>(), "Specify output root filename")
       ("slurm-array-num", po::value<int>(&slurmArrayLimit), "set slurm array limit ");

      po::options_description indLevel("Individual-level model options:");
      indLevel.add_options()
      // Individual-level model
      // gwas genotype data parameters
      ("bfile", po::value<string>(&bedFile), "prefix of PLINK .fam, .bim, .bed files (see PLINK user manual for details).")
      ("chr", po::value<unsigned>(&includeChr), "Include SNPs on a specific chromosome in the analysis, e.g. chromosome 1.")
      ("extract", po::value<string>(&includeSnpFile), "Specify a list of SNPs to be included in the analysis.")
      ("exclude", po::value<string>(&excludeSnpFile), "Specify a list of SNPs to be excluded from the analysis.")
  
      // phenotype and covariate data parameters
      ("pheno", po::value<string>(&phenotypeFile),
       "GWAS trait phenotype file (header required; FID IID must be first two columns)")
      ("keep", po::value<string>(&keepIndFile), "Specify a list of individuals to be included in the GWAS analysis.")
      ("keep-gene-ind", po::value<string>(&keepIndGeneFile), "Specify a list of individuals with at least one gene to be included in the analysis.")
      ("mphone", po::value<unsigned>(&mphen), "GWAS trait phenotype column header")
      ("covar", po::value<string>(&covariateFile),
       "Input quantitative covariates from a plain text file, e.g. test.qcovar."
       "Each quantitative covariate is recognized as a continuous variable.")
      ("gene", po::value<string>(&geneInfoFile),"molecuar phenotype file")
      ("flist-path", po::value<string>(&subGenePath),"Specify path for esd  in filst ")
      ("plist-path", po::value<string>(&subGenePath),"Specify path for gene pheno in plist ");

      po::options_description sumLevel("Summary-level model options:");
      sumLevel.add_options()
      // Summary-level model
      // LD marix management
      ("make-full-ldm", " full LD matrix")
      ("make-band-ldm", " banded LD matrix")
      ("make-shrunk-ldm", " shrunk LD matrix")
      ("make-sparse-ldm", " sparse LD matrix")
      ("mldm", po::value<string>(&ldmatrixFile),"multiple LD matrix file")
      // way 1. generate blockwise low-rank LD matrix
      ("make-block-ldm", " blocked LD matrix")
      ("block-info", po::value<string>(&ldBlockInfoFile),"LD matrix file")
      ("ldm", po::value<string>(&ldmatrixFile),"input LD matrix floder")
      ("merge-block-ldm-info", "merge blocked LD matrix")
      ("block", po::value<unsigned>(&includeBlock), "block ")
      ("keep-block",po::value<string>(),"keep block string")
      ("make-ldm-eigen","a low-rank LD matrix")
      // way 2. generate LD block from genotype directly.
      ("make-eigen","Generate low-rank LD matrix for gene region")

      // Generate gene low-rank LD matrix
      ("make-eigen-gene","Generate low-rank LD matrix for gene region")
      ("gene-snp-pair",po::value<string>(),"Generate low-rank LD matrix for gene region using gene-snp-pairs.")
      ("merge-eigen-gene","merge multiple low-rank LD matrix.")
      ("eigen-list", po::value<string>(&geneListFile), "read file lists.")

      // read blockwise low-rank LD matrix via ldm-eigen and gene low-rank LD matrix via ldm-eigen-gene
      ("ldm-eigen", po::value<string>(&eigenMatrixFile),"low-rank LD matrix file")
      ("ldm-eigen-gene",po::value<string>(&geneEigenMatrixFile),"read a low-rank matrix for gene region")
      ("ldm-eigen-cutoff", po::value<string>(), "prefix of PLINK .fam, .bim, .bed files")
      ("ldm-eigen-gene-cutoff",po::value<string>() ,"read gene region low-rank matrix")

      // gwas summary
      ("gwas-summary",po::value<string>(&gwasSummaryFile),"read gwas summary statistics.")
      ("impute-summary", "read gwas summary dataset")
      ("merge-block-gwas-summary", "merge block gwas summary dataset")
      
      // eqtl summary
      ("beqtl-summary", po::value<string>(&eqtlSummaryFile), "read eqtl BESD format dataset")
      ("beqtl-summary-gz",po::value<string>(&eqtlSummaryQueryFile),"read query format")
      ("add-gene-n", po::value<string>(&geneSamSizeFile),"read gene sample size file.")
      ("keep-gene", po::value<string>(),"Specifiy a list of genes to be included in the analysis.")
      ("remove-gene", po::value<string>(),"Specifiy a list of genes to be excluded in the analysis.")
      ("keep-one-gene", po::value<string>(&includeSpecificGeneID),"Specifiy one gene to be included in the analysis.")
      ("eqtl-flist",po::value<string>(),"read flist file ")

      // data management for ukb-ppp-protein data
      ("rsid", po::value<string>(&rsidMapFiles),"protein rsid")
      ("promap", po::value<string>(&proteinMapFile),"protein map")
      ("proPath", po::value<string>(&proteinPath),"protein pqtl path")
      ("grch",po::value<string>(),"GrCh38 (hg38) or GrCh37 (hg19)")
      ("make-besd","generate besd format.")
      // ("make-besd-from-query","generate besd format based on query format")
      ("make-anno","read besd and generate annotation file")
      ("merge-besd","merge multiple besd formatfiles.")
      ("make-query","Query moQTL detail info based on various conditions.")
      ("make-query-ma","Query moQTL detail info ma format on various conditions.")
      ("make-query-sumstats","Query moQTL detail info ma format on various conditions.")
      ("resample-xqtl-N",po::value<int>(), "resample xQTL effect and se given new sample size")
      ("resample-method-type",po::value<int>(), "resample xQTL effect and se given new sample size")
      ("beqtl-list", po::value<string>(&eqtlSummaryFile), "read file lists.")
      ("anno-binary","annotation is continuous or binary")
      ("pvalue", po::value <double> (&pValueThreshold), "Specified QTL/moQTL p-value threshold. The default value is 1.")
      ("extract-snp-p", po::value <double> (), "Extract a subset of SNPs with p < a threshold.")
      ;

    po::options_description additional("Additional options:");
    additional.add_options()
      // error-checking
      ("noMapCheck", "disable automatic check of genetic map scale")
       ;

    // phenotype simulation parameters
    po::options_description hidden("Hidden options:");
    hidden.add_options()
      ("debug","debug mode")
      ("format-data","output formatted dataset.")
      ("allowX", "enable chrX analysis (now always enabled)");

    ///////////////////////////////////////////////////////////////////
    // Step2: combine mentioned options
    po::options_description visible("Options");
    visible.add(typical).add(indLevel).add(sumLevel).add(additional);

    po::options_description all("All options");
    all.add(typical).add(indLevel).add(sumLevel).add(additional).add(hidden);
    all.add_options()
      ("bad-args", po::value< vector <string> >(), "bad args")
      ;

    po::positional_options_description positional_desc;
    positional_desc.add("bad-args", -1); // for error-checking command line

    // Step3. boost cmd setting
    // po::variables_map vm;
    po::command_line_parser cmd_line(argc, argv);
    cmd_line.options(all);
    cmd_line.style(po::command_line_style::default_style ^ po::command_line_style::allow_guessing);
    cmd_line.positional(positional_desc);
    po::store(cmd_line.run(), vm);
    try {
      po::store(cmd_line.run(), vm);
      
      // if (vm.count("help") || vm.count("helpFull")) {
      //   if (vm.count("helpFull"))
      //     LOGGER << visible << endl;
	    //   else
	    //     LOGGER << typical << endl;
      //   exit(0);
      // }
      if (vm.count("help")) {
	        LOGGER << typical << endl;
          LOGGER << indLevel << endl;
          LOGGER << sumLevel << endl;
        exit(0);
      }
      po::notify(vm);  // throws an error if there are any problems
      //////////////////////////////////////////////////////
      //////////////////////////////////////////////////////
      //// The following is flags of main part
      //////////////////////////////////////////////////////
      //////////////////////////////////////////////////////
      /// Basic option
      if(vm.count("out")){
        title = vm["out"].as<string>();
      } 

      if(vm.count("debug")){
        diagnosticMode = true;
      }
      if(vm.count("format-data")){
        dataMode = true;
      }
      // Phenotype-related settings
      if(vm.count("pheno")) {
        phenotypeFile = vm["pheno"].as<string>();
      }
      // if(vm.count("keep")){
      //   excludeSnpFile = vm["keep"].as<string>();
      // }
      //// data management
      // individual --gene
      // summary --beqtl-summary --beql-summary-gz
      if( vm.count("gene") || vm.count("beqtl-summary") || vm.count("beqtl-summary-gz") ){
        haveXqtlDataBool = true;
      }


      if(vm.count("make-anno")){
        analysisType = "DataManagement";
        makeAnnoBool = true;
      }
      if(vm.count("make-besd")){
        analysisType = "DataManagement";
        makeBesdBool = true;
      }
      // if(vm.count("make-besd-from-query")){
      //   analysisType = "DataManagement";
      //   makeBesdBool = true;
      // }
      if(vm.count("merge-besd")){
        analysisType = "DataManagement";
        mergeBesdBool = true;
      }
      if(vm.count("merge-eigen-gene")){
        analysisType = "DataManagement";
        mergeEigenGeneBool = true;
      }
      if(vm.count("make-query")){
        analysisType = "DataManagement";
        makeQueryBool = true;
      }

      if(vm.count("make-query-ma")){
        analysisType = "DataManagement";
        makeQueryCojoMaBool = true;
      }
      if(vm.count("make-query-sumstats")){
        analysisType = "DataManagement";
        makeSumstatsFormatBool = true;
      }

      if(vm.count("resample-xqtl-N")){
        analysisType = "DataManagement";
        makeXqtlNBool = true;
        reEffectSampleSize = vm["resample-xqtl-N"].as<int>();
      }
      if(vm.count("resample-method-type")){
        analysisType = "DataManagement";
        reSamType = vm["resample-method-type"].as<int>();
      }

      if(vm.count("extract-snp-p")){
        analysisType = "DataManagement";
        pValueThreshold = vm["extract-snp-p"].as<double>();
      }
      if(vm.count("eqtl-flist")){
        analysisType = "DataManagement";
        geneInfoFile = vm["eqtl-flist"].as<string>();
      }

      if(vm.count("write-mcmc-txt")){
        writeTxtPosterior = true;
      }

      if(vm.count("anno-binary")){
        isAnnoBinary = true;
      }
      if(vm.count("grch")){
        grchType = vm["grch"].as<string>();
        if(grchType != "hg19" && grchType != "hg38"){
          LOGGER.e(0,"Error: grchType " + grchType + ", but grch type should be either [hg19] or [hg38].");
        }
      }
      /// Genotype-related settings
      if(vm.count("bfile")) {
        bedFile = vm["bfile"].as<string>();
      }
      if(vm.count("bfile-gene")){
        bedGeneFile = vm["bfile-gene"].as<string>();
      }
      
      if(vm.count("wind")) {
        windowWidth = unsigned((vm["wind"].as <double> ()) * Megabase);
      }
      if(vm.count("cis-wind")) {
        // For cis-analysis the distance between the probe midpoint and SNP genomic location 
        // was up to 1 Mb and for trans-analysis the distance was more than 5 Mb or 
        //  the probe and SNP were on different chromosomes.
        // Kasela S, Kisand K, Tserel L, Kaleviste E, Remm A, Fischer K, Esko T, Westra HJ, 
        // Fairfax BP, Makino S, Knight JC. Pathogenic implications for autoimmune mechanisms 
        // derived by comparative eQTL analysis of CD4+ versus CD8+ T cells. PLoS genetics. 2017 Mar 1;13(3):e1006643.
        cisRegionWind = ((vm["cis-wind"].as <double> ()) * Megabase);
        useCisWinBool = true;
      }
      /// LD matrix management
      if(vm.count("make-full-ldm")){
        analysisType = "LDmatrix";
        outLDmatType = "full";
      }
      if(vm.count("make-band-ldm")){
        analysisType = "LDmatrix";
        outLDmatType = "band";
      }
      if(vm.count("make-shrunk-ldm")){
        analysisType = "LDmatrix";
        outLDmatType = "shrunk";
      }
      if(vm.count("make-sparse-ldm")){
        analysisType = "LDmatrix";
        outLDmatType = "sparseshrunk";
      }
      if(vm.count("make-block-ldm")){
        analysisType = "LDmatrix";
        outLDmatType = "block";
      }

      //// ldm
      if(vm.count("ldm")){
        ldmatrixFile = vm["ldm"].as<string>();
      }
      if(vm.count("mldm")){
        ldmatrixFile = vm["mldm"].as<string>();
        multiLDmat = true;
      }
      if(vm.count("ldm-eigen")){
        eigenMatrixFile = vm["ldm-eigen"].as<string>();
      }
      if(vm.count("ldm-eigen-cutoff")){
        Gadget::Tokenizer strvec;
        strvec.getTokens(vm["ldm-eigen-cutoff"].as<string>(), " ,");
        eigenCutoff.resize(strvec.size());
        for (unsigned j=0; j<strvec.size(); ++j) eigenCutoff[j] = stod(strvec[j]);
        std::sort(eigenCutoff.data(), eigenCutoff.data()+eigenCutoff.size());
      }
      if(vm.count("block-info")){
        ldBlockInfoFile = vm["block-info"].as<string>();
      }
      if(vm.count("block")){
        includeBlock = vm["block"].as <unsigned> ();
      }
      if(vm.count("keep-block")){
        includeBlockID = vm["keep-block"].as<string> ();
      }
      if(vm.count("merge-block-ldm-info")){
        analysisType = "LDmatrix";
        mergeLdm = true;
        outLDmatType = "block";
      }

      // make low-rank LD matrix using genotype directly
      if(vm.count("make-ldm-eigen")){
        analysisType = "LDmatrixEigen";
      }
      if(vm.count("make-eigen")){
        analysisType = "LDmatrixEigen";
      }
      //// gene information
      if(vm.count("make-eigen-gene")){
        analysisType = "LDmatrixEigen";
      }

      if(vm.count("gene-snp-pair")){
        geneListFile = vm["gene-snp-pair"].as<string>();
      }

      if(vm.count("keep-gene")){
        includeGeneFile = vm["keep-gene"].as<string>();
      }
      if(vm.count("remove-gene")){
        excludeGeneFile = vm["remove-gene"].as<string>();
      }

      //// gwas summary
      // if(vm.count("gwas-summary")){
      //   gwasSummaryFile = vm["gwas-summary"].as<string>();
      // }

      if(vm.count("impute-summary")){
        analysisType = "ImputeSumStats";
        imputeSummary = true;
      }
      if(vm.count("merge-block-gwas-summary")){
        analysisType = "MergeGwasSummary";
        outLDmatType = "block";
      }

      /// eqtl summary
      if(vm.count("ldm-eigen-gene")){
        geneEigenMatrixFile = vm["ldm-eigen-gene"].as<string>();
      }

      if(vm.count("ldm-eigen-gene-cutoff")){
        Gadget::Tokenizer strvec;
        strvec.getTokens(vm["ldm-eigen-gene-cutoff"].as<string>(), " ,");
        geneEigenCutoff.resize(strvec.size());
        for (unsigned j=0; j<strvec.size(); ++j) geneEigenCutoff[j] = stod(strvec[j]);
        std::sort(geneEigenCutoff.data(), geneEigenCutoff.data() + geneEigenCutoff.size());
      }



      setThread();

      //////////////////////////////////////////////////////
      //////////////////////////////////////////////////////
    // }catch (po::error &e) {
    //   cerr << "ERROR: " << e.what() << endl << endl;
    //   cerr << visible << endl;
    //   return false;
    // }
    } catch (const string &err_msg) {
        LOGGER.e(0, err_msg);
        return false;
    } catch (const char *err_msg) {
        LOGGER.e(0, string(err_msg));
        return false;
    }
    return true;
  } // End of definition Option::boostProgramOptionsRoutine
