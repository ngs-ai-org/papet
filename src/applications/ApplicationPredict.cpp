#include <applications/ApplicationPredict.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <iomanip>
#include <future>                          // std::promise, std::future
#include <boost/program_options.hpp>       // variable_map, options_descriptions
#include <boost/archive/text_iarchive.hpp> // boost::archive::text_iarchive
#include <boost/serialization/utility.hpp> // std::pair serialization
#include <pbbam/BamFile.h>                  // PacBio::BAM::BamFile
#include <pbbam/CompositeBamReader.h>       // PacBio::BAM::GenomicIntervalCompositeBamReader

#include <ngsaipp/utility/string_utility.hpp>       // ngsai::split()
#include <ngsaipp/epigenetics/KineticModel.hpp>
#include <ngsaipp/epigenetics/RawKineticModel.hpp>
#include <ngsaipp/epigenetics/NormalizedKineticModel.hpp>
#include <ngsaipp/epigenetics/PairWiseKineticModel.hpp>
#include <ngsaipp/epigenetics/PairWiseNormalizedKineticModel.hpp>
#include <ngsaipp/epigenetics/DiPositionKineticModel.hpp>
#include <ngsaipp/epigenetics/DiPositionNormalizedKineticModel.hpp>
#include <ngsaipp/epigenetics/KineticClassifier.hpp>
#include <ngsaipp/io/bed_io.hpp>
#include <ngsaipp/genome/CpGRegion.hpp>
#include <ngsaipp/genome/constants.hpp>
#include <ngsaipp/parallel/ThreadPool.hpp>          // ngsai::ThreadPool


namespace po = boost::program_options ;


ngsai::app::ApplicationPredict::ApplicationPredict(
                        int argc,
                        char** argv)
    : ApplicationInterface(argc, argv),
      m_paths_bam(),
      m_classifier(),
      m_cpgs(),
      m_prob_meth(0.),
      m_threads_n(0)
{   int parsing = this->parseOptions() ;
    if(parsing == this->getExitCodeSuccess())
    {   m_is_runnable = true ; }
    else
    {   m_is_runnable = false ; }
}


ngsai::app::ApplicationPredict::~ApplicationPredict()
{ ; }


int
ngsai::app::ApplicationPredict::run()
{   
    if(not this->isRunnable())
    {   return this->getExitCodeError() ; }

    // sub-sets of CpGs to treat by each thread
    std::vector<std::pair<size_t,size_t>> slices = 
                    ngsai::ThreadPool::split_range(0,
                                            m_cpgs.size(),
                                            m_threads_n) ;

    // get promises and futures
    std::vector<std::promise<std::list<double>>> 
                                    promises(m_threads_n) ;
    std::vector<std::future<std::list<double>>> 
                                    futures(m_threads_n) ;
    for(size_t i=0; i<m_threads_n; i++)
    {   futures[i] = promises[i].get_future() ; }

    // thread pool
    ngsai::ThreadPool threads(m_threads_n) ;

    // distribute to threads
    for(size_t i=0; i<m_threads_n; i++)
    {   auto slice = slices[i] ;
        threads.addJob(
            std::move(
                std::bind(
                    &ApplicationPredict::predictRoutine,
                    this,
                    slice.first,
                    slice.second,
                    std::ref(promises[i])))) ;
    }

    // wait until all thread is done
    threads.join() ;
    
    // print results
    for(size_t n=0; n<m_threads_n; n++)
    {   
        // CpG indices treated by the thread
        const auto& slice = slices[n] ;
        size_t from = slice.first ;
        size_t to   = slice.second ;
        
        // prob of meth. computed by the thread
        auto& future = futures[n] ;
        std::list<double> probs = future.get() ;
        
        // point to 1st CpG treated by thread and 
        // corresponding prob
        std::list<double>::iterator prob = probs.begin() ;
        std::list<ngsai::genome::CpGRegion>::iterator 
            cpg = m_cpgs.begin() ;
        std::advance(cpg, from) ;
        
        for(size_t i=from; i<to; i++)
        {   std::cout 
                << cpg->chrom << '\t'
                << cpg->start << '\t'
                << cpg->end   << '\t'
                << ""         << '\t'
                << std::setprecision(4) << *prob  << '\t'
                << ngsai::genome::strand_to_char(cpg->strand)
                << std::endl ;
            cpg++ ;
            prob++ ;
        }
    }

    return this->getExitCodeSuccess() ;
}


int
ngsai::app::ApplicationPredict::parseOptions()
{
    // check arguments were given
    if(m_argc == 1)
    {   std::cerr << "Error ! no options given"
                  << std::endl ; 
        return this->getExitCodeError() ;
    }

    // help messages
    std::string desc_msg =  "\n"
                            "Usage : predict [options] > [FILE]"
                            "\n"
                            "\tPredicts the CpG methylation status and \n"
                            "\treturns the results on stdout in BED 6 format.\n"
                            "\tThe score field contains the methylation\n"
                            "\tprobability.\n"
                            "\tWritten by Romain Groux, October 2022\n\n" ;
    std::string opt_help_msg  = "Produces this help message." ;
    std::string opt_bam_msg = "A coma separated list of paths to the bam\n"
                              "files containing the mapped PacBio CCS of\n"
                              "interest." ;
    std::string opt_bed_msg = "The path to a bed file containing the\n" 
                              "coordinates of the CpGs interest.";
    std::string opt_model_meth_msg     = "The path to the file containing the\n"
                                        "methylated kinetic model to use." ;
    std::string opt_model_unmeth_msg  = "The path to the file containing the\n"
                                        "unmethylated kinetic model to use." ;
    std::string opt_prob_msg  = "The prior probability of methylation\n"
                                 "for any CpG. It must belong to [0,1].\n" 
                                 "0.5 by default.";
    std::string opt_thread_msg = "The number of threads, by default 1." ;


    // option parser
    std::string path_bam("") ;
    std::string path_bed("") ;
    std::string path_mod_m("") ;
    std::string path_mod_u("") ;
    double prob_meth = 0.5 ;
    size_t n_threads(1) ;

    po::variables_map vm ;
    po::options_description desc(desc_msg) ;
    desc.add_options()
        ("help,h",      opt_help_msg.c_str())
        ("bam",         po::value<std::string>(&(path_bam)), 
                        opt_bam_msg.c_str())
        ("bed",         po::value<std::string>(&(path_bed)), 
                        opt_bed_msg.c_str())
        ("modelMeth",   po::value<std::string>(&(path_mod_m)), 
                        opt_model_meth_msg.c_str())
        ("modelUnmeth", po::value<std::string>(&(path_mod_u)), 
                        opt_model_unmeth_msg.c_str())
        ("prob",        po::value<double>(&(prob_meth)), 
                        opt_prob_msg.c_str())
        ("thread",      po::value<size_t>(&(n_threads)), 
                        opt_thread_msg.c_str()) ;
    
    // parse
    try
    {   po::store(po::parse_command_line(m_argc, 
                                         m_argv, 
                                         desc), vm) ;
        po::notify(vm) ;
    }
    catch(std::invalid_argument& e)
    {   std::string msg = std::string("Error! Invalid "
                                      "option given\n") + 
                          std::string(e.what()) ;
        return this->getExitCodeError() ;
    }
    catch(...)
    {   std::cerr << "Error! an unknown error occured "
                     "while parsing the options" 
                  << std::endl ; 
        return this->getExitCodeError() ;
    }

    // display help if needed
    bool help = vm.count("help") ;
    if(help)
    {   std::cout << desc << std::endl ;
        return this->getExitCodeError() ;
    }

    // check options
    if(path_bam == "")
    {   std::cerr <<"Error! no bam file given (--bam)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(path_bed == "")
    {   std::cerr <<"Error! no bed file given (--bed)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(path_mod_m == "")
    {   std::cerr <<"Error! no methylated kinetic model "
                    "file given (--modelMeth)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(path_mod_u == "")
    {   std::cerr <<"Error! no unmethylated kinetic model "
                    "file given (--modelUnmeth)"
                  << std::endl ;
        return getExitCodeError() ;
    }
    else if(prob_meth < 0. or 
            prob_meth > 1.)
    {   std::cerr << "Error! prior methylation probability "
                     "must belong to [0,1] (--prob)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(n_threads == 0)
    {   std::cerr << "Error! number of threads must by > 0 "
                     "(--thread)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }

    // load models and transform them into log densities
    if(this->loadModels(path_mod_m, path_mod_u) !=
       this->getExitCodeSuccess())
    {   return this->getExitCodeError() ; }
   
    // load the CpG BED regions
    if(this->loadBed(path_bed) !=
       this->getExitCodeSuccess())
    {   return this->getExitCodeError() ; }

    // check bam files
    std::vector<std::string> paths_bam = 
                            ngsai::split(path_bam, ',') ;
    for(const auto& path_bam : paths_bam)
    {   if(this->checkBamFile(path_bam) !=
            this->getExitCodeSuccess()) 
        {   return this->getExitCodeError() ; }
    }

    // set remaining fields
    m_paths_bam = paths_bam ;
    m_prob_meth = prob_meth ;
    m_threads_n = n_threads ;

    return this->getExitCodeSuccess() ;
}


int
ngsai::app::ApplicationPredict::loadModels(
                const std::string& path_model_meth,
                const std::string& path_model_unmeth)
{  
    ngsai::KineticModel* model_meth(nullptr) ;
    ngsai::KineticModel* model_unmeth(nullptr) ; 

    // load methylated kinetic model
    if(ngsai::endswith(path_model_meth,
                    ".rawkineticmodel"))
    {   model_meth = new ngsai::RawKineticModel() ;
        model_meth->load(path_model_meth) ;
    }
    else if(ngsai::endswith(path_model_meth, 
                    ".normalizedkineticmodel"))
    {   model_meth = new ngsai::NormalizedKineticModel() ;
        model_meth->load(path_model_meth) ;
    }
    else if(ngsai::endswith(path_model_meth,
                    ".pairwisekineticmodel"))
    {   model_meth = new ngsai::PairWiseKineticModel() ;
        model_meth->load(path_model_meth) ;
    }
    else if(ngsai::endswith(path_model_meth, 
                    ".pairwisenormalizedkineticmodel"))
    {   model_meth = 
            new ngsai::PairWiseNormalizedKineticModel() ;
        model_meth->load(path_model_meth) ;
    }
    else if(ngsai::endswith(path_model_meth,
                    ".dipositionkineticmodel"))
    {   model_meth = new ngsai::DiPositionKineticModel() ;
        model_meth->load(path_model_meth) ;
    }
    else if(ngsai::endswith(path_model_meth, 
                    ".dipositionnormalizedkineticmodel"))
    {   model_meth = 
            new ngsai::DiPositionNormalizedKineticModel() ;
        model_meth->load(path_model_meth) ;
    }
    else
    {   std::cerr 
            << "Error! could not determine the type of "
                "the methylated kinetic model stored in "
            << path_model_meth << std::endl ;
        return this->getExitCodeError() ; 
    }

    // transform methylated kinetic model to log density
    if(not model_meth->isInit())
    {   std::cerr 
            << "Error! methylated kinetic signal model is " 
               " not initialised"
            << std::endl ;
        return this->getExitCodeError() ;
    }
    if(not model_meth->isDensity())
    {   model_meth->density() ; }
    if(not model_meth->isLog())
    {   model_meth->log() ; }

    // load unmethylated kinetic model
    if(ngsai::endswith(path_model_unmeth,
                    ".rawkineticmodel"))
    {   model_unmeth = new ngsai::RawKineticModel() ;
        model_unmeth->load(path_model_unmeth) ;
    }
    else if(ngsai::endswith(path_model_unmeth, 
                    ".normalizedkineticmodel"))
    {   model_unmeth = new ngsai::NormalizedKineticModel() ;
        model_unmeth->load(path_model_unmeth) ;
    }
    else if(ngsai::endswith(path_model_unmeth,
                    ".pairwisekineticmodel"))
    {   model_unmeth = new ngsai::PairWiseKineticModel() ;
        model_unmeth->load(path_model_unmeth) ;
    }
    else if(ngsai::endswith(path_model_unmeth, 
                    ".pairwisenormalizedkineticmodel"))
    {   model_unmeth = 
            new ngsai::PairWiseNormalizedKineticModel() ;
        model_unmeth->load(path_model_unmeth) ;
    }
    else if(ngsai::endswith(path_model_unmeth,
                    ".dipositionkineticmodel"))
    {   model_unmeth = new ngsai::DiPositionKineticModel() ;
        model_unmeth->load(path_model_unmeth) ;
    }
    else if(ngsai::endswith(path_model_unmeth, 
                    ".dipositionnormalizedkineticmodel"))
    {   model_unmeth = 
            new ngsai::DiPositionNormalizedKineticModel() ;
        model_unmeth->load(path_model_unmeth) ;
    }
    else
    {   std::cerr 
            << "Error! could not determine the type of "
                "the methylated kinetic model stored in "
            << path_model_unmeth << std::endl ;
        return this->getExitCodeError() ; 
    }

    // transform unmethylated kinetic model to log density
    if(not model_unmeth->isInit())
    {   std::cerr 
            << "Error! unmethylated kinetic signal model " 
               "is not initialised"
            << std::endl ;
        return this->getExitCodeError() ;
    }
    if(not model_unmeth->isDensity())
    {   model_unmeth->density() ; }
    if(not model_unmeth->isLog())
    {   model_unmeth->log() ; }

    // build classifier
    try
    {    m_classifier.setModels(model_meth,
                                model_unmeth) ;
    }
    catch(const std::exception& e)
    {   std::cerr << "Error! could not construct a "
                     "kinetic signal classifier:" 
                  << std::endl 
                  << e.what() << std::endl ;
        return this->getExitCodeError() ;
    }

    return this->getExitCodeSuccess() ;
}


int 
ngsai::app::ApplicationPredict::loadBed(
    const std::string& path_bed)
{   
   try
    {   ngsai::BedRecord bed_record ;
        ngsai::BedReader bed_reader(path_bed) ;
        while(bed_reader.getNext(bed_record))
        {   // only keep fw strand since rv have same coords 
            if(bed_record.strand != 
                    ngsai::genome::strand::FORWARD)
            {   continue ; }

            bed_record.strand = 
                ngsai::genome::strand::UNORIENTED ;
            m_cpgs.push_back(
                ngsai::genome::CpGRegion(bed_record)) ;
        }
    }
    catch(const std::exception& e)
    {   std::cerr << "Error! could not load BED regions:"
                  << std::endl 
                  << e.what() << std::endl ;
        return this->getExitCodeError() ;
    }
    return this->getExitCodeSuccess() ;
}


void 
ngsai::app::ApplicationPredict::predictRoutine(
            size_t from,
            size_t to,
            std::promise<std::list<double>>& promise) const
{   
    double prob_unmeth = 1. - m_prob_meth ;

    // will contain the prob of methylation of each CpG treated
    std::list<double> probs ;

    PacBio::BAM::BamRecord record_bam ;
    PacBio::BAM::GenomicIntervalCompositeBamReader reader_bam(m_paths_bam) ;

    std::list<ngsai::genome::CpGRegion>::const_iterator cpg = m_cpgs.begin();
    std::advance(cpg, from) ;
    for(size_t i=from; i<to; i++)
    {   // extract overlapping CCSs
        std::list<PacBio::BAM::BamRecord> ccss ;
        PacBio::BAM::GenomicInterval interval(cpg->chrom, 
                                              cpg->start,
                                              cpg->end) ;
        reader_bam.Interval(interval) ;
        while(reader_bam.GetNext(record_bam))
        {   ccss.push_back(record_bam) ; }
        std::pair<double,double> prob = 
                    m_classifier.classify(*cpg, 
                                          ccss,
                                          m_prob_meth,
                                          prob_unmeth) ;
        probs.push_back(prob.first) ;
        cpg++ ;
    }

    promise.set_value(std::move(probs)) ;
}
