#include <applications/ApplicationModelKinetic.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <random>
#include <functional>
#include <typeinfo>
#include <fstream>
#include <limits>
#include <thread>                               // std::thread
#include <boost/program_options.hpp>            // variable_map, options_descriptions
#include <boost/archive/text_oarchive.hpp>                   // boost::archive::text_oarchive
#include <boost/serialization/utility.hpp>                   // std::pair serialization
#include <boost/archive/text_iarchive.hpp>                   // boost::archive::text_oarchive

#include <ngsaipp/io/bed_io.hpp>                                     // ngsai::BedReader, ngsai::BedRecord
#include <ngsaipp/epigenetics/RawKineticModel.hpp>                   // ngsai::RawKineticModel
#include <ngsaipp/epigenetics/NormalizedKineticModel.hpp>            // ngsai::NormalizedKineticModel
#include <ngsaipp/epigenetics/DiPositionKineticModel.hpp>            // ngsai::DiPositionKineticModel
#include <ngsaipp/epigenetics/DiPositionNormalizedKineticModel.hpp>  // ngsai::DiPositionNormalizedKineticModel
#include <ngsaipp/epigenetics/PairWiseKineticModel.hpp>              // ngsai::PairWiseKineticModel
#include <ngsaipp/epigenetics/PairWiseNormalizedKineticModel.hpp>    // ngsai::PairWiseNormalizedKineticModel
#include <ngsaipp/epigenetics/model_utility.hpp>                     // train_KineticModdel()
#include <ngsaipp/utility/string_utility.hpp>                        // ngsai::split()
#include <ngsaipp/parallel/ThreadPool.hpp>                           // ngsai::ThreadPool::split_range()

namespace po = boost::program_options ;


ngsai::app::ApplicationModelKinetic::
                ApplicationModelKinetic(
                                    int argc,
                                    char** argv)
    : ApplicationInterface(argc, argv),
      m_mode(modes::undefined),
      m_paths_bam(),
      m_path_out(),
      m_path_kmermap(),
      m_size(0),
      m_nb_bins(0),
      m_xmin(std::numeric_limits<double>::min()),
      m_xmax(std::numeric_limits<double>::max()),
      m_pseudo_counts(0.),
      m_nb_threads(1),
      m_cpgs(),
      m_kmermap(nullptr),
      m_models()
{   int parsing = this->parseOptions() ;
    if(parsing == this->getExitCodeSuccess())
    {   m_is_runnable = true ; }
    else
    {   m_is_runnable = false ; }
}


ngsai::app::ApplicationModelKinetic::
                ~ApplicationModelKinetic()
{
    if(m_kmermap != nullptr)
    {   delete m_kmermap ;
        m_kmermap = nullptr ;
    }
    this->freeKineticModels() ;
}


int
ngsai::app::ApplicationModelKinetic::run()
{   
    if(not this->isRunnable())
    {   return this->getExitCodeError() ; }

    // threads
    std::vector<std::thread> threads;

    // sub-sets of regions to train models on
    std::vector<std::pair<size_t,size_t>> slices = 
                ngsai::ThreadPool::split_range(0, 
                                            m_cpgs.size(),
                                            m_nb_threads) ;
    // start all threads
    // -------------- threads start --------------
    for(size_t i=0; i<m_nb_threads; i++)
    {   
        // models have been allocated and parameters set
        // already
        void (*pf)(ngsai::KineticModel*,
                   const std::vector<ngsai::BedRecord>&,
                   size_t,
                   size_t,
                   const std::vector<std::string>&) = 
                                ngsai::train_KineticModel ;
        threads.push_back(
                std::thread(pf,
                            m_models[i],
                            std::ref(m_cpgs),
                            size_t(slices[i].first),
                            size_t(slices[i].second),
                            std::ref(m_paths_bam))) ;
    }
    for(auto& thread : threads)
    {   if(thread.joinable())
        {   thread.join() ; }
    }
    // -------------- threads end --------------

    // aggregate models
    for(size_t i=1; i<m_models.size(); i++)
    {   m_models[0]->add(*(m_models[i])) ; }

    // serialize model
    m_models[0]->save(m_path_out) ;

    // free memory
    this->freeKineticModels() ;
    
    return this->getExitCodeSuccess() ;
}


int
ngsai::app::ApplicationModelKinetic::parseOptions()
{
    // check arguments were given
    if(m_argc == 1)
    {   std::cerr << "Error! no options given"
                  << std::endl ; 
        return this->getExitCodeError() ;
    }

    // help messages
    std::string desc_msg =  "\n"
                            "Usage : model-kinetic [type] [options]"
                            "\n"
                            "\tTrains a kinetic signal model. [type] defines\n"
                            "\tthe type of kinetic signal model. The possible\n"
                            "\tvalues of [type] are:\n"
                            "\t-'raw': simply computes the distribution of IPD\n"
                            "\t and PWD at each position in the window.\n"
                            "\t-'raw-norm': computes the normalized IPD and PWD'\n"
                            "\t distributions at each position in the window.\n"
                            "\t-'diposition': computes the distributions of IPD\n"
                            "\t and PWD signal at each position in the window\n"
                            "\t as a function of its direct neighbour.\n"
                            "\t-'diposition-norm': computes the distributions of\n"
                            "\t normalized IPD and PWD signal at each position\n"
                            "\t in the window as a function of its direct neighbor\n"
                            "\t-'pairwise': computes all distributions of IPD and\n"
                            "\t PWD as a function of the signal at another\n"
                            "\tposition in the window.\n"
                            "\t-'pairwise-norm': computes all distributions of\n"
                            "\t normalized IPD and PWD as a function of the\n"
                            "\tsignal at another position in the window.\n"
                            "\tThe trained model is serialized in the given\n"
                            "\tfile.\n"
                            "\tWritten by Romain Groux, November 2022\n\n" ;
    std::string opt_help_msg = "Produces this help message." ;
    std::string opt_bam_msg  = "A coma separated list of paths to the bam "
                               "files containing the mapped PacBio CCS of "
                               "interest." ;
    std::string opt_bckg_msg = "Only for normalized models: the path to a "
                               "file containing the background model to use. "
                               "It must contain a serialized KmerMap." ;
    std::string opt_bed_msg  = "The path to the bed file containing the " 
                               "genomic coordinates of the CpGs of interest.";
    std::string opt_out_msg  = "The path to file in which the kinetic model "
                               "will be saved." ;
    std::string opt_size_msg = "The size of the model, the length of the "
                               "signal window to model, in bp." ;
    std::string opt_nbin_msg = "The number of bins in each histogram." ;
    std::string opt_xmin_msg = "The lower limit of the lower bin in each "
                               "histogram." ;
    std::string opt_xmax_msg = "The upper limit of the upper bin in each "
                               "histogram." ;
    std::string opt_pcnt_msg = "A number of counts that will be added to "
                               "each bin in each histogram, by default 0." ;
    std::string opt_thread_msg = "The number of threads, by default 1." ;

    // option parser
    std::string path_bam("") ;
    std::string path_bed("") ;
    std::string path_out("") ;
    std::string path_kmermap("") ;
    size_t size(0) ;
    size_t nb_bins(0) ;
    double xmin(std::numeric_limits<double>::min()) ;
    double xmax(std::numeric_limits<double>::max()) ;
    double pseudo_counts(0.) ;
    size_t n_threads(1) ;

    po::variables_map vm ;
    po::options_description desc(desc_msg) ;
    desc.add_options()
        ("help,h",     opt_help_msg.c_str())
        ("bam",        po::value<std::string>(&(path_bam)), 
                       opt_bam_msg.c_str())
        ("bed",        po::value<std::string>(&(path_bed)), 
                       opt_bed_msg.c_str())
        ("out",        po::value<std::string>(&(path_out)), 
                       opt_out_msg.c_str())
        ("background", po::value<std::string>(&(path_kmermap)), 
                       opt_bckg_msg.c_str())
        ("size",       po::value<size_t>(&(size)), 
                       opt_size_msg.c_str())
        ("nbin",       po::value<size_t>(&(nb_bins)), 
                       opt_size_msg.c_str())
        ("xmin",       po::value<double>(&(xmin)), 
                       opt_xmin_msg.c_str())
        ("xmax",       po::value<double>(&(xmax)), 
                       opt_xmax_msg.c_str())
        ("pseudocount",  
                    po::value<double>(&(pseudo_counts)), 
                    opt_pcnt_msg.c_str())
        ("thread",  
                    po::value<size_t>(&(n_threads)), 
                    opt_thread_msg.c_str());

    // parse
    try
    {   po::store(po::parse_command_line(m_argc, 
                                         m_argv,
                                         desc), vm) ;
        po::notify(vm) ;
    }
    catch(const std::invalid_argument& e)
    {   std::cerr << "Error! Invalid option given:"
                  << std::endl  
                  << e.what()
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    catch(const std::exception& e)
    {   std::cerr << "Error! an unknown error occured "
                     "while parsing the options:" 
                  << std::endl
                  << e.what()
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
    {   std::cerr <<"no bam file given (--bam)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(path_bed == "")
    {   std::cerr <<"no bed file given (--bed)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(path_out == "")
    {   std::cerr <<"no output file given (--out)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(size == 0)
    {   std::cerr <<"invalid model size given (--size)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(nb_bins == 0)
    {   std::cerr <<"invalid number of bins given (--nbin)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(xmin == std::numeric_limits<double>::min())
    {   std::cerr <<"invalid x-axis minimum given (--xmin)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(xmax == std::numeric_limits<double>::max())
    {   std::cerr <<"invalid x-axis maximum given (--xmax)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(xmin >= xmax)
    {   std::cerr <<"xmin must be smaller than xmax "
                    "(--xmin --xmax)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(pseudo_counts < 0.)
    {   std::cerr <<"pseudo counts must be >= 0 "
                    "(--pseudocount)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(n_threads == 0)
    {   std::cerr <<"number of threads must by > 0 "
                    "(--thread)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }

    // type of model to train
    std::string opt_mode(m_argv[1]) ;
    if(opt_mode == "raw")
    {   m_mode = modes::raw ; }
    else if(opt_mode == "raw-norm")
    {   m_mode = modes::raw_norm ; }
    else if(opt_mode == "diposition")
    {   m_mode = modes::diposition ; }
     else if(opt_mode == "diposition-norm")
    {   m_mode = modes::diposition_norm ; }
    else if(opt_mode == "pairwise")
    {   m_mode = modes::pairwise ; }
    else if(opt_mode == "pairwise-norm")
    {   m_mode = modes::pairwise_norm ; }
    else
    {   std::cerr << "Error! " << m_argv[1] << " "
                  << "does not indicate a model type. The "
                  << "accepted values are: raw, raw-norm "
                  << "diposition, diposition-norma, "
                  << "pairwise, pairwise-norm"
                  << std::endl ;
        return this->getExitCodeError() ;
    }

    if(((m_mode == modes::raw_norm) or
             (m_mode == modes::diposition_norm) or
             (m_mode == modes::pairwise_norm)) and 
            (path_kmermap == ""))
    {   std::cerr <<"no background model file given "
                    "(--background)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }

    // check bam files
    std::vector<std::string> paths_bam = 
                            ngsai::split(path_bam, ',') ;
    for(const auto& path: paths_bam)
    {   if(this->checkBamFile(path) != 
                    this->getExitCodeSuccess())
        {   return this->getExitCodeError() ; }
    }

    // set fields
    m_paths_bam = paths_bam ;
    m_path_out = path_out ; 
    m_path_kmermap = path_kmermap ;
    m_size = size ; 
    m_nb_bins = nb_bins ; 
    m_xmin = xmin ;
    m_xmax = xmax ;
    m_pseudo_counts = pseudo_counts ;
    m_nb_threads = n_threads ;

    // load BED file
    if(this->loadBed(path_bed))
    {   return this->getExitCodeError() ; }

    // allocate model memory
    if(this->allocateKineticModels())
    {   return this->getExitCodeError() ; }

    return this->getExitCodeSuccess() ;
}


int
ngsai::app::ApplicationModelKinetic::
    loadKmerMap(const std::string& path)
{   
    if(m_kmermap != nullptr)
    {   delete m_kmermap ;
        m_kmermap = nullptr ;
    }

    try
    {
        std::ifstream f_in(path) ;
        boost::archive::text_iarchive arch_in(f_in) ;
        m_kmermap = new ngsai::KmerMap(1) ;
        arch_in >> *(m_kmermap) ;
        f_in.close() ;
    }
    catch(const std::exception& e)
    {
        std::cerr << "Error! could not load the KmerMap "
                     "from "
                  << path
                  << " :"
                  << std::endl 
                  << e.what() 
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    return this->getExitCodeSuccess() ;
}

int
ngsai::app::ApplicationModelKinetic::allocateKineticModels()
{   // ensured no kinetic model allocated
    if(this->freeKineticModels() != 
            this->getExitCodeSuccess())
    {   return this->getExitCodeError() ; }

    try
    {
        m_models = std::vector<ngsai::KineticModel*>
                    (m_nb_threads, nullptr) ;

        if(m_mode == modes::raw)
        {   
            for(auto& ptr : m_models)
            {   ptr = new ngsai::RawKineticModel() ; }
        }
        else if(m_mode == modes::raw_norm)
        {   if(this->loadKmerMap(m_path_kmermap) != 
               this->getExitCodeSuccess())
            {   return this->getExitCodeError() ; }
            
            for(auto& ptr : m_models)
            {   ptr = new ngsai::NormalizedKineticModel(
                                            *m_kmermap) ;
            }
            delete m_kmermap ;
            m_kmermap = nullptr ;
        }
        else if(m_mode == modes::diposition)
        {   for(auto& ptr: m_models)
            {   ptr = 
                    new ngsai::DiPositionKineticModel() ;
            }
        }
        else if(m_mode == modes::diposition_norm)
        {   if(this->loadKmerMap(m_path_kmermap) != 
               this->getExitCodeSuccess())
            {   return this->getExitCodeError() ; }
            
            for(auto& ptr: m_models)
            {   ptr = new 
                    ngsai::DiPositionNormalizedKineticModel(
                                            *m_kmermap) ;
            }
            delete m_kmermap ;
            m_kmermap = nullptr ;
        }
        else if(m_mode == modes::pairwise)
        {   for(auto& ptr: m_models)
            {   ptr = new ngsai::PairWiseKineticModel() ; }
        }
        else if(m_mode == modes::pairwise_norm)
        {   if(this->loadKmerMap(m_path_kmermap) != 
               this->getExitCodeSuccess())
            {   return this->getExitCodeError() ; }
            
            for(auto& ptr : m_models)
            {   ptr = new 
                    ngsai::PairWiseNormalizedKineticModel(
                                            *m_kmermap) ;
            }
            delete m_kmermap ;
            m_kmermap = nullptr ;
        }
        else
        {   std::cerr << "Error! could not determine the "
                         "type of kinetic signal model "
                         "to train"
                      << std::endl ;
            return this->getExitCodeError() ;
        }

        // set model parameters
        for(size_t i=0; i<m_models.size(); i++)
        {   
            // only add pseudo counts to 1st partial model 
            // because after when summing them we want to  
            // have pseudo counts added only once, not 
            // once per thread
            double pc = (i==0) ? m_pseudo_counts : 0 ;
            
            m_models[i]->setParameters(m_size,
                                       m_xmin,
                                       m_xmax,
                                       m_nb_bins,
                                       pc) ;
        }
    }
    catch(const std::exception& e)
    {
        std::cerr << "Error! Something occured while "
                     "allocating the KineticModel: "
                  << std::endl
                  << e.what()
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    return this->getExitCodeSuccess() ;
}


int
ngsai::app::ApplicationModelKinetic::freeKineticModels()
{   try
    {
        for(auto& ptr : m_models)
        {   if(ptr != nullptr)
            {   delete ptr ;
                ptr = nullptr ;
            }
        }
        m_models.clear() ;
    }
    catch(const std::exception& e)
    {   std::cerr << "Error! something occured while "
                     "free the KineticModel: "
                  << std::endl 
                  << e.what() 
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    return this->getExitCodeSuccess() ;
}


int 
ngsai::app::ApplicationModelKinetic::loadBed(
    const std::string& path_bed)
{   
   try
    {   ngsai::BedRecord bed_record ;
        ngsai::BedReader bed_reader(path_bed) ;
        while(bed_reader.getNext(bed_record))
        {   // only keep fw strand since rv have same coords 
            // if(bed_record.strand != 
            //         ngsai::genome::strand::FORWARD)
            // {   continue ; }

            // bed_record.strand = 
            //     ngsai::genome::strand::UNORIENTED ;
            m_cpgs.push_back(bed_record) ;
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
