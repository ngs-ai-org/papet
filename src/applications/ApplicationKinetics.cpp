#include <applications/ApplicationKinetics.hpp>

#include <iostream>
#include <iomanip>
#include <string>
#include <boost/program_options.hpp>       // variable_map, options_descriptions
#include <boost/archive/text_iarchive.hpp> // boost::archive::text_iarchive
#include <boost/serialization/utility.hpp> // std::pair serialization
#include <pbbam/BamFile.h>                 // BamFile
#include <pbbam/CompositeBamReader.h>      // CompositeBamReader
#include <pbbam/BamRecord.h>               // BamRecord

#include <ngsaipp/io/bed_io.hpp>                        // ngsai::BedReader, ngsai::BedRecord 
#include <ngsaipp/io/utility.hpp>
#include <ngsaipp/utility/string_utility.hpp>           // ngsai::split()
#include <ngsaipp/epigenetics/CcsKineticExtractor.hpp>  // ngsai::CcsKineticExtractor
#include <ngsaipp/epigenetics/KmerMap.hpp>              // ngsai::KmerMap
#include <ngsaipp/epigenetics/model_utility.hpp>        // ngsai::normalize_kinetics
#include <ngsaipp/genome/constants.hpp>                 // ngsai::genome::strand
#include <applications/utilities.hpp>                   // ngsai::app::print_vector


namespace po = boost::program_options ;


ngsai::app::ApplicationKinetics::ApplicationKinetics(
                        int argc,
                        char** argv)
    : ApplicationInterface(argc, argv),
      m_bam_files(),
      m_path_bed(),
      m_win_size(0),
      m_kmermap(nullptr)
{   int parsing = this->parseOptions() ;
    if(parsing == this->getExitCodeSuccess())
    {   m_is_runnable = true ; }
    else
    {   m_is_runnable = false ; }
}


ngsai::app::ApplicationKinetics::~ApplicationKinetics()
{
    if(m_kmermap != nullptr)
    {   delete m_kmermap ;
        m_kmermap = nullptr ;
    }
}


int
ngsai::app::ApplicationKinetics::run()
{   if(not this->isRunnable())
    {   return this->getExitCodeError() ; }

    // number of bases to consider around the C of a CpG to get window
    size_t win_size_half = m_win_size / 2 ;

    // set precision for floating values printing
    std::cout << std::setprecision(5) ;

    // headers
    this->printHeader(std::cout, '\t') ;

    // read features overlapping CpG
    ngsai::BedReader bed_reader(m_path_bed) ;
    ngsai::BedRecord cpg ;
    ngsai::CcsKineticExtractor extractor ;
    PacBio::BAM::BamRecord ccs ;
    while(bed_reader.getNext(cpg))
    {   // consider only regular chromosomes
        if(cpg.chrom.find("chr") != 0)
        {  continue ; } 
        
        PacBio::BAM::GenomicInterval interval(cpg.chrom, 
                                              cpg.start,
                                              cpg.end) ;
        PacBio::BAM::GenomicIntervalCompositeBamReader reader(interval,
                                                              m_bam_files);

        // coordinates of the window, on the reference
        // region centered on the C of the CpG
        // CpG are encoded as [start,end) with respect to FORWARD strand
        //      start  end
        //        |     |
        //  ... N C p G N ... forward strand
        //  ... N G p C N ... reverse strand
        //        |     |
        //      start  end
        ngsai::BedRecord window(cpg) ;
        // on + strand C corresponds to start pos of bed entry
        if(cpg.strand == ngsai::genome::FORWARD)
        {   window.start -= win_size_half ;
            window.end   += win_size_half - 1 ;
        }
        // on - strand C corresponds to end-1 pos of bed entry
        else if(cpg.strand == ngsai::genome::REVERSE)
        {   window.start -= win_size_half - 1 ;
            window.end   += win_size_half ;
        }
        // no orientation -> cannot extract feature
        else
        {   continue ; }

        while(reader.GetNext(ccs))
        {    try
            {   if(extractor.extract(ccs, window))
                {   std::vector<uint16_t> ipd = extractor.getIPD() ;
                    std::vector<uint16_t> pwd = extractor.getPWD() ;
                    std::string           seq = extractor.getSequence() ;
                    char strand = ngsai::genome::strand_to_char(window.strand) ;
                    // print read features
                    std::cout << window.chrom  << '\t'
                              << window.start  << '\t'
                              << window.end    << '\t'
                              << strand        << '\t'
                              << seq           << '\t';
                    if(m_kmermap != nullptr)
                    {   auto ratios = 
                            ngsai::normalize_kinetics(seq, 
                                                      ipd, 
                                                      pwd, 
                                                      *m_kmermap) ;
                        print_vector(std::cout, ratios.first,'\t')  << '\t' ;
                        print_vector(std::cout, ratios.second, '\t') << std::endl ;
                    }
                    else
                    {   print_vector(std::cout, ipd, '\t') << '\t' ;
                        print_vector(std::cout, pwd, '\t') << std::endl ;
                    }
                }
            }
            catch(const std::exception& e)
            {
                std::cerr << "Error! something occured "
                             "while treating "
                          << std::endl
                          << "bed region : " 
                          << cpg 
                          << std::endl
                          << "read name : " 
                          << ccs.FullName()     
                          << std::endl
                          << "read mapping start  : " 
                          << ccs.ReferenceStart() 
                          << std::endl
                          << "read mapping end    : "
                          << ccs.ReferenceEnd() 
                          << std::endl 
                          << "read mapping strand : "
                          << ccs.AlignedStrand()
                          << std::endl 
                          << "error message : "
                          << e.what()
                          << std::endl ;
                return this->getExitCodeError() ;
            }
        }
    }
    return this->getExitCodeSuccess() ;
}


int
ngsai::app::ApplicationKinetics::parseOptions()
{   
    // check arguments were given
    if(m_argc == 1)
    {   std::cerr << "Error ! no options given"
                  << std::endl ; 
        return this->getExitCodeError() ;
    }

    // help messages
    std::string desc_msg =  "\n"
                            "Usage : kinetics [options] > [FILE]"
                            "\n"
                            "\tkinetics is an application to extract\n"
                            "\tinterpulse duration (IPDs) and pulse widths\n"
                            "\t(PWDs) from mapped PacBio CCS reads that\n"
                            "\toverlap a given set of genomic regions\n"
                            "\tspecified in a BED file. Only reads that align\n"
                            "\tperfectly over the regions are reported.\n"
                            "\tThe results are printed on stdout in tsv \n"
                            "\tformat. The first row is a header. Then, each\n"
                            "\tline contains per read sequence, IPDs and PWDs.\n"
                            "\tWritten by Romain Groux, November 2022\n\n" ;
    std::string opt_help_msg  = "Produces this help message." ;
    std::string opt_bam_msg = "A coma separated list of paths to the bam "
                              "files containing the mapped PacBio CCS of "
                              "interest." ;
    std::string opt_bed_msg = "The path to the bed file containing the " 
                              "genomic regions of interest.";
    std::string opt_map_msg = "The path to a file containing a kmerMap to " 
                              "normalize the IPD and PWD kinetic signal. If "
                              "none is given, the raw kinetic values will be "
                              "returned." ;
    std::string opt_win_msg = "The size in bp of the windows from which the "
                              "CCS features will be extracted. It must be "
                              "odd. The window will be centered on the center "
                              "of the genomic regions of interst." ;

    // option parser
    std::string path_bam("") ;
    std::string path_bed("") ;
    std::string path_map("") ;
    bool normalization = false ;  // will become true if a map is given
    int win_size = -1 ;
    po::variables_map vm ;
    po::options_description desc(desc_msg) ;
    desc.add_options()
        ("help,h",  opt_help_msg.c_str())
        ("bam",     po::value<std::string>(&(path_bam)), 
                    opt_bam_msg.c_str())
        ("bed",     po::value<std::string>(&(path_bed)), 
                    opt_bed_msg.c_str())
        ("model",   po::value<std::string>(&(path_map)), 
                    opt_map_msg.c_str())
        ("winSize", po::value<int>(&(win_size)), 
                    opt_win_msg.c_str()) ;

    // parse
    try
    {   po::store(po::parse_command_line(m_argc, 
                                         m_argv, 
                                         desc), 
                                         vm) ;
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
    {   std::cerr <<"Error! no bam files given (--bam)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(path_bed == "")
    {   std::cerr <<"Error! no bed file given (--bed)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(path_map != "")
    {   normalization = true ; }
    else if(win_size <= 0)
    {   std::cerr <<"Error! window size must be > 0 "
                    "(--winSize)" 
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(win_size % 2 == 0)
    {   std::cerr << "Error! window size must be odd "
                     "(--winSize)" 
                  << std::endl ;
        return this->getExitCodeError() ;
    }

    // check the bam files
    std::vector<std::string> paths_bam = 
                        ngsai::split(path_bam, ',') ;
    // check bam file indexes
    for(const auto& path_bam : paths_bam)
    {
        if(this->checkBamFile(path_bam) != 
            this->getExitCodeSuccess())
        {   return this->getExitCodeError() ; }
        m_bam_files.push_back(
            PacBio::BAM::BamFile(path_bam)) ;
    }

    // load KmerMap
    if(normalization)
    {   if(this->loadKmerMap(path_map) != 
           this->getExitCodeSuccess())
        {   return this->getExitCodeError() ; } 
    }

    // sets fields
    m_path_bed = path_bed ;
    m_win_size = win_size ;

    return this->getExitCodeSuccess() ;
}


int
ngsai::app::ApplicationKinetics::loadKmerMap(
                            const std::string& path)
{   
    if(m_kmermap != nullptr)
    {   delete m_kmermap ; 
        m_kmermap = nullptr ;
    }

    m_kmermap = new ngsai::KmerMap(1) ;
    try
    {   // read dump
        std::ifstream f_in(path) ;
        boost::archive::text_iarchive arch_in(f_in) ;
        arch_in >> (*m_kmermap) ;
        f_in.close() ; 
    }
    catch(std::exception& e)
    {   std::cerr << "Error! could not load the KmerMap "
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


void
ngsai::app::ApplicationKinetics::printHeader(
                                std::ostream& stream,
                                char separator) const
{   std::vector<std::string> headers = {"chrom", 
                                        "start",
                                        "end",
                                        "strand",
                                        "seq"} ;
    // IPD columns
    for(size_t i=0; i<m_win_size; i++)
    {   char ipd[64] ;
        sprintf(ipd, "ipd_%zu", i+1) ;
        headers.push_back(ipd) ;
    }

    // PWD columns
    for(size_t i=0; i<m_win_size; i++)
    {   char ipd[64] ;
        sprintf(ipd, "pwd_%zu", i+1) ;
        headers.push_back(ipd) ;
    }

    print_vector(stream, headers, separator) << std::endl ;
}
