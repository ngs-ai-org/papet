#include <applications/ApplicationKineticsKmer.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>                    // std::reverse()
#include <pbbam/BamReader.h>            // PacBio::BAM::BamReader
#include <pbbam/BamRecord.h>            // PacBio::BAM::BamRecord
#include <boost/program_options.hpp>    // PacBio::BAM::variable_map, options_descriptions
#include <unordered_map>                // std::unordered_map

#include <ngsaipp/dna/dna_utility.hpp>         // ngsai::get_reverse_complement()
#include <ngsaipp/utility/string_utility.hpp>  // ngsai::split()
#include <applications/utilities.hpp>          // ngsai::app::print_vector()

namespace po = boost::program_options ;


ngsai::app::ApplicationKineticsKmer::
                ApplicationKineticsKmer(int argc,
                                       char** argv)
    : ApplicationInterface(argc, argv),
      m_paths_bam(),
      m_kmer_size(0)
{   int parsing = this->parseOptions() ;
    if(parsing == this->getExitCodeSuccess())
    {   m_is_runnable = true ; }
    else
    {   m_is_runnable = false ; }
}


ngsai::app::ApplicationKineticsKmer::
                ~ApplicationKineticsKmer()
{ ; }


int
ngsai::app::ApplicationKineticsKmer::run()
{   
    if(not this->isRunnable())
    {   return this->getExitCodeError() ; }


    // half the kmer size
    size_t kmer_size_half = m_kmer_size / 2 ;
    
    // highest possible value for decoded IPD/PWD
    size_t frame_max_score = 952 ;

    // headers
    ngsai::app::print_vector(std::cout, 
                             this->generateHeaders(
                                     frame_max_score,
                                     frame_max_score),
                             '\t') ;
    std::cout << std::endl ;

    // the map to store the kmer score distributions
    // for each kmer, there is a vector counting the 
    // number of occurences of each possible scores.
    // There are frame_max_score + 1 possible values :
    // 0, 1, 2, ..., 951, 952 
    std::unordered_map<std::string,
                       std::vector<uint32_t>> map_ipd ;
    std::unordered_map<std::string,
                       std::vector<uint32_t>> map_pwd ;

    // extract kmer and kinetics
    PacBio::BAM::BamRecord ccs ;
    for(const auto& path_bam : m_paths_bam)
    {   
        // to hold CCS kinetics 
        std::vector<uint16_t> ipd ;
        std::vector<uint16_t> pwd ;
        std::string seq ;
        
        // to hold kmer kinetic value
        uint16_t ipd_kmer ;
        uint16_t pwd_kmer ;
        std::string kmer ;
        PacBio::BAM::BamReader bam_reader(path_bam) ;

        while(bam_reader.GetNext(ccs))
        {   
            try
            { 
                // forward strand
                // decoded IPDs
                ipd = ccs.ForwardIPD(
                    PacBio::BAM::Orientation::NATIVE).Data() ;
                // CCS with too few passes have no IPD
                if(ipd.size() == 0)
                {   continue ; }
                // PWDs
                pwd = ccs.ForwardPulseWidth(
                    PacBio::BAM::Orientation::NATIVE).Data() ;            
                // sequence
                seq = ccs.Sequence(
                    PacBio::BAM::Orientation::NATIVE) ;
                // [from,to) interval
                for(size_t from=0; 
                    from<seq.size() - m_kmer_size + 1;
                    from++)
                {   size_t to = from + m_kmer_size ;
                    size_t center = from + kmer_size_half ;
                    kmer = std::string(seq.begin() + from, 
                                       seq.begin() + to) ;
                    ipd_kmer = ipd[center] ;
                    pwd_kmer = pwd[center] ;

                    // store in map
                    if(map_ipd.find(kmer) == map_ipd.end())
                    {   // there are frame_max_score + 1 
                        // possible values:
                        // 0, 1, 2, ..., 951, 952
                        map_ipd.emplace(
                            kmer, 
                            std::vector<uint32_t>(
                                frame_max_score + 1, 0)) ;
                        map_pwd.emplace(
                            kmer,
                            std::vector<uint32_t>(
                                frame_max_score + 1 , 0)) ;
                    }
                    map_ipd.at(kmer)[ipd_kmer] += 1 ;
                    map_pwd.at(kmer)[pwd_kmer] += 1 ;
                }

                // reverse strand
                // decoded IPDs
                ipd = ccs.ReverseIPD(
                    PacBio::BAM::Orientation::NATIVE).Data() ;
                // PWDs
                pwd = ccs.ReversePulseWidth(
                    PacBio::BAM::Orientation::NATIVE).Data() ;
                // sequence
                seq = ngsai::dna::get_reverse_complement(seq) ;
                // [from,to) interval
                for(size_t from=0; 
                    from<seq.size() - m_kmer_size + 1;
                    from++)
                {   size_t to = from + m_kmer_size ;
                    size_t center = from + kmer_size_half ;
                    kmer = std::string(seq.begin() + from, 
                                       seq.begin() + to) ;
                    ipd_kmer = ipd[center] ;
                    pwd_kmer = pwd[center] ;

                    // store in map
                    if(map_ipd.find(kmer) == map_ipd.end())
                    {   // there are frame_max_score + 1 
                        // possible values:
                        // 0, 1, 2, ..., 951, 952
                        map_ipd.emplace(
                            kmer,
                            std::vector<uint32_t>(
                                frame_max_score + 1, 0)) ;
                        map_pwd.emplace(
                            kmer,
                            std::vector<uint32_t>(
                                frame_max_score + 1, 0)) ;
                    }
                    map_ipd.at(kmer)[ipd_kmer] += 1 ;
                    map_pwd.at(kmer)[pwd_kmer] += 1 ;
                }
            }
            catch(const std::exception& e)
            {
                std::cerr << "Error! something occured "
                             "while processing CCS "
                          << ccs.FullName() 
                          << " in "
                          << path_bam
                          << ":"
                          << std::endl
                          << e.what() 
                          << std::endl ;
                return this->getExitCodeError() ;
            }
        }
    }

    // print
    for(const auto pair: map_ipd)
    {   // kmer
        std::cout << pair.first << '\t' ;
        // IPD count distribution
        ngsai::app::print_vector(std::cout,
                                 pair.second,
                                 '\t') ;
        std::cout << '\t' ;
        // PWD count distribution
        ngsai::app::print_vector(std::cout, 
                                 map_pwd.at(pair.first),
                                 '\t') ;
        std::cout << std::endl ;
    }

    return this->getExitCodeSuccess() ; 
}


int
ngsai::app::ApplicationKineticsKmer::parseOptions()
{
    // check arguments were given
    if(m_argc == 1)
    {   std::cerr << "Error! no options given"
                  << std::endl ; 
        return this->getExitCodeError() ;
    }

    // help messages
    std::string desc_msg =  "\n"
                            "Usage : kinetics-kmer [options] > [FILE]\n"
                            "\n"
                            "\tkinetics-kmer computes the per kmer\n"
                            "\tdistribution of IPD and PWD values. The\n"
                            "\tresults are printed on stdout in tsv format.\n"
                            "\tThe first row contains the headers and the\n"
                            "\tfollowing rows, the data. Each row contains\n"
                            "\tthe kmer sequence followed by 953 IPD counts -\n"
                            "\tcorresponding to the counts of all possible\n"
                            "\tIPD values from 0 to 952 - and by 953 PWD\n" 
                            "\tcounts - corresponding to the counts of all\n"
                            "\tpossible PWD values from 0 to 952.\n"
                            "\tWritten by Romain Groux, November 2022\n\n" ;
    std::string opt_help_msg  = "Produces this help message." ;
    std::string opt_bam_msg  = "A coma separated list of paths to the bam "
                               "files containing the mapped PacBio CCS of "
                               "interest." ;
    std::string opt_kmer_msg = "The kmer length in base pairs." ;

    // option parser
    std::string path_bam("") ;
    size_t kmer_size(0) ;
    po::variables_map vm ;
    po::options_description desc(desc_msg) ;
    desc.add_options()
        ("help,h",  opt_help_msg.c_str())
        ("bam",     po::value<std::string>(&(path_bam)), 
                    opt_bam_msg.c_str())
        ("kmer,k",  po::value<size_t>(&(kmer_size)), 
                    opt_kmer_msg.c_str()) ;

    // parse
    try
    {   po::store(po::parse_command_line(m_argc, 
                                         m_argv, desc), vm) ;
        po::notify(vm) ;
    }
    catch(std::invalid_argument& e)
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
    else if(kmer_size <= 0)
    {   std::cerr <<"kmer size must be > 0 (--kmer)" 
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(kmer_size % 2 == 0)
    {   std::cerr <<"kmer size must be odd (--kmer)" 
                  << std::endl ;
        return this->getExitCodeError() ;
    }

    // split the bam paths
    std::vector<std::string> paths_bam = 
                        ngsai::split(path_bam, ',') ;
    std::vector<PacBio::BAM::BamFile> bam_files ;
    for(const auto& path_bam : paths_bam)
    {
        if(this->checkBamFile(path_bam) != 
            this->getExitCodeSuccess())
        {   return this->getExitCodeError() ; }
    }

    // set fields
    m_paths_bam = paths_bam ;
    m_kmer_size = kmer_size ;

    return this->getExitCodeSuccess() ;
}


std::vector<std::string> 
ngsai::app::ApplicationKineticsKmer::generateHeaders(
                                        size_t max_ipd,
                                        size_t max_pwd)
{
    std::vector<std::string> headers ;
    
    // kmer sequence column
    headers.push_back("seq") ;

    // IPD columns
    for(size_t i=0; i<=max_ipd; i++)
    {   char header[64] ;
        sprintf(header, "ipd_%zu", i) ;
        headers.push_back(header) ;
    }

    // PWD columns
    for(size_t i=0; i<=max_pwd; i++)
    {   char header[64] ;
        sprintf(header, "pwd_%zu", i) ;
        headers.push_back(header) ;
    }

    return headers ;
}
