#include <applications/ApplicationModelSequence.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <utility>       // std::make_pair
#include <stdexcept>

#include <pbbam/BamFile.h>
#include <pbbam/BamReader.h>
#include <pbbam/CompositeBamReader.h>
#include <pbbam/BamRecord.h>

#include <boost/program_options.hpp>       //variable_map, options_descriptions
#include <boost/archive/text_oarchive.hpp>  // boost::archive::text_oarchive

#include <ngsaipp/utility/string_utility.hpp>       // ngsai::split()
#include <ngsaipp/dna/dna_utility.hpp>              // ngsai::get_reverse_complement()
#include <ngsaipp/epigenetics/KmerMap.hpp>


namespace po = boost::program_options ;


ngsai::app::ApplicationModelSequence::
    ApplicationModelSequence(int argc,char** argv)
    : ApplicationInterface(argc, argv),
      m_paths_bam(),
      m_path_out(),
      m_kmermap(nullptr)
{   int parsing = this->parseOptions() ;
    if(parsing == this->getExitCodeSuccess())
    {   m_is_runnable = true ; }
    else
    {   m_is_runnable = false ; }
}


ngsai::app::ApplicationModelSequence::
    ~ApplicationModelSequence()
{
    if(m_kmermap != nullptr)
    {   delete m_kmermap ;
        m_kmermap = nullptr ;
    }    
}


int
ngsai::app::ApplicationModelSequence::
    run()
{   if(not this->isRunnable())
    {   return this->getExitCodeError() ; }

    // construct model
    for(auto& path_bam : m_paths_bam)
    {   int code = this->updateKmerMap(path_bam) ;
        if(code != this->getExitCodeSuccess())
        {   return this->getExitCodeError() ; }
    }
    // average per nb of occurences
    for(auto iter=m_kmermap->begin(); 
             iter!=m_kmermap->end(); 
             iter++)
    {   uint32_t n = iter->first ;
        // otherwise div by 0
        if(n != 0)
        {   for(size_t j=0; j<iter->second.ipd.size(); j++)
            {   iter->second.ipd[j] /= n ;
                iter->second.pwd[j] /= n ;
            }
        }
    }
    // dump
    std::ofstream f_out(m_path_out);
    boost::archive::text_oarchive arch_out(f_out);
    arch_out << *m_kmermap ;
    f_out.close() ;

    return this->getExitCodeSuccess() ;
}


int
ngsai::app::ApplicationModelSequence::
    parseOptions()
{   
    if(m_kmermap != nullptr)
    {   delete m_kmermap ; 
        m_kmermap = nullptr ;
    }

    // check arguments were given
    if(m_argc == 1)
    {   std::cerr << "Error ! no options given"
                  << std::endl ; 
        return this->getExitCodeError() ;
    }

    // help messages
    std::string desc_msg =  "\n"
                            "Usage : model-sequence [options]"
                            "\n"
                            "\tmodel-sequence computes the expected average\n"
                            "\tkinetic signal (IPD and PWD) for each kmer \n"
                            "\tfrom the given CCS data. The results are\n" 
                            "\twritten to a file as a serialized KmerMap.\n"
                            "\tWritten by Romain Groux, November 2022\n\n" ;
    std::string opt_help_msg  = "Produces this help message." ;
    std::string opt_bam_msg = "A coma separated list of paths to the bam "
                              "files containing the PacBio CCS of "
                              "interest." ;
    std::string opt_out_msg = "The path to the file in which the KmerMap  "
                              "will be dumped." ; 
    std::string opt_win_msg = "The size of the kmers." ;

    // option parser
    std::string path_bam("") ;
    std::string path_out("") ;
    int kmer_size = -1 ;

    po::variables_map vm ;
    po::options_description desc(desc_msg) ;
    desc.add_options()
        ("help,h",  opt_help_msg.c_str())
        ("bam",     po::value<std::string>(&(path_bam)), 
                    opt_bam_msg.c_str())
        ("out",     po::value<std::string>(&(path_out)), 
                    opt_out_msg.c_str())
        ("kmer",    po::value<int>(&(kmer_size)), 
                    opt_win_msg.c_str()) ;

    // parse
    try
    {   po::store(po::parse_command_line(m_argc,
                                         m_argv, 
                                         desc), vm) ;
        po::notify(vm) ;
    }
    catch(std::invalid_argument& e)
    {   std::cerr << "Error! Invalid option given: " 
                  << std::endl 
                  << e.what() 
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    catch(const std::exception& e)
    {   std::cerr << "Error! an unknown error occured "
                     "while parsing the options: "
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
    {   std::cerr <<"Error! no bam file given (--bam)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(path_out == "")
    {   std::cerr <<"Error! no file to dump background "
                    "model given (--out)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(kmer_size <= 0)
    {   std::cerr <<"Error! kmer size must be > 0 (--kmer)" 
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(kmer_size % 2 == 0)
    {   std::cerr << "Error! kmer size must be odd (--kmer)" 
                  << std::endl ;
        return this->getExitCodeError() ;
    }

    // check bam files
    std::vector<std::string> paths_bam = 
                            ngsai::split(path_bam, ',') ;
    for(const auto& path_bam : paths_bam)
    {   if(this->checkBamFile(path_bam) !=
            this->getExitCodeSuccess()) 
        {   return this->getExitCodeError() ; }
    }

    // set fields
    m_paths_bam = paths_bam ;
    m_path_out = path_out ;
    m_kmermap = new ngsai::KmerMap(kmer_size) ;

    return this->getExitCodeSuccess() ;
}


int
ngsai::app::ApplicationModelSequence::
    updateKmerMap(const std::string& path_bam)
{
    try
    {
        // size_t n_seq = 0 ;
        PacBio::BAM::BamRecord record ;
        PacBio::BAM::BamReader reader(path_bam) ;
        while(reader.GetNext(record))
        {    
            // forward strand of the CCS
            // sequence
            std::string seq = 
                record.Sequence(
                    PacBio::BAM::Orientation::NATIVE) ;
            // IPDs
            auto ipds_16 = 
                record.ForwardIPD(
                    PacBio::BAM::Orientation::NATIVE).
                        Data() ;
            if(ipds_16.size() == 0)
            {   continue ; }
            std::vector<uint32_t> ipds(ipds_16.begin(), 
                                       ipds_16.end()) ;
            // PWDs
            auto pwds_16 = 
                record.ForwardPulseWidth(
                    PacBio::BAM::Orientation::NATIVE).
                        Data() ;
            std::vector<uint32_t> pwds(pwds_16.begin(), 
                                       pwds_16.end()) ;
            // insert in kmermap
            m_kmermap->insert(
                seq,
                ipds,
                pwds, 
                ngsai::KmerMap::insert_mode::SUM) ;

            // reverse strand of the CCS
            // sequence
            seq = ngsai::dna::get_reverse_complement(seq) ;
            // IPDs
            ipds_16 = 
                record.ReverseIPD(
                    PacBio::BAM::Orientation::NATIVE).
                    Data() ;
            if(ipds_16.size() == 0)
            {   continue ; }
            ipds = std::vector<uint32_t>(ipds_16.begin(), 
                                         ipds_16.end()) ;
            // PWDs
            pwds_16 = 
                record.ReversePulseWidth(
                    PacBio::BAM::Orientation::NATIVE).
                        Data() ;
            pwds = std::vector<uint32_t>(pwds_16.begin(), 
                                         pwds_16.end()) ;
            // insert in map
            m_kmermap->insert(
                seq, 
                ipds, 
                pwds, 
                ngsai::KmerMap::insert_mode::SUM) ;
        }
    }
    catch(const std::exception& e)
    {   std::cerr << "Error! something occured while "
                     "parsing "
                  << path_bam << " :" 
                  << std::endl 
                  << e.what()
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    return this->getExitCodeSuccess() ;
}
