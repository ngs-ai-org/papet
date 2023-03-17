#include <applications/ApplicationKineticsWig.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <pbbam/CompositeBamReader.h>           // GenomicIntervalCompositeBamReader
#include <boost/program_options.hpp>            // variable_map, options_descriptions


#include <ngsaipp/utility/string_utility.hpp>           // ngsai::split()
#include <ngsaipp/io/bed_io.hpp>                        // ngsai::BedReader, ngsai::BedRecord
#include <ngsaipp/epigenetics/CcsKineticExtractor.hpp>  // ngsai::CcsKineticExtractor
#include <ngsaipp/genome/constants.hpp>                 // ngsai::genome::strand


namespace po = boost::program_options ;


ngsai::app::ApplicationKineticsWig::
                ApplicationKineticsWig(int argc,
                                       char** argv)
    : ApplicationInterface(argc, argv),
      m_path_bed(),
      m_paths_bam(),
      m_prefix_out(),
      m_win_size()
{   int parsing = this->parseOptions() ;
    if(parsing == this->getExitCodeSuccess())
    {   m_is_runnable = true ; }
    else
    {   m_is_runnable = false ; }
}


ngsai::app::ApplicationKineticsWig::
                ~ApplicationKineticsWig()
{ ; }


int
ngsai::app::ApplicationKineticsWig::run()
{   
    if(not this->isRunnable())
    {   return this->getExitCodeError() ; }

    try
    {   this->createWigTracks() ; }
    catch(const std::exception& e)
    {   std::cerr << "Error! something occured while " 
                     "creating the tracks:"
                  << std::endl
                  << e.what() 
                  << std::endl ;
    }

    return this->getExitCodeSuccess() ;
}


int
ngsai::app::ApplicationKineticsWig::parseOptions()
{
    // check arguments were given
    if(m_argc == 1)
    {   std::cerr << "Error! no options given"
                  << std::endl ; 
        return this->getExitCodeError() ;
    }

    // help messages
    std::string desc_msg =  "\n"
                            "Usage : kinetics-wig [options]"
                            "\n"
                            "\tkinetics-wig creates 4 WIG track displaying the\n"
                            "\tper position average IPD or PWD signal around\n"
                            "\tCpGs. For each CpG, the average signal is\n"
                            "\tcomputed from the pileup of reads mapping in\n"
                            "\tthe region. Both strands are treated\n"
                            "\tindependently.\n"
                            "\tThe track is printed on stdout.\n"
                            "\tWritten by Romain Groux, November 2022\n\n" ;
    std::string opt_help_msg = "Produces this help message." ;
    std::string opt_bam_msg  = "A coma separated list of paths to the bam "
                               "files containing the mapped PacBio CCS from "
                               "which the IPD/PWD track must be computed." ;
    std::string opt_bed_msg  = "The path to the bed file containing the " 
                               "genomic coordinates of the CpGs of interest.";
    std::string opt_out_msg  = "A path prefix to use to write the results "
                               "files. In total, 4 resulting files will be "
                               "created with this prefix : "
                               "<prefix>_IPDfw.wig, <prefix>_IPDrv.wig, "
                               "<prefix>_PWDfw.wig and <prefix>_PWDrv.wig "
                               "containing the IPD and PWD forward and "
                               "reverse track respectively." ;
    std::string opt_win_msg  = "The size of the window (in bp) around the "
                               "CpGs in which the average kinetic signal will "
                               "be computed." ;

    // option parser
    std::string path_bam("") ;
    std::string path_bed("") ;
    std::string path_out("") ;
    size_t win_size(0) ;
    po::variables_map vm ;
    po::options_description desc(desc_msg) ;
    desc.add_options()
        ("help,h",  opt_help_msg.c_str())
        ("bam",     po::value<std::string>(&(path_bam)), 
                    opt_bam_msg.c_str())
        ("bed",     po::value<std::string>(&(path_bed)), 
                    opt_bed_msg.c_str())
        ("out",     po::value<std::string>(&(path_out)),
                    opt_out_msg.c_str())
        ("winSize", po::value<size_t>(&(win_size)), 
                    opt_win_msg.c_str()) ;

    // parse
    try
    {   po::store(po::parse_command_line(m_argc,
                                         m_argv, 
                                         desc), vm) ;
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
    else if(path_bed == "")
    {   std::cerr <<"no bed file given (--bed)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(path_out == "")
    {   std::cerr <<"no prefix for results file given "
                    "(--out)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(win_size <= 0)
    {   std::cerr <<"window size must be > 0 (--winSize)" 
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    else if(win_size % 2 == 0)
    {   std::cerr <<"window size must be odd (--winSize)" 
                  << std::endl ;
        return this->getExitCodeError() ;
    }    

    // split bam paths and check the files
    std::vector<std::string> paths_bam = 
                        ngsai::split(path_bam, ',') ;
    for(const auto& path_bam : paths_bam)
    {
        if(this->checkBamFile(path_bam) != 
            this->getExitCodeSuccess())
        {   return this->getExitCodeError() ; }
    }

    // check BED file
    try
    {   ngsai::BedReader reader(path_bed) ; }
    catch(const std::exception& e)
    {
        std::cerr << "Error! could not open " 
                  << path_bed 
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    
    // set fields
    m_path_bed = path_bed ;
    m_paths_bam = paths_bam ;
    m_prefix_out = path_out ;
    m_win_size = win_size ;

    return this->getExitCodeSuccess() ;
}


void
ngsai::app::ApplicationKineticsWig::createWigTracks() const
{
    // file track paths
    std::string path_ipd_fw = 
        std::string(m_prefix_out).append("_IPDfw.wig") ;
    std::string path_ipd_rv = 
        std::string(m_prefix_out).append("_IPDrv.wig") ;
    std::string path_pwd_fw = 
        std::string(m_prefix_out).append("_PWDfw.wig") ;
    std::string path_pwd_rv = 
        std::string(m_prefix_out).append("_PWDrv.wig") ;

    // open streams
    std::ofstream f_ipd_fw(path_ipd_fw) ;
    std::ofstream f_ipd_rv(path_ipd_rv) ;
    std::ofstream f_pwd_fw(path_pwd_fw) ;
    std::ofstream f_pwd_rv(path_pwd_rv) ;

    // track definition lines
    f_ipd_fw << "track type=wiggle_0 "
                "name=\"IPD fw\" "
                "description=\"fw IPD averages\" "
                "visibility=full "
                "autoScale=off "
                "color=50,150,255 "
                "priority=1" 
             << std::endl ;
    f_ipd_rv << "track type=wiggle_0 "
                "name=\"IPD rv\" "
                "description=\"rv IPD averages\" "
                "visibility=full "
                "autoScale=off "
                "color=50,150,255 "
                "priority=3" 
             << std::endl ;
    f_pwd_fw << "track type=wiggle_0 "
                "name=\"PWD fw\" "
                "description=\"fw PWD averages\" "
                "visibility=full "
                "autoScale=off "
                "color=0,200,100 "
                "priority=2" 
              << std::endl ;
    f_pwd_rv << "track type=wiggle_0 "
                "name=\"PWD rv\" "
                "description=\"rv PWD averages\" "
                "visibility=full "
                "autoScale=off "
                "color=0,200,100 "
                "priority=4" 
             << std::endl ;

    // half the window size
    size_t win_size_half = m_win_size / 2 ;
    
    // bam readers
    PacBio::BAM::GenomicIntervalCompositeBamReader 
                                reader_bam(m_paths_bam) ;

    // extract kinetics at all regions and fill histograms
    ngsai::BedRecord cpg ;
    ngsai::BedReader bed_reader(m_path_bed) ;
    ngsai::CcsKineticExtractor extractor ;
    PacBio::BAM::BamRecord ccs ;
    std::string chrom_current("") ;
    while(bed_reader.getNext(cpg))
    {  
        // only keep CpG on + and compute coordinates of 
        // + and - strand CpGs
        if(cpg.strand == ngsai::genome::REVERSE)
        {   continue ; }

        // include one declaration line per chromosome 
        if(cpg.chrom != chrom_current)
        {   chrom_current = cpg.chrom ;
            f_ipd_fw << "variableStep chrom=" 
                     << cpg.chrom << " "
                     << "span=1" 
                     << std::endl ;
            f_ipd_rv << "variableStep chrom=" 
                     << cpg.chrom << " "
                     << "span=1" 
                     << std::endl ;
            f_pwd_fw << "variableStep chrom=" 
                     << cpg.chrom << " "
                     << "span=1" 
                     << std::endl ;
            f_pwd_rv << "variableStep chrom=" 
                     << cpg.chrom << " "
                     << "span=1" << std::endl ;
        }

        // CpG window on + strand
        ngsai::BedRecord window_p(cpg) ;
        window_p.start -= win_size_half ;
        window_p.end   += win_size_half - 1 ;
        // CpG window on - strand
        ngsai::BedRecord window_m(cpg) ;
        window_m.strand = ngsai::genome::REVERSE ;
        window_m.start -= win_size_half - 1 ;
        window_m.end   += win_size_half ;

        PacBio::BAM::GenomicInterval interval(cpg.chrom, 
                                              cpg.start,
                                              cpg.end) ;

        // compute mean IPD and PWD from CCS for this region
        double n_p = 0. ;
        double n_m = 0. ;
        std::vector<double> ipds_m_p(m_win_size, 0.) ;
        std::vector<double> ipds_m_m(m_win_size, 0.) ;
        std::vector<double> pwds_m_p(m_win_size, 0.) ;
        std::vector<double> pwds_m_m(m_win_size, 0.) ;
        reader_bam.Interval(interval) ;
        while(reader_bam.GetNext(ccs))
        {   
            // extract kinetics on + strand
            if(extractor.extract(ccs, window_p))
            {   
                std::vector<uint16_t> ipds = 
                                    extractor.getIPD() ;
                std::vector<uint16_t> pwds = 
                                    extractor.getPWD() ;
                for(size_t i=0; i<ipds_m_p.size(); i++)
                {   ipds_m_p[i] += 
                            static_cast<double>(ipds[i]) ;
                    pwds_m_p[i] += 
                            static_cast<double>(pwds[i]) ;
                }
                n_p += 1. ;
            }

            // extract kinetics on - strand
            if(extractor.extract(ccs, window_m))
            {   
                std::vector<uint16_t> ipds = 
                            extractor.getIPD() ;
                std::vector<uint16_t> pwds = 
                            extractor.getPWD() ;
                for(size_t i=0; i<ipds_m_m.size(); i++)
                {   ipds_m_m[i] += 
                            static_cast<double>(ipds[i]) ;
                    pwds_m_m[i] += 
                            static_cast<double>(pwds[i]) ;
                }
                n_m += 1. ;
            }
        }
        for(size_t i=0; i<ipds_m_p.size(); i++)
        {   if(n_p > 0.)
            {   ipds_m_p[i] /= n_p ;
                pwds_m_p[i] /= n_p ;
            }
            if(n_m > 0.)
            {   ipds_m_m[i] /= n_m ;
                pwds_m_m[i] /= n_m ;
            }
        }

        // print IPD + strand track

        if(n_p > 0.)
        {   for(size_t i=0, pos=window_p.start; 
                i<ipds_m_p.size(); 
                i++, pos++)
            {   f_ipd_fw << pos + 1 << ' ' 
                         << ipds_m_p[i] 
                         << std::endl ;
            }

            // print PWD + strand track
            for(size_t i=0, pos=window_p.start; 
                i<pwds_m_p.size(); 
                i++, pos++)
            {   f_pwd_fw << pos + 1 << ' ' 
                         << pwds_m_p[i] 
                         << std::endl ;
            }
        }
        if(n_m > 0.)
        {   // print IPD - strand track
            for(int i=ipds_m_m.size()-1, pos=window_p.start; 
                i>=0; 
                i--, pos++)
            {   f_ipd_rv << pos + 2 << ' ' 
                         << ipds_m_m[i] 
                         << std::endl ;
            }

            // print PWD - strand track
            for(int i=pwds_m_m.size()-1, pos=window_p.start; 
                i>=0; 
                i--, pos++)
            {   f_pwd_rv << pos + 2 << ' ' 
                         << pwds_m_m[i] 
                         << std::endl ; }
        }
    }
    f_ipd_fw.close() ;
    f_ipd_rv.close() ;
    f_pwd_fw.close() ;
    f_pwd_rv.close() ;
}
