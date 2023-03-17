#include <applications/ApplicationModelSequenceTxt.hpp>

#include <iostream>
#include <string>
#include <fstream>
#include <vector>

#include <boost/program_options.hpp>        //variable_map, options_descriptions
#include <boost/archive/text_iarchive.hpp>  // boost::archive::text_oarchive

#include <ngsaipp/epigenetics/KmerMap.hpp>
#include <applications/utilities.hpp>       // ngsai::app::print_vector


namespace po = boost::program_options ;


ngsai::app::ApplicationModelSequenceTxt::
    ApplicationModelSequenceTxt(int argc,char** argv)
    : ApplicationInterface(argc, argv),
      m_kmermap(nullptr)
{   int parsing = this->parseOptions() ;
    if(parsing == this->getExitCodeSuccess())
    {   m_is_runnable = true ; }
    else
    {   m_is_runnable = false ; }
}


ngsai::app::ApplicationModelSequenceTxt::
    ~ApplicationModelSequenceTxt()
{
    if(m_kmermap != nullptr)
    {   delete m_kmermap ;
        m_kmermap = nullptr ;
    }    
}


int
ngsai::app::ApplicationModelSequenceTxt::run()
{
    if(not this->isRunnable())
    {   return this->getExitCodeError() ; }

    // pretty print
    this->printKmerMap(std::cout) ;

    return this->getExitCodeSuccess() ;
}


int
ngsai::app::ApplicationModelSequenceTxt::
    parseOptions()
{
    // check arguments were given
    if(m_argc == 1)
    {   std::cerr << "Error ! no options given"
                  << std::endl ; 
        return this->getExitCodeError() ;
    }

    // help messages
    std::string desc_msg =  "\n"
                            "Usage : model-sequence-txt [options] > [FILE]"
                            "\n"
                            "\tConverts a sequence model (KmerMap) file in \n"
                            "\ttxt format."
                            "\tWritten by Romain Groux, November 2022\n\n" ;
    std::string opt_help_msg = "Produces this help message." ;
    std::string opt_mod_msg  = "The path to the file containing the KmerMap "
                               "to convert in tsv format." ;

    // option parser
    std::string path_mod("") ;

    po::variables_map vm ;
    po::options_description desc(desc_msg) ;
    desc.add_options()
        ("help,h",  opt_help_msg.c_str())
        ("model",   po::value<std::string>(&(path_mod)), 
                    opt_mod_msg.c_str()) ;
    
    // parse
    try
    {   po::store(po::parse_command_line(m_argc, 
                                         m_argv, desc), vm) ;
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
    if(path_mod == "")
    {   std::cerr <<"Error! no model file given (--model)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }

    // load KmerMap
    if(this->loadKmerMap(path_mod) != 
        this->getExitCodeSuccess())
    {   return this->getExitCodeError() ; }

    return this->getExitCodeSuccess() ;
}


int
ngsai::app::ApplicationModelSequenceTxt::
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


void
ngsai::app::ApplicationModelSequenceTxt::
    printKmerMap(std::ostream& stream) const
{   
    char s = '\t' ;
    size_t kmer_size = m_kmermap->getKmerSize() ;
    std::vector<std::string> headers = {"seq", 
                                        "n_occurence"} ;
    // IPD columns
    for(size_t i=0; i<kmer_size; i++)
    {   char ipd[64] ;
        sprintf(ipd, "ipd_%zu", i+1) ;
        headers.push_back(ipd) ;
    }
    // PWD columns
    for(size_t i=0; i<kmer_size; i++)
    {   char ipd[64] ;
        sprintf(ipd, "pwd_%zu", i+1) ;
        headers.push_back(ipd) ;
    }
    // print headers
   print_vector(stream, headers, s) << std::endl ;

    // print map
    for(auto iter=m_kmermap->begin(); 
        iter!=m_kmermap->end(); 
        iter++)
    {   stream << iter->second.sequence           << s
               << iter->first                     << s ;
        print_vector(stream, iter->second.ipd, s) << s ;
        print_vector(stream, iter->second.pwd, s) << std::endl ;
    }
}
