#include <applications/ApplicationModelKineticTxt.hpp>

#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <boost/program_options.hpp>        // variable_map, options_descriptions

#include <ngsaipp/epigenetics/RawKineticModel.hpp>                  // ngsai::RawKineticModel
#include <ngsaipp/epigenetics/NormalizedKineticModel.hpp>           // ngsai::NormalizedKineticModel
#include <ngsaipp/epigenetics/PairWiseKineticModel.hpp>             // ngsai::PairWiseKineticModel
#include <ngsaipp/epigenetics/PairWiseNormalizedKineticModel.hpp>   // ngsai::PairWiseNormalizedKineticModel
#include <ngsaipp/epigenetics/DiPositionKineticModel.hpp>           // ngsai::DiPositionKineticModel
#include <ngsaipp/epigenetics/DiPositionNormalizedKineticModel.hpp> // ngsai::DiPositionNormalizedKineticModel
#include <ngsaipp/utility/string_utility.hpp>                       // ngsai::split(), ngsai::endswith()


namespace po = boost::program_options ; 


ngsai::app::ApplicationModelKineticTxt::
                ApplicationModelKineticTxt(
                                    int argc,
                                    char** argv)
    : ApplicationInterface(argc, argv),
      m_model(nullptr)
{   int parsing = this->parseOptions() ;
    if(parsing == this->getExitCodeSuccess())
    {   m_is_runnable = true ; }
    else
    {   m_is_runnable = false ; }
}


ngsai::app::ApplicationModelKineticTxt::
                ~ApplicationModelKineticTxt()
{
    if(m_model != nullptr)
    {   delete m_model ;
        m_model = nullptr ;
    }
}


int
ngsai::app::ApplicationModelKineticTxt::run()
{   
    if(not this->isRunnable())
    {   return this->getExitCodeError() ; }

    // print model
    std::cout << m_model->toString() << std::endl ;
    
    // free memory
    if(m_model != nullptr)
    {   delete m_model ; 
        m_model = nullptr ;
    }

    return this->getExitCodeSuccess() ;
}


int 
ngsai::app::ApplicationModelKineticTxt::parseOptions()
{   
    // check arguments were given
    if(m_argc == 1)
    {   std::cerr << "Error! no options given"
                  << std::endl ; 
        return this->getExitCodeError() ;
    }
    
    // help messages
    std::string desc_msg =  "\n"
                            "Usage : model-kinetic-txt [options] > [FILE]"
                            "\n"
                            "\tConverts a kineticsignal  model file in tsv \n"
                            "\tformat and prints the results on stdout. If the \n"
                            "\tkinetic model is a normalized one, the KmerMap is \n"
                            "\tnot included in the tsv conversion.\n"
                            "\tWritten by Romain Groux, November 2022\n\n" ;
    std::string opt_help_msg   = "Produces this help message." ;
    std::string opt_model_msg  = "The path to the file containing the kinetic "
                                 "model to convert in tsv format." ;


    // option parser
    std::string path_model("") ;

    po::variables_map vm ;
    po::options_description desc(desc_msg) ;
    desc.add_options()
        ("help,h",  opt_help_msg.c_str())
        ("model",   po::value<std::string>(&(path_model)), 
                    opt_model_msg.c_str()) ;
    
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
    {   std::cerr << "Error! an unknown error occured while "
                      "parsing the options" 
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
    if(path_model == "")
    {   std::cerr <<"no model file given (--model)"
                  << std::endl ;
        return this->getExitCodeError() ;
    }

    // load model
    if(this->loadKineticModel(path_model) != 
            this->getExitCodeSuccess())
    {   return this->getExitCodeError() ; }

    return this->getExitCodeSuccess() ;
}


int
ngsai::app::ApplicationModelKineticTxt::loadKineticModel(
                                    const std::string& path)
{   // ensure no model is loaded
    if(m_model != nullptr)
    {   delete m_model ; 
        m_model = nullptr ;
    }

    if(ngsai::endswith(path, ".rawkineticmodel"))
    {   m_model = new ngsai::RawKineticModel() ;
        m_model->load(path) ;
    }
    else if(ngsai::endswith(path, 
                    ".normalizedkineticmodel"))
    {   m_model = new ngsai::NormalizedKineticModel() ;
        m_model->load(path) ;
    }
    else if(ngsai::endswith(path,
                    ".pairwisekineticmodel"))
    {   m_model = new ngsai::PairWiseKineticModel() ;
        m_model->load(path) ;
    }
    else if(ngsai::endswith(path,
                    ".pairwisenormalizedkineticmodel"))
    {   m_model = 
            new ngsai::PairWiseNormalizedKineticModel() ;
        m_model->load(path) ;
    }
    else if(ngsai::endswith(path, 
                    ".dipositionkineticmodel"))
    {   m_model = new ngsai::DiPositionKineticModel() ;
        m_model->load(path) ;
    }
    else if(ngsai::endswith(path, 
                    ".dipositionnormalizedkineticmodel"))
    {   m_model = 
            new ngsai::DiPositionNormalizedKineticModel() ;
        m_model->load(path) ;
    }
    else
    {   std::cerr << "Error! Could not load the "
                  << "KineticModel in "
                  << path
                  << "because could not assert its type"
                  << std::endl ;
        return this->getExitCodeError() ;
    }
    return this->getExitCodeSuccess() ;
}
