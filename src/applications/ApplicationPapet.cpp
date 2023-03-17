#include <applications/ApplicationPapet.hpp>

#include <iostream>
#include <boost/program_options.hpp>    // PacBio::BAM::variable_map, options_descriptions
#include <unordered_map>                // std::unordered_map

#include <applications/ApplicationModelKinetic.hpp>
#include <applications/ApplicationModelKineticTxt.hpp>
#include <applications/ApplicationKinetics.hpp>
#include <applications/ApplicationKineticsWig.hpp>
#include <applications/ApplicationKineticsKmer.hpp>
#include <applications/ApplicationModelSequence.hpp>
#include <applications/ApplicationModelSequenceTxt.hpp>
#include <applications/ApplicationPredict.hpp>


std::string recepe = "\n"
                     "\tPAPET VAUDOIS:\n"
                     "\t--------------\n\n"
                     "\tRecepe for 4 persons:\n"
                     "\t- 2 pork cabage saussages\n"
                     "\t- 800g of leak\n"
                     "\t- 500g of firm-flesh potatoes \n"
                     "\t- 1 table spoon of butter\n"
                     "\t- 1 tea spoon of salt\n"
                     "\t- 0.5 dl of white wine (preferably Chasselas)\n"
                     "\t- 2dl of half fat cream\n"
                     "\t- 50g parlsey\n"
                     "\n\n"
                     "\tIMPORTANT DISCLAIMER:\n"
                     "\tPick the saussages with tootpicks at both of their edges!\n"
                     "\tReally do it! Otherwise the saussages will get filled with\n"
                     "\twith cooking water, won't get as tasty and are likely to \n"
                     "\tlitteraly pop or explode whenever you will try to cut them\n"
                     "\tReally! whathever you are told about this, do it!!\n"
                     "\n"
                     "\tSimmer the saussages in boiling water for 40 minutes.\n"
                     "\t\n"
                     "\tSlice the leaks in 1cm wide slices and wash the slices to\n"
                     "\tremove any earth. Peel the potatoes and chope them in ~2cm\n"
                     "\twide cubes\n"
                     "\t\n"
                     "\tIn a pot, melt the butter, add in the leak slices and saute\n"
                     "\tthem for 5 minutes. Add the salt, and pepper if you like it.\n"
                     "\tPut the potaotes in and cover the pot with a lid. Let\n"
                     "\teverythingcook for 3 minutes.\n"
                     "\t\n"
                     "\tAdd the wine in the pot, close the lid again let everything\n"
                     "\tcook for 30 minutes at least to 1h. Add some water in the pot\n"
                     "\tif needed such that the leak and potatoes don't stick and burn\n"
                     "\tin the pot.\n"
                     "\t\n"
                     "\tWhen everything is cooked, chope the parsley and add it inside\n"
                     "\tthe pot. Add the cream and 2 to 4 table spoons of red wine\n"
                     "\tvinegar. Mix everything gently. Taste and add salt, pepper,\n"
                     "\tcream or vinegar as you like. Excess of vinegar can be counter\n"
                     "\tbalanced with extra cream\n"
                     "\t\n"
                     "\tBon appetit !\n" ;

namespace po = boost::program_options ;


ngsai::app::ApplicationPapet::
                ApplicationPapet(int argc,
                                 char** argv)
    : ApplicationInterface(argc, argv),
      m_app_map{{app_types::undefined, "undefined"},
                {app_types::model_kinetic,"model-kinetic"},
                {app_types::model_kinetic_txt, "model-kinetic-txt"},
                {app_types::kinetics, "kinetics"},
                {app_types::kinetics_wig, "kinetics-wig"},
                {app_types::kinetics_kmer, "kinetics-kmer"},
                {app_types::model_sequence, "model-sequence"},
                {app_types::model_sequence_txt, "model-sequence-txt"},
                {app_types::predict, "predict"}},
      m_app_cmd(),
      m_app(nullptr)
{   int parsing = this->parseOptions() ;
    if(parsing == this->getExitCodeSuccess())
    {   m_is_runnable = true ; }
    else
    {   m_is_runnable = false ; }
}


ngsai::app::ApplicationPapet::
                ~ApplicationPapet()
{   this->freeApp() ; }


int
ngsai::app::ApplicationPapet::run()
{   
    if(not this->isRunnable())
    {   return this->getExitCodeError() ; }

    return m_app->run() ;
}


int
ngsai::app::ApplicationPapet::parseOptions()
{   // in case
    this->freeApp() ;

    // check arguments were given
    if(m_argc == 1)
    {   std::cerr << "Error! no options given"
                  << std::endl ; 
        return this->getExitCodeError() ;
    }

    // help messages
    char desc_msg[4096] ;
    sprintf(desc_msg,
            "\n"
            "Usage : papet [command] [options]\n"
            "\n"
            "\tpapet is a suite of softwares to model and the\n"
            "\tPacBio kinetic signal and predict epigenetics\n"
            "\tDNA modifications.\n\n"
            "\tpapet contains the commands:\n\n"
            "\t-h --help            Displays this help message\n"
            "\t   --vaudois         Surprise\n\n"
            "\t%s        Creates kinetic signal models from CCSs\n\n"
            "\t%s    Dumps a kinetic signal model in txt format\n\n"
            "\t%s             Extracts CCS kinetic information in txt format.\n\n"
            "\t%s         Creates WIG tracks from CCSs\n\n"
            "\t%s        Computes the per-kmer distribution of kinetic signal from CCSs\n\n"
            "\t%s       Creates DNA sequence kinetic signal models from CCSs\n\n"
            "\t%s   Dumps a DNA sequence kinetic signal model in txt format\n\n"
            "\t%s              Predicts the presence of epignetic modifications from CCSs\n\n"
            "\tWritten by Romain Groux, November 2022\n\n",
            m_app_map.at(app_types::model_kinetic).c_str(),
            m_app_map.at(app_types::model_kinetic_txt).c_str(),
            m_app_map.at(app_types::kinetics).c_str(),
            m_app_map.at(app_types::kinetics_wig).c_str(),
            m_app_map.at(app_types::kinetics_kmer).c_str(),
            m_app_map.at(app_types::model_sequence).c_str(),
            m_app_map.at(app_types::model_sequence_txt).c_str(),
            m_app_map.at(app_types::predict).c_str()) ;

    // check if help invoked
    if((std::string(m_argv[1]) == "-h") or 
       (std::string(m_argv[1]) == "--help"))
    {   std::cout << desc_msg << std::endl ; 
        return this->getExitCodeError() ;
    }

   if((std::string(m_argv[1]) == "--vaudois"))
    {   std::cout << recepe << std::endl ; 
        return this->getExitCodeError() ;
    }

    // get command invoked and update command line
    std::string cmd(m_argv[1]) ;
    for(int i=1; i<m_argc-1; i++)
    {   m_argv[i] = m_argv[i+1] ; }
    m_argc-- ;

    if(cmd == m_app_map.at(app_types::model_kinetic))
    {   m_app_cmd  = cmd ; 
        m_app = 
            new ngsai::app::ApplicationModelKinetic(
                                        m_argc, m_argv) ;
    }
    else if(cmd == m_app_map.at(
                app_types::model_kinetic_txt))
    {   m_app_cmd  = cmd ;
        m_app = new ngsai::app::ApplicationModelKineticTxt(m_argc, m_argv) ;
    }
    else if(cmd == m_app_map.at(app_types::kinetics))
    {   m_app_cmd  = cmd ; 
        m_app = 
            new ngsai::app::ApplicationKinetics(
                                        m_argc, m_argv) ;
    }
    else if(cmd == m_app_map.at(app_types::kinetics_wig))
    {   m_app_cmd  = cmd ; 
        m_app = 
            new ngsai::app::ApplicationKineticsWig(
                                        m_argc, m_argv) ;
    }
    else if(cmd == m_app_map.at(app_types::kinetics_kmer))
    {   m_app_cmd  = cmd ; 
        m_app = 
            new ngsai::app::ApplicationKineticsKmer(
                                        m_argc, m_argv) ;
    }
    else if(cmd == m_app_map.at(app_types::model_sequence))
    {   m_app_cmd  = cmd ; 
        m_app = 
            new ngsai::app::ApplicationModelSequence(
                                        m_argc, m_argv) ;
    }
    else if(cmd == m_app_map.at(
                app_types::model_sequence_txt))
    {   m_app_cmd  = cmd ; 
        m_app = 
            new ngsai::app::ApplicationModelSequenceTxt(
                                        m_argc, m_argv) ;
    }
    else if(cmd == m_app_map.at(app_types::predict))
    {   m_app_cmd  = cmd ; 
        m_app = 
            new ngsai::app::ApplicationPredict(
                                        m_argc, m_argv) ;
    }
    else
    {
        std::cerr << "Error! invalid command : " 
                  << cmd 
                  << std::endl ;
        std::cout << desc_msg << std::endl ; 
        return this->getExitCodeError() ;
    }

    return this->getExitCodeSuccess() ;
}


int
ngsai::app::ApplicationPapet::freeApp()
{
    if(m_app != nullptr)
    {   delete m_app ; 
        m_app = nullptr ;
    }
    return this->getExitCodeSuccess() ;
}


int main(int argc, char** argv)
{
    ngsai::app::ApplicationPapet app(argc, argv) ;
    return app.run() ;
}

