#ifndef NGSAI_APP_APPLICATIONMODELSEQUENCE_HPP
#define NGSAI_APP_APPLICATIONMODELSEQUENCE_HPP

#include <applications/ApplicationInterface.hpp>

#include <string>
#include <vector>

#include <ngsaipp/epigenetics/KmerMap.hpp>   // ngsai::KmerMap

namespace ngsai
{
    namespace app
    {   
        /*!
        * \brief The ApplicationModelSequence class creates  
        * a standalone applications to train a sequence 
        * kinetic model from CCSs of interest.
        */
        class ApplicationModelSequence : 
            public ngsai::app::ApplicationInterface
        {   
            public:
                /*!
                * \brief Constructor.
                * Saves the argc and argv values and sets 
                * the app as not runnable.
                * \param argc the number of command line 
                * argument.
                * \param argv the command line argument 
                * vector.
                */
                ApplicationModelSequence(int argc, 
                                         char** argv) ;

                /*!
                * \brief Destructor.
                */
                virtual 
                ~ApplicationModelSequence() override ;
                
                /*!
                * \brief Runs the application, with all its
                * functionalities.
                * \return the exit code to return to the OS.
                */
                virtual 
                int
                run() override ;
            
            protected:
                /*!
                 * \brief Parses the command line options 
                 * and sets the fields.
                 * \return an exit code, 
                 * getExitCodeSuccess() if it went well.
                 */
                virtual
                int
                parseOptions() override ;

                /*!
                 * \brief Updates the KmerMap with the 
                 * CCSs contained in the given file.
                 * \param path_bam the path to a BAM file 
                 * to use.
                 * \return an exit code, 
                 * getExitCodeSuccess() if it went well.
                 */
                int
                updateKmerMap(const std::string& path_bam) ;
            
            protected:
                /*!
                 * \brief The paths to the BAM files 
                 * containing the CCSs of interest.
                 */
                std::vector<std::string> m_paths_bam ;
                /*!
                 * \brief The path to the file in which 
                 * the KmerMap will be dumped.
                 */
                std::string m_path_out ;
                /*!
                 * \brief A pointer to the KmerMap to 
                 * construct.
                 */
                ngsai::KmerMap* m_kmermap ;

        } ;

    }  // namespace app

}  // namespace ngsai

#endif // NGSAI_APP_APPLICATIONMODELSEQUENCE_HPP