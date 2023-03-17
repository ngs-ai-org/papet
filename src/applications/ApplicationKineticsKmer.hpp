#ifndef NGSAI_APP_APPLICATIONKINETICSKMER_HPP
#define NGSAI_APP_APPLICATIONKINETICSKMER_HPP

#include <applications/ApplicationInterface.hpp>

#include <string>
#include <vector>


namespace ngsai
{
    namespace app
    {   
        /*!
        * \brief The ApplicationKineticsKmer class creates  
        * a standalone applications to compute per the per
        * kmer average IPD and PWD kinetic signal.
        */
        class ApplicationKineticsKmer : 
            public ngsai::app::ApplicationInterface
        {   
            public:
                /*!
                 * \brief Generate a vector of headers for 
                 * the table that will be printed.
                 * \param max_ipd the maximum possible IPD 
                 * value.
                 * \param max_pwd the maximum possible PWD 
                 * value. 
                 * \returns a vector of header that can be 
                 * printed.
                 */
                static
                std::vector<std::string> 
                generateHeaders(size_t max_ipd, 
                                size_t max_pwd) ;
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
                ApplicationKineticsKmer(int argc, 
                                         char** argv) ;

                /*!
                * \brief Destructor.
                */
                virtual 
                ~ApplicationKineticsKmer() override ;
                
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

            protected:
                /*!
                 * \brief The list of path to the BAM files 
                 * of interest.
                 */
                std::vector<std::string> m_paths_bam ;

                /*!
                 * \brief The kmer size in bp.
                 */
                size_t m_kmer_size ;
        } ;
    }  // namespace app

}  // namespace ngsai

#endif  // NGSAI_APP_APPLICATIONKINETICSKMER_HPP