#ifndef NGSAI_APP_APPLICATIONKINETICSWIG_HPP
#define NGSAI_APP_APPLICATIONKINETICSWIG_HPP

#include <applications/ApplicationInterface.hpp>

#include <string>
#include <vector>


namespace ngsai
{
    namespace app
    {   
        /*!
        * \brief The ApplicationKineticsWig class creates  
        * a standalone applications to create a WIG track
        * from kinetic data.
        */
        class ApplicationKineticsWig : 
            public ngsai::app::ApplicationInterface
        {   
            public:
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
                ApplicationKineticsWig(int argc, 
                                         char** argv) ;

                /*!
                * \brief Destructor.
                */
                virtual 
                ~ApplicationKineticsWig() override ;
                
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
                 * \brief Reads the data and creates the 
                 * 4 tracks.
                 */
                void
                createWigTracks() const ;

            protected:
                /*!
                 * The path to the BED file of interest.
                 */
                std::string m_path_bed ;

                /*!
                 * \brief The list of path to the BAM files 
                 * of interest.
                 */
                std::vector<std::string> m_paths_bam ;

                /*!
                 * \brief The prefix for the result files 
                 * paths that will be created.
                 */
                std::string m_prefix_out ;

                /*!
                 * \brief The size in bp of the windows 
                 * around the BED entries in which a 
                 * kinetic track will be created.
                 */
                size_t m_win_size ;
        } ;
    
    }  // namespace app

}  // namespace ngsai

#endif // NGSAI_APP_APPLICATIONKINETICSWIG_HPP