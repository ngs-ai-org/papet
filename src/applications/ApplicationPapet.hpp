#ifndef NGSAI_APP_APPLICATIONPAPET_HPP
#define NGSAI_APP_APPLICATIONPAPET_HPP


#include <applications/ApplicationInterface.hpp>

#include <string>
#include <unordered_map>


namespace ngsai
{
    namespace app
    {   
        /*!
        * \brief The ApplicationPapet class creates  
        * a standalone applications to run PAPET.
        */
        class ApplicationPapet : 
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
                ApplicationPapet(int argc, 
                                 char** argv) ;

                /*!
                * \brief Destructor.
                */
                virtual 
                ~ApplicationPapet() override ;
                
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
                 * \brief Allocates the application to 
                 * be run.
                 * \return an exit code, 
                 * getExitCodeSuccess() if it went well.
                 */
                int
                allocateApp() ;

                /*!
                 * \brief Frees the application to 
                 * be run.
                 * \return an exit code, 
                 * getExitCodeSuccess() if it went well.
                 */
                int
                freeApp() ;
            
            protected:
                /*!
                 * \brief An enumeration indicating the 
                 * the application to be run.
                 */
                enum class app_types {undefined,
                                      model_kinetic,
                                      model_kinetic_txt,
                                      kinetics,
                                      kinetics_wig,
                                      kinetics_kmer,
                                      model_sequence,
                                      model_sequence_txt,
                                      predict} ;

            protected:
                /*!
                 * \brief Maps the app type with the 
                 * command name to call the corresponding 
                 * app.
                 */
                const std::unordered_map<app_types,
                                         std::string> 
                                            m_app_map ;
                
                /*!
                 * \brief The command to run.
                 */
                std::string m_app_cmd ;

                /*!
                 * \brief A pointer to the app to run.
                 */
                ngsai::app::ApplicationInterface* m_app ;
        } ;
    
    }  // namespace app

}  // namespace ngsai


#endif  // NGSAI_APP_APPLICATIONPAPET_HPP