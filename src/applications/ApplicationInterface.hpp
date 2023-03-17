#ifndef NGSAI_APP_APPLICATIONINTERFACE_HPP
#define NGSAI_APP_APPLICATIONINTERFACE_HPP

#include <string>


namespace ngsai
{
    namespace app
    {   
        /*!
        * \brief The ApplicationInterface class provides an 
        * interface for classes deisnged to create 
        * standalone applications which execution only 
        * requires to run a run() method.
        */
        class ApplicationInterface
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
                ApplicationInterface(int argc, 
                                     char** argv) ;

                /*!
                * \brief Destructor.
                */
                virtual
                ~ApplicationInterface() ;
                
                /*!
                * \brief Runs the application, with all its
                * functionalities.
                * \return the exit code to return to the OS.
                */
                virtual
                int
                run() = 0 ;

                /*!
                 * \brief Returns the success exit code.
                 * \return the success exit code.
                 */
                int
                getExitCodeSuccess() const ;

                /*!
                 * \brief Returns the error exit code.
                 * \return the error exit code.
                 */
                int
                getExitCodeError() const ;

                /*!
                 * \brief Indicates whether the app can 
                 * be run.
                 */
                bool
                isRunnable() const ;

            protected:
                /*!
                 * \brief Parses the command line options 
                 * and sets the fields.
                 * \return an exit code, 
                 * getExitCodeSuccess() if it went well.
                 */
                virtual
                int
                parseOptions() = 0;

                /*!
                 * \brief Checks that a given BAM file 
                 * has both a PacBio and a BAI index file 
                 * names <path>.pb and <path>.pbi
                 * \return an exit code, 
                 * getExitCodeSuccess() if it went well.
                 */
                int 
                checkBamFile(
                    const std::string& path) const ;

            protected:
                /*!
                * \brief the number of command line 
                argument.
                */
                int m_argc ;
                /*!
                * \brief the command line arguments.
                */
                char** m_argv ;
                /*!
                * \brief whether the app can run. 
                */
                bool m_is_runnable ;

        } ; // ApplicationInterface

    } // namespace app

} // namespace ngsai

#endif // NGSAI_APP_APPLICATIONINTERFACE_HPP

