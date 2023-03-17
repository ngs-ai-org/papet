#ifndef NGSAI_APP_APPLICATIONMODELKINETICTXT_HPP
#define NGSAI_APP_APPLICATIONMODELKINETICTXT_HPP

#include <applications/ApplicationInterface.hpp>

#include <string>
#include <ngsaipp/epigenetics/KineticModel.hpp>


namespace ngsai
{
    namespace app
    {
        class ApplicationModelKineticTxt : 
                    public ApplicationInterface
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
                ApplicationModelKineticTxt(int argc, 
                                           char** argv) ;

                /*!
                * \brief Destructor.
                */
                virtual
                ~ApplicationModelKineticTxt() override ;
                
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
                 * \brief Loads the KineticModel stored in 
                 * the given file.
                 * \param path the path to the file to 
                 * load.
                 * \return an exit code, 
                 * getExitCodeSuccess() if it went well.
                 */
                int
                loadKineticModel(const std::string& path) ;

            protected:
                KineticModel* m_model ;
        } ;
    }  // namespace app

}  // namespace ngsai



#endif // NGSAI_APP_APPLICATIONMODELKINETICTXT_HPP