#ifndef NGSAI_APP_APPLICATIONMODELSEQUENCETXT_HPP
#define NGSAI_APP_APPLICATIONMODELSEQUENCETXT_HPP

#include <applications/ApplicationInterface.hpp>

#include <string>
#include <vector>

#include <ngsaipp/epigenetics/KmerMap.hpp>   // ngsai::KmerMap

namespace ngsai
{
    namespace app
    {   
        /*!
        * \brief The ApplicationModelSequenceTxt class 
        * creates  a standalone applications to convert a 
        * sequence kinetic model into a structured text 
        * format.
        */
        class ApplicationModelSequenceTxt : 
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
                ApplicationModelSequenceTxt(int argc, 
                                         char** argv) ;

                /*!
                * \brief Destructor.
                */
                virtual 
                ~ApplicationModelSequenceTxt() override ;
                
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
                 * \brief Loads the KmerMap.
                 * \param path the path to the file 
                 * containing the KmerMap.
                 * \return an exit code, 
                 * getExitCodeSuccess() if it went well.
                 */
                int
                loadKmerMap(const std::string& path) ;

                /*!
                 * \brief Prints to KmerMap in text format 
                 * on the given stream.
                 * \param stream the stream of interest.
                 */
                void
                printKmerMap(std::ostream& stream) const ;
            
            protected:
                /*!
                 * \brief A pointer to the KmerMap to 
                 * convert.
                 */
                ngsai::KmerMap* m_kmermap ;
        } ;

    }  // namespace app

} // namespace ngsai

#endif // NGSAI_APP_APPLICATIONMODELSEQUENCETXT_HPP