#ifndef NGSAI_APP_APPLICATIONPREDICT_HPP
#define NGSAI_APP_APPLICATIONPREDICT_HPP

#include <applications/ApplicationInterface.hpp>

#include <string>
#include <vector>
#include <future>
#include <ngsaipp/epigenetics/KineticModel.hpp>
#include <ngsaipp/epigenetics/KineticClassifier.hpp>
#include <ngsaipp/genome/CpGRegion.hpp>

namespace ngsai
{
    namespace app
    {   
        /*!
        * \brief The ApplicationPredict class creates a 
        * standalone applications to compute CpG 
        * methylation prediction.
        */
        class ApplicationPredict : 
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
                ApplicationPredict(int argc, char** argv) ;

                /*!
                * \brief Destructor.
                */
                virtual 
                ~ApplicationPredict() override ;
                
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
                 * \brief Loads the kinetic signal models 
                 * and builds the classifier.
                 * \param path_model_meth the path to the 
                 * file containt the methylated kinetic 
                 * signal model to load.
                 * \param path_model_unmeth the path to the 
                 * file containt the unmethylated kinetic 
                 * signal model to load.
                 * \return an exit code, 
                 * getExitCodeSuccess() if it went well.
                 */
                int
                loadModels(
                    const std::string& path_model_meth,
                    const std::string& path_unmodel_meth) ;
                
                /*!
                 * \brief Loads the content of the BED 
                 * file.
                 * \return an exit code, 
                 * getExitCodeSuccess() if it went well.
                 */ 
                int
                loadBed(const std::string& path_bed) ;
                
                /*!
                 * \brief The prediction routine ran by 
                 * worker threads. The predictions will 
                 * be computed on the CpG slice [from,to).
                 * \param from the index of the 1st CpG in 
                 * the list to compute a prediction for.
                 * \param to the index of the past last CpG 
                 * in the list to compute a prediction for.
                 * \param promise a promise to collect a 
                 * list of methylation probabilities 
                 * corresponding to each CpG in the range 
                 * [from,to). The promise is filled when 
                 * the worker thread is done.
                 */
                void
                predictRoutine(
                    size_t from,
                    size_t to,
                    std::promise<std::list<double>>& promise) 
                    const ;

            protected:
                /*!
                 * \brief the paths to the bam files.
                 */
                std::vector<std::string> m_paths_bam ;
                /*!
                 * \brief the signal classifier.
                 */
                ngsai::KineticClassifier m_classifier ;
                /*!
                 * \brief The list of CpGs to compute 
                 * predictions from.
                 */
                std::list<ngsai::genome::CpGRegion> m_cpgs ;
                /*!
                 * \brief the prior probability of 
                 * methylation.
                 */
                double m_prob_meth ;
                /*!
                 * \brief the number of worker threads
                 */
                size_t m_threads_n ;
        } ;
    }
}
#endif // NGSAI_APP_APPLICATIONINTERFACE_HPP

