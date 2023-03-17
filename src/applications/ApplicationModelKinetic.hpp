#ifndef NGSAI_APP_APPLICATIONMODELKINETIC_HPP
#define NGSAI_APP_APPLICATIONMODELKINETIC_HPP

#include <applications/ApplicationInterface.hpp>

#include <iostream>
#include <vector>

#include <ngsaipp/io/BedRecord.hpp>              // ngsai::BedRecord
#include <ngsaipp/epigenetics/KmerMap.hpp>       // ngsai::KmerMap
#include <ngsaipp/epigenetics/KineticModel.hpp>  // ngsai::KineticModel


namespace ngsai
{
    namespace app
    {   
        /*!
        * \brief The ApplicationModelKinetic class creates  
        * a standalone applications to train a kinetic 
        * signal model from CCSs of interest.
        */
        class ApplicationModelKinetic : 
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
                ApplicationModelKinetic(int argc, 
                                         char** argv) ;

                /*!
                * \brief Destructor.
                */
                virtual 
                ~ApplicationModelKinetic() override ;
                
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
                 * \brief Loads the content of the BED 
                 * file.
                 * \return an exit code, 
                 * getExitCodeSuccess() if it went well.
                 */ 
                int
                loadBed(const std::string& path_bed) ;

                /*!
                 * \brief Allocates the KineticModel memory 
                 * and sets the model parameters.
                 * \return an exit code, 
                 * getExitCodeSuccess() if it went well.
                 */
                int
                allocateKineticModels() ;

                /*!
                 * \brief Frees the KineticModel memory.
                 * \return an exit code, 
                 * getExitCodeSuccess() if it went well.
                 */
                int
                freeKineticModels() ;

            protected:
                /*!
                 * \brief An enumeration indicating the 
                 * exact type of model trained. 
                 */
                enum class modes {undefined,
                                  raw,
                                  raw_norm,
                                  diposition,
                                  diposition_norm,
                                  pairwise,
                                  pairwise_norm} ;
            protected:
                /*!
                 * \brief The type of model that needs to 
                 * be trained.
                 */
                modes m_mode ;
                /*!
                 * \brief the paths to the training data.
                 */
                std::vector<std::string> m_paths_bam ;
                /*!
                 * \brief the path to the file in which the 
                 * KineticModel will be serialized.
                 */
                std::string m_path_out ;
                /*!
                 * \brief the path to the file containing 
                 * the background model to use for 
                 * normalization.
                 */
                std::string m_path_kmermap ;
                /*!
                 * \brief the size of the model in bp.
                 */
                size_t m_size ;
                /*!
                 * \brief the number of bins in the 
                 * histograms.
                 */
                size_t m_nb_bins ;
                /*!
                 * \brief the lower limit of the lower bin
                 * in the histograms.
                 */
                double m_xmin ;
                /*!
                 * \brief the upper limit of the upper bin
                 * in the histograms.
                 */
                double m_xmax ;
                /*!
                 * \brief a pseudo count to add to each 
                 * histogram bin.
                 */
                double m_pseudo_counts ;
                /*!
                 * \brief the number of worker threads.
                 */
                size_t m_nb_threads ;
                /*!
                 * \brief the CpGs from which the training 
                 * should be performed. 
                 */
                std::vector<ngsai::BedRecord> m_cpgs ;
                /*!
                 * \brief the background model to use 
                 * for normalization.
                 */
                ngsai::KmerMap* m_kmermap ;
                /*!
                 * \brief the partial models that will be
                 * trained by each thread.
                 */
                std::vector<ngsai::KineticModel*> m_models ;
        } ;
    
    }  // namespace app

} // namespace ngsai

#endif // NGSAI_APP_APPLICATIONMODELKINETIC_HPP