#ifndef NGSAI_APP_APPLICATIONKINETICS_HPP
#define NGSAI_APP_APPLICATIONKINETICS_HPP

#include <applications/ApplicationInterface.hpp>

#include <string>
#include <vector>
#include <pbbam/BamFile.h>              // BamFile

#include <ngsaipp/epigenetics/KmerMap.hpp>   // ngsai::KmerMap


namespace ngsai
{
    namespace app
    {   
        /*!
        * \brief The ApplicationKinetics class creates a 
        * standalone applications to extract kinetics from 
        * CCSs around genomic coordinates of interest.
        */
        class ApplicationKinetics : 
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
                ApplicationKinetics(int argc, char** argv) ;

                /*!
                * \brief Destructor.
                */
                virtual 
                ~ApplicationKinetics() override ;
                
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
                 * \brief Generates and prints the headers 
                 * for the features that will be extracted.
                 * \param stream the stream on which to 
                 * print the headers.
                 * \param separator a separator to use to 
                 * separate each header value.
                 */
                void printHeader(std::ostream& stream,
                                 char separator) const ;
            
            protected:
                /*!
                 * \brief The BAM files containing the CCSs 
                 * from which the kinetic must be extracted.
                 */
                std::vector<PacBio::BAM::BamFile> m_bam_files ;
                /*!
                 * \brief The path to the BED file 
                 * containing the regions from which the 
                 * CCS kinetic must be extracted.
                 */
                std::string m_path_bed ;
                /*!
                 * \brief The size of the window in which 
                 * the kinetics will be extracted.
                 */
                size_t m_win_size ;
                /*!
                 * \brief A pointer to the KmerMap to use 
                 * to normalize the kinetic signal.
                 */
                ngsai::KmerMap* m_kmermap ;
        } ;

    }  // namespace app

}  // namespace ngsai

#endif // NGSAI_APP_APPLICATIONKINETICS_HPP