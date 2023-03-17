#include <applications/ApplicationInterface.hpp>


#include <iostream>
#include <string>
#include <pbbam/BamFile.h>


ngsai::app::ApplicationInterface::ApplicationInterface(
                                int argc,
                                char** argv)
    : m_argc(argc),
      m_argv(argv),
      m_is_runnable(false)
{ ; }


ngsai::app::ApplicationInterface::~ApplicationInterface()
{ ; }


int 
ngsai::app::ApplicationInterface::getExitCodeSuccess() const
{   return 0 ; }


int
ngsai::app::ApplicationInterface::getExitCodeError() const
{   return 1  ;}


bool
ngsai::app::ApplicationInterface::isRunnable() const
{   return m_is_runnable ; }


int
ngsai::app::ApplicationInterface::checkBamFile(
                        const std::string& path) const
{   try
    { 
      PacBio::BAM::BamFile file_bam(path) ;
      if (not file_bam.IsPacBioBAM())
      {   std::cerr << path 
                  << " is not a PacBio BAM file" 
                  << std::endl ;
          return this->getExitCodeError() ;
      }
      if(not file_bam.PacBioIndexExists())
      {   std::cerr << "Error! "
                    << path 
                    << " has no PacBio BAM index file" 
                    << std::endl ;
          return this->getExitCodeError() ;
      }
      else if (not file_bam.StandardIndexExists())
      {   std::cerr << "Error! "
                    << path 
                    << " has no a BAM index file" 
                    << std::endl ;
          return this->getExitCodeError() ;
      }

    }
    catch(const std::exception& e)
    {   std::cerr << "Error! something occured while "
                     "trying to check "
                  << path 
                  << " :" << std::endl
                  << e.what() << std::endl ; 
        return this->getExitCodeError() ;
    }
    return this->getExitCodeSuccess() ;
}