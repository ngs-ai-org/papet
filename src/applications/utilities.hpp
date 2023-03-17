#ifndef NGSAI_APP_UTILITIES_HPP
#define NGSAI_APP_UTILITIES_HPP

namespace ngsai
{
    namespace app
    {
        /*!
        * \brief Writes a vector in tab-separated format on 
        * the given stream.
        * \param stream the stream to print on.
        * \param v the vector of interest.
        * \param separator a separator to use to separate 
        * each header value.
        * \returns a reference to the stream to which the 
        * data were sent.
        */
        template<class T>
        std::ostream& print_vector(std::ostream& stream,
                                   const std::vector<T>& v,
                                   char separator)
        {
            for(size_t i=0; i<v.size(); i++)
            {   stream << v[i] ;
                if(i == (v.size() - 1))
                {   ; }
                else
                {   stream << separator ; }
            }
            return stream ;
        }

    }  // namespace app
    
}  // namespace ngsai

#endif // NGSAI_APP_UTILITIES_HPP