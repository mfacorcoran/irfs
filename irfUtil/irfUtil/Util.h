/**
 * @file Util.h
 * @brief Utility functions used by response function classes and others.
 * @author J. Chiang
 *
 * $Header$
 */

#ifndef irfUtil_Util_h
#define irfUtil_Util_h

#include <string>
#include <vector>

namespace irfUtil {

/**
 * @class Util
 *
 * $Header$
 */

class Util {

public:

   virtual ~Util() {}

   static void getCaldbFile(const std::string &detName, 
                            const std::string &respName,
                            const std::string &version,
                            std::string &filename, long &extnum,
                            const std::string & telescope = "GLAST",
                            const std::string & instrument = "LAT",
                            const std::string & filter = "NONE",
                            const std::string & date = "2003-01-01",
                            const std::string & time = "00:00:00");

protected:

   Util() {}

};

} // namespace irfUtil

#endif // irfUtil_Util_h