/*  -*- c++ -*-  */
#ifndef OUTPUTFILE_H
#define OUTPUTFILE_H

#include <protomol/output/Output.h>
#include <protomol/base/Exception.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/base/Timer.h>

#ifdef HAVE_LIBFAH
#include <cbang/os/File.h>

#else
#include <fstream>
#endif

namespace ProtoMol {
  class Configuration;

  //____ OutputFile

  /**
     OutputFile is a base Output class to write and append data (text)
     to a  file. It provides size and frequency based and caching and
     also whether the file should be closed and opened between two
     writes. OutputFile implements doRun(), but provides an
     equivalent, cached method doRunCached(). The cached (buffered)
     data is written to file if the number of calls of doRun() is at
     least the cache frequency and also the buffer size is
     satisfied. The open/close of the file depends on the give time
     interval between two writes. Furthermore,     getParameters()
     defines the file name, output frequency, cache frequency, cache
     size,     and the minimal close time keyword.
   */
  class OutputFile : public Output {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OutputFile();
    OutputFile(const std::string &filename, int freq, int cacheFreq,
               int cacheSize, Real closeTime);
    virtual ~OutputFile();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OutputFile
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    /// open and clear file
    bool open();

    /// close file
    bool close();

    /// new implementation of run method
    virtual void doRunCached(int step) = 0;

    /// optional implemenatation for cache flush
    virtual void doFlushCache() {};
  private:
    /// open file for append
    bool reopen();

    /// close file if requested
    bool reclose();

    /// flush cache (buffer) to file
    void flushCache();

    /// clear cache (buffer)
    void clearCache();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Output
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
  private:
    virtual void doRun(int step);
    virtual void doFinalize(int step);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual void getParameters(std::vector<Parameter> &parameter) const;
    virtual bool adjustWithDefaultParameters(std::vector<Value> &values,
                                             const Configuration *config) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  protected:
    std::string myFilename;
    std::stringstream myBuffer;
  private:
    int myCacheFreq;
    int myCount;
    int myCacheSize;
    Real myCloseTime;
    Timer myTimer;

#ifdef HAVE_LIBFAH
    cb::File myFile;
#else
    std::ofstream myFile;
#endif
  };
}
#endif
