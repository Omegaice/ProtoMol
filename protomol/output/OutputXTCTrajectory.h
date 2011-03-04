/*  -*- c++ -*-  */
#ifndef OUTPUTXTCTRAJECTORY_H
#define OUTPUTXTCTRAJECTORY_H

#include <protomol/output/Output.h>

namespace ProtoMol {
  class DCDTrajectoryWriter;

  //____ OutputXTCTrajectory
  class OutputXTCTrajectory : public Output {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OutputXTCTrajectory();
    OutputXTCTrajectory(const std::string &filename, int freq,
                            bool minimal, int frameoffs);
    virtual ~OutputXTCTrajectory();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class Output
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  From class Output
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual Output *doMake(const std::vector<Value> &values) const;
    virtual void doInitialize();
    virtual void doRun(int step);
    virtual void doFinalize(int step);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {return keyword;}
    virtual void getParameters(std::vector<Parameter> &parameter) const;
    virtual bool adjustWithDefaultParameters(std::vector<Value> &values,
                                             const Configuration *config) const;

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

  private:

    void *fxtc;//t_fileio
    
    bool myMinimalImage;
    int myFrameOffset;
    string myFilename;
    
  };
}
#endif
