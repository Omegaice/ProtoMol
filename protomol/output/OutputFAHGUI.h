/*  -*- c++ -*-  */
#ifndef OUTPUT_FAH_GUI_H
#define OUTPUT_FAH_GUI_H

//Updated for standalone GUI
#if defined (HAVE_GUI) || defined (HAVE_LIBFAH)

#include <protomol/output/Output.h>
#include <protomol/base/Timer.h>
#include <string>

#ifdef HAVE_LIBFAH
namespace FAH {
  class GUIServer;
}
#else
#include <protomol/output/GUIServer.h>
#endif

namespace ProtoMol {
  class Configuration;

  //____ OutputFAHGUI
  class OutputFAHGUI : public Output {
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OutputFAHGUI();
    OutputFAHGUI(const std::string &name, int freq, int port, 
                 int prange, const std::string &projn, Real timeout,
                 bool pause);

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    //  From class Output
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    virtual Output *doMake(const std::vector<Value> &values) const;
    virtual void doInitialize();
    virtual void doRun(int step);
    virtual void doFinalize(int);
    virtual bool isIdDefined(const Configuration *config) const;
    virtual bool addDoKeyword() const {return false;}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // From class Makeable
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    virtual std::string getIdNoAlias() const {return keyword;}
    virtual void getParameters(std::vector<Parameter> &) const;
    virtual bool adjustWithDefaultParameters(std::vector<Value> &values,
                                             const Configuration *config) const;

  private:
    void setCoords();
    void setBonds();
    void setAtoms();

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    static const std::string keyword;

  private:
    std::string name;
    int myPort, myPortRange;
    std::string myProjName;
    Real myTimeout;
    bool myPause;
    Timer guiTimer;

#ifdef HAVE_LIBFAH
    FAH::GUIServer *server;
#else
    GUIServer *server;
#endif
  };
}

#endif // HAVE_GUI
#endif
