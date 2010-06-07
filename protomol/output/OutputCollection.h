/*  -*- c++ -*-  */
#ifndef OUTPUTCOLLECTION_H
#define OUTPUTCOLLECTION_H

#include <list>

namespace ProtoMol {
  class Output;
  class OutputFactory;
  class ProtoMolApp;
  //____ OutputCollection

  /**
     Container class for Output objects invoked at application level.
   */
  class OutputCollection  {
    friend class OutputFactory;

    typedef std::list<Output *> Container;
    typedef std::list<Output *>::iterator iterator;

  public:
    typedef std::list<Output *>::const_iterator const_iterator;
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // Constructors, destructors, assignment
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    OutputCollection();
    ~OutputCollection();
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // New methods of class OutputCollection
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  public:
    /// Initialize all Output objects
    void initialize(const ProtoMolApp *app);

    /// Invoke all Output objects with run().  Returns true if an Output ran.
    bool run(int step);

    /// Finalize all Outout objects
    void finalize(int step);

    /// Add new Output object to the collection
    int getNext() const;

    void adoptOutput(Output *output);

    /// Iterators, const
    const_iterator begin() const {return myOutputList.begin();}
    const_iterator end()   const {return myOutputList.end();}

  private:
    /// Iterators
    iterator begin()       {return myOutputList.begin();}
    iterator end()         {return myOutputList.end();}

    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    // My data members
    //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  private:
    Container myOutputList;

    const ProtoMolApp *app;
  };
}
#endif
