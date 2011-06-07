#include <protomol/output/OutputXTCTrajectory.h>
#include <protomol/config/Configuration.h>
#include <protomol/output/OutputCache.h>
#include <protomol/module/MainModule.h>
#include <protomol/base/StringUtilities.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/io/DCDTrajectoryWriter.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/base/Exception.h>

// GROMACS headers
#ifdef HAVE_GROMACS
extern "C" {
#include <gromacs/xtcio.h>
}
#endif

using namespace std;
using namespace ProtoMol::Report;
using namespace ProtoMol;

//____ OutputXTCTrajectory
const string OutputXTCTrajectory::keyword("XTCFile");

OutputXTCTrajectory::OutputXTCTrajectory() :
  Output(), myMinimalImage(false), myFrameOffset(0), myFilename("") {

#if defined (HAVE_GROMACS)
    fxtc = NULL;
#endif
}

OutputXTCTrajectory::OutputXTCTrajectory(const string &filename, int freq,
                                         bool minimal, int frameoffs) :
  Output(freq), //new DCDTrajectoryWriter(filename)),
  myMinimalImage(minimal), myFrameOffset(frameoffs), myFilename(filename) {

#if defined (HAVE_GROMACS)
    fxtc = NULL;
#endif
    
    //flag value
    report << plain << "XTC FrameOffset parameter set to " << myFrameOffset << "." << endr;
    
}

OutputXTCTrajectory::~OutputXTCTrajectory() {
}

void OutputXTCTrajectory::doInitialize() {
  
#if defined (HAVE_GROMACS)
    
  //Get first frame (must exist or error)
  const int firstframe = toInt(app->config["firststep"]);
  
  report << debug(2) << "Firstframe " << firstframe << "." << endr;

  //open file
  //if myFrameOffset is zero default to overwrite data
  //note: now include "firstframe" data.
  if( myFrameOffset == 0){
      //new file
      fxtc = open_xtc(myFilename.c_str(),"w");
  }else{
      //append
      fxtc = open_xtc(myFilename.c_str(),"a");
  }

  //test opened
  if( fxtc == NULL ){
    THROW(string("Can not open '") + myFilename +
          "' for " + getId() + ".");

  }

#else

    THROW("GROMACS XTC format output not available.");

#endif

}

void OutputXTCTrajectory::doRun(int) {

#if defined (HAVE_GROMACS)

  const Vector3DBlock *pos =
    (myMinimalImage ? app->outputCache.minimalPositions() : &app->positions);

    //size
    const unsigned int possize = pos->size();

    //number of atoms
    int natoms = possize;
    
    //step
    int step = app->currentStep;

    //real time
    real time = app->outputCache.time() ;
    
    //defines precision, 1000 is the GROMACS default
    //can be read from TPR file
    real prec = 1000;

    //box containing model
    // The computational box which is stored as a set of three basis vectors,
    // to allow for triclinic PBC. For a rectangular box the box edges are
    // stored on the diagonal of the matrix.
    
    //bounding box
    Vector3D a, b;
    app->topology->getBoundingbox(*pos, a, b);

    //set box
    matrix box = { { a.c[0],a.c[1],a.c[2]}, {b.c[0],b.c[1],b.c[2]},
                    {b.c[0] - a.c[0],b.c[1] - a.c[1],b.c[2] - a.c[2]} };

    //Gromacs xyz data struct
    rvec *x = new rvec[possize];

    //copy data
    for (int i = 0; i < possize; i++) {
        x[i][0] = (*pos)[i][0] * Constant::ANGSTROM_NM;
        x[i][1] = (*pos)[i][1] * Constant::ANGSTROM_NM;
        x[i][2] = (*pos)[i][2] * Constant::ANGSTROM_NM;
    }

    //write to file
    int bOK = write_xtc((t_fileio*)fxtc, natoms, step, time,
                      box, x, prec);

    //delete temp data
    delete[] x;

    //test write OK
    if( bOK != 1 ){

        THROW(string("Could not write ") + getId() + " '" +
              myFilename + "'.");

    }

#endif

}

void OutputXTCTrajectory::doFinalize(int) {

#if defined (HAVE_GROMACS)
    //close XTC file
    close_xtc( (t_fileio*)fxtc );
#endif

}

Output *OutputXTCTrajectory::doMake(const vector<Value> &values) const {
  return new OutputXTCTrajectory(values[0], values[1], values[2], values[3]);
}

void OutputXTCTrajectory::getParameters(vector<Parameter> &parameter) const {
  parameter.push_back
    (Parameter(getId(), Value(myFilename,
                              ConstraintValueType::NotEmpty())));
  parameter.push_back
    (Parameter(keyword + "OutputFreq",
               Value(getOutputFreq(), ConstraintValueType::Positive()) ));
  parameter.push_back
    (Parameter(keyword + "MinimalImage", Value(myMinimalImage),
               Text("whether the coordinates should be transformed to minimal "
                    "image or not")));
  parameter.push_back
    (Parameter(keyword + "FrameOffset", 
                Value(myFrameOffset, ConstraintValueType::NotNegative()), 0 ));
}

bool OutputXTCTrajectory::adjustWithDefaultParameters(
  vector<Value> &values, const Configuration *config) const {
  if (!checkParameterTypes(values)) return false;

  if (config->valid(InputOutputfreq::keyword) && !values[1].valid())
    values[1] = (*config)[InputOutputfreq::keyword];

  if (config->valid(InputMinimalImage::keyword) && !values[2].valid())
    values[2] = (*config)[InputMinimalImage::keyword];
  
  return checkParameters(values);
}

