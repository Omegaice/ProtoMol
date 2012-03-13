#include <protomol/force/ForceGroup.h>

#include <protomol/force/system/SystemForce.h>
#include <protomol/type/Vector3DBlock.h>
#include <protomol/topology/GenericTopology.h>
#include <protomol/type/ScalarStructure.h>
#include <protomol/force/extended/ExtendedForce.h>
#include <protomol/force/MollyForce.h>
#include <protomol/force/MetaForce.h>
#include <protomol/base/TimerStatistic.h>
#include <protomol/ProtoMolApp.h>
#include <protomol/parallel/Parallel.h>

//____#define DEBUG_OUTSTANDING_MSG

#ifdef DEBUG_OUTSTANDING_MSG
#include <mpi.h>
#endif

#ifdef HAVE_MPI
#include <mpi.h>
#endif

using namespace std;
using namespace ProtoMol;
using namespace ProtoMol::Report;

//____ ForceGroup
ForceGroup::ForceGroup() {}

ForceGroup::~ForceGroup() {
  for (list<SystemForce *>::iterator currentForce =
         mySystemForcesList.begin(); currentForce != mySystemForcesList.end();
       ++currentForce)
    delete (*currentForce);

  for (list<ExtendedForce *>::iterator currentForce =
         myExtendedForcesList.begin();
       currentForce != myExtendedForcesList.end();
       ++currentForce)
    delete (*currentForce);

  for (list<MollyForce *>::iterator currentForce = myMollyForcesList.begin();
       currentForce != myMollyForcesList.end();
       ++currentForce)
    delete (*currentForce);

  for (list<MetaForce *>::iterator currentForce = myMetaForcesList.begin();
       currentForce != myMetaForcesList.end();
       ++currentForce)
    delete (*currentForce);
}

void ForceGroup::evaluateSystemForces(ProtoMolApp *app,
                                      Vector3DBlock *forces) const {
	if (mySystemForcesList.empty()) return;
	app->topology->uncacheCellList();

  //not parallel forces?~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if( !Parallel::isParallel() ){
    
		TimerStatistic::timer[TimerStatistic::FORCES].start();
		list<SystemForce *>::const_iterator currentForce;
    
		for (currentForce = mySystemForcesList.begin(); currentForce != mySystemForcesList.end(); ++currentForce){
      (*currentForce)->evaluate(app->topology, &app->positions, forces, &app->energies);
    }

    //post process
    for (currentForce = mySystemForcesList.begin(); currentForce != mySystemForcesList.end(); ++currentForce){
      (*currentForce)->postProcess();
    }

    TimerStatistic::timer[TimerStatistic::FORCES].stop();

    //return;
    
  }else{
    //parallel code here~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  
    Parallel::distribute(&app->energies, forces);
    
    //find forces that require post-parallel processing
    list<SystemForce *>::const_iterator startForce = mySystemForcesList.begin();
    
    list<SystemForce *>::const_iterator stopAtForce;
    
    //loop until all forces complete
    while( startForce != mySystemForcesList.end() ){
      
      //int rank;
      //MPI_Comm_rank(MPI_COMM_WORLD,&rank); 
      //std::cout << "Force iter " << (*startForce)->getId() << ", " << rank << std:: endl;
      
      int doPostParallel = 0;
      
      //find end of list OR next post-parallel process
      for (stopAtForce = startForce; stopAtForce != mySystemForcesList.end(); ++stopAtForce) {
      
        if( (*stopAtForce)->getId().find( "BornRadii" ) != std::string::npos ){		
          ++stopAtForce;
          doPostParallel = 1;
          break;
        }
        
        if( (*stopAtForce)->getId().find( "BornSelf" ) != std::string::npos ){		
          ++stopAtForce;
          doPostParallel = 2;
          break;
        }

      }
      
      //
      if (Parallel::isDynamic()) {
        // Collecting the number of blocks of each force.
        vector<int> blocks;
        list<SystemForce *>::const_iterator currentForce;
        for (currentForce = startForce; currentForce != stopAtForce; ++currentForce) {
          blocks.push_back((*currentForce)->numberOfBlocks(app->topology, &app->positions));
        }

        Parallel::resetNext(blocks);
      }

      if (Parallel::iAmSlave()) {
        Parallel::resetNext();

        TimerStatistic::timer[TimerStatistic::FORCES].start();
        list<SystemForce *>::const_iterator currentForce;
        
        for (currentForce = startForce; currentForce != stopAtForce; ++currentForce){

          (*currentForce)->parallelEvaluate(app->topology, &app->positions, forces, &app->energies);
				
				}//do forces
        
        //local reduce here
        if( doPostParallel ){
          
          const unsigned int atomnumber = app->topology->atoms.size();
          
          if( doPostParallel == 1 ){
            // Copy Radii
            Real *radii = new Real[ atomnumber ];

            //put radii (minus zeta) into array
            for( unsigned int i = 0; i < atomnumber; i++ ){
              radii[i] = app->topology->atoms[i].mySCPISM_A->bornRadius - app->topology->atoms[i].mySCPISM_A->zeta;
            }
            
            //sum accross nodes
            Parallel::reduce(radii, radii + atomnumber); //reduceSlaves only?
            
            //put radii back and add in zeta
            for( unsigned int i = 0; i < atomnumber; i++ ){
              app->topology->atoms[i].mySCPISM_A->bornRadius = radii[i] + app->topology->atoms[i].mySCPISM_A->zeta;
            }
            
            delete [] radii;
              
          }else{
          
            //find self energy count
            // Copy self energy count
            Real *selfcount = new Real[ atomnumber ];
            
            //put self energy count into array
            for( unsigned int i = 0; i < atomnumber; i++ ){
              selfcount[i] = (Real)app->topology->atoms[i].mySCPISM_A->energySum;
            }
            
            //sum accross nodes
            Parallel::reduce(selfcount, selfcount + atomnumber); //reduceSlaves only?
            
            //find self energies
            // Copy self energies
            Real *selfenergy = new Real[ atomnumber ];
            
            //put self energy into array
            for( unsigned int i = 0; i < atomnumber; i++ ){
              selfenergy[i] = app->topology->atoms[i].mySCPISM_A->selfEnergy;
            }
            
            //sum accross nodes
            Parallel::reduce(selfenergy, selfenergy + atomnumber); //reduceSlaves only?
            
            //put corrected self energy back
            for( unsigned int i = 0; i < atomnumber; i++ ){
              
              if( selfcount[i] != 0.0 )
                app->topology->atoms[i].mySCPISM_A->selfEnergy = selfenergy[i] / selfcount[i] / (Real)Parallel::getNum();
            }
            
            delete [] selfenergy;
            
            delete [] selfcount;
          }
        }
        
        TimerStatistic::timer[TimerStatistic::FORCES].stop();
        
			}//if slave
      
      //point to next steps
      startForce = stopAtForce;
            
		}//stop while
    
    Parallel::reduce(&app->energies, forces);
    
    //post process
    list<SystemForce *>::const_iterator currentForce;
    
    for (currentForce = mySystemForcesList.begin(); currentForce != mySystemForcesList.end(); ++currentForce){
      (*currentForce)->postProcess();
    }


	}//lel-serial test

}

//____ Evaluate all system forces in this group.

void ForceGroup::evaluateExtendedForces(ProtoMolApp *app, Vector3DBlock *forces) const {
  if (myExtendedForcesList.empty()) return;

  app->topology->uncacheCellList();

  Parallel::distribute(&app->energies, forces);

  if (Parallel::isDynamic()) {
    // Collecting the number of blocks of each force.
    vector<int> blocks;
    list<ExtendedForce *>::const_iterator currentForce;
    for (currentForce = myExtendedForcesList.begin();
         currentForce != myExtendedForcesList.end(); ++currentForce)
      blocks.push_back
        ((*currentForce)->numberOfBlocks(app->topology, &app->positions));

    Parallel::resetNext(blocks);
  }

  if (Parallel::iAmSlave()) {
    Parallel::resetNext();

  TimerStatistic::timer[TimerStatistic::FORCES].start();
  list<ExtendedForce *>::const_iterator currentForce;
  for (currentForce = myExtendedForcesList.begin();
       currentForce != myExtendedForcesList.end(); ++currentForce)
    if (Parallel::isParallel())
      (*currentForce)->parallelEvaluate(app->topology, &app->positions,
                                        &app->velocities, forces,
                                        &app->energies);
    else
      (*currentForce)->evaluate(app->topology, &app->positions,
                                &app->velocities, forces, &app->energies);

  TimerStatistic::timer[TimerStatistic::FORCES].stop();
}

#ifdef DEBUG_OUTSTANDING_MSG
  report << allnodes << plain << "Node " << Parallel::getId() << " done."
         << endr;

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Status status;
  int test = 0;
  MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &test, &status);
  if (test != 0)
    report << plain << allnodes << "Node " << Parallel::getId() <<
      " outstanding msg from " << status.MPI_SOURCE << endr;
#endif

  Parallel::reduce(&app->energies, forces);
}

void ForceGroup::evaluateMollyForces(GenericTopology *topo,
                                     const Vector3DBlock *positions,
                                     vector<ReducedHessAngle> *angleFilter)
const {
  if (myMollyForcesList.empty())
    return;

  topo->uncacheCellList();

  TimerStatistic::timer[TimerStatistic::FORCES].start();

  for (list<MollyForce *>::const_iterator currentForce =
         myMollyForcesList.begin();
       currentForce != myMollyForcesList.end();
       ++currentForce)
    (*currentForce)->evaluate(topo, positions, angleFilter);

  TimerStatistic::timer[TimerStatistic::FORCES].stop();
}

void ForceGroup::addSystemForce(SystemForce *force) {
  if (force != NULL)
    mySystemForcesList.push_back(force);
}

void ForceGroup::addExtendedForce(ExtendedForce *force) {
  if (force != NULL)
    myExtendedForcesList.push_back(force);
}

void ForceGroup::addMollyForce(MollyForce *force) {
  if (force != NULL)
    myMollyForcesList.push_back(force);
}

void ForceGroup::addMetaForce(MetaForce *force) {
  if (force != NULL)
    myMetaForcesList.push_back(force);
}

void ForceGroup::addForce(Force *force) {
  force->addToForceGroup(this);
}

void ForceGroup::getDefinition(vector<MakeableDefinition> &forces) const {
  for (list<SystemForce *>::const_iterator currentForce =
         mySystemForcesList.begin();
       currentForce != mySystemForcesList.end();
       ++currentForce)
    forces.push_back((*currentForce)->getDefinition());

  for (list<ExtendedForce *>::const_iterator currentForce =
         myExtendedForcesList.begin();
       currentForce != myExtendedForcesList.end();
       ++currentForce)
    forces.push_back((*currentForce)->getDefinition());

  for (list<MollyForce *>::const_iterator currentForce =
         myMollyForcesList.begin();
       currentForce != myMollyForcesList.end();
       ++currentForce)
    forces.push_back((*currentForce)->getDefinition());

  for (list<MetaForce *>::const_iterator currentForce =
         myMetaForcesList.begin();
       currentForce != myMetaForcesList.end();
       ++currentForce)
    forces.push_back((*currentForce)->getDefinition());
}

void ForceGroup::uncache() {
  for (list<SystemForce *>::iterator currentForce = mySystemForcesList.begin();
       currentForce != mySystemForcesList.end();
       ++currentForce)
    (*currentForce)->uncache();

  for (list<ExtendedForce *>::iterator currentForce =
         myExtendedForcesList.begin();
       currentForce != myExtendedForcesList.end();
       ++currentForce)
    (*currentForce)->uncache();

  for (list<ExtendedForce *>::iterator currentForce =
         myExtendedForcesList.begin();
       currentForce != myExtendedForcesList.end();
       ++currentForce)
    (*currentForce)->uncache();

  for (list<MetaForce *>::iterator currentForce = myMetaForcesList.begin();
       currentForce != myMetaForcesList.end();
       ++currentForce)
    (*currentForce)->uncache();
}

vector<Force *> ForceGroup::getForces() const {
  vector<Force *> res;
  for (list<SystemForce *>::const_iterator currentForce =
         mySystemForcesList.begin();
       currentForce != mySystemForcesList.end();
       ++currentForce)
    res.push_back(*currentForce);

  for (list<ExtendedForce *>::const_iterator currentForce =
         myExtendedForcesList.begin();
       currentForce != myExtendedForcesList.end();
       ++currentForce)
    res.push_back(*currentForce);

  for (list<ExtendedForce *>::const_iterator currentForce =
         myExtendedForcesList.begin();
       currentForce != myExtendedForcesList.end();
       ++currentForce)
    res.push_back(*currentForce);

  for (list<MetaForce *>::const_iterator currentForce =
         myMetaForcesList.begin();
       currentForce != myMetaForcesList.end();
       ++currentForce)
    res.push_back(*currentForce);

  return res;
}

vector<Force *> ForceGroup::getDeepMetaForces() const {
  vector<Force *> res;
  for (list<MetaForce *>::const_iterator currentForce =
         myMetaForcesList.begin();
       currentForce != myMetaForcesList.end();
       ++currentForce)
    (*currentForce)->getDeepForces(res);

  return res;
}

