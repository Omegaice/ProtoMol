/*  -*- c++ -*-  */
#ifndef OUTPUTCHECKPIINT_H
#define OUTPUTCHECKPIINT_H

#include <protomol/base/StringUtilities.h>
#include <protomol/output/Output.h>
#include <protomol/base/Timer.h>

namespace ProtoMol {
    class Configuration;

    //____ OutputCheckpoint
    class OutputCheckpoint : public Output {
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // Constructors, destructors, assignment
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        public:
            OutputCheckpoint();
            OutputCheckpoint( const std::string &name, int freq, int start,
                              const std::string& posbase, const std::string& velbase );
        private:
            void WritePositions( int step );
            void WriteVelocities( int step );
            void WriteConfig( int step );

        public:
            void doIt( int step );

            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            //  From class Output
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        private:
            virtual Output *doMake( const std::vector<Value> &values ) const;
            virtual void doInitialize();
            virtual void doRun( int step );
            virtual void doFinalize( int );
            virtual bool isIdDefined( const Configuration *config ) const;
            virtual bool addDoKeyword() const {return false;}

            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // From class Makeable
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        public:
            virtual std::string getIdNoAlias() const {return keyword;}
            virtual void getParameters( std::vector<Parameter> & ) const;
            virtual bool adjustWithDefaultParameters( std::vector<Value> &values,
                    const Configuration *config ) const;

            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // My data members
            //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        public:
            static const std::string keyword;

        private:
            int mCurrent;
            std::string mName;
            std::string mPosBase, mVelBase;
    };
}

#endif
