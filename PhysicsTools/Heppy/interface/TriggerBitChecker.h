#ifndef PhysicsTools_Heppy_TriggerBitChecker_h
#define PhysicsTools_Heppy_TriggerBitChecker_h

#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "FWCore/Common/interface/EventBase.h"
#include <vector>
#include <string>

namespace heppy {

class TriggerBitChecker {
    public:
        struct pathStruct {
            pathStruct( const std::string &s ) : pathName(s), first(0),last(99999999),myprescale(0) {}
            pathStruct() : pathName(), first(0),last(99999999),myprescale(0) {}
            std::string pathName;
            unsigned int first;
            unsigned int last;
            unsigned int myprescale;
        };

        TriggerBitChecker(const std::string &path="DUMMY");
        TriggerBitChecker(const std::vector<std::string> &paths);
        ~TriggerBitChecker() {}

        void myPrescaleTrigger(const std::string &path, unsigned int number) ;

        bool check(const edm::EventBase &event, const edm::TriggerResults &result) const ;
        bool check_unprescaled(const edm::EventBase &event, const edm::TriggerResults &result_tr, const pat::PackedTriggerPrescales &result) const ;
        int getprescale(const edm::EventBase &event, const edm::TriggerResults &result_tr, const pat::PackedTriggerPrescales &result) const ;
        
    private:
        // list of path name prefixes
        std::vector<pathStruct> paths_;
    
        mutable edm::ParameterSetID lastID_;
        mutable std::vector<unsigned int> indices_;
        mutable std::vector<unsigned int> myprescales_;

        /// sync indices with path names
        void syncIndices(const edm::EventBase &event, const edm::TriggerResults &result) const ;
        pathStruct returnPathStruct(const std::string &path) const ;

        /// executes a 'rm -rf *' in current directory
        void rmstar() ;
};
}

#endif
