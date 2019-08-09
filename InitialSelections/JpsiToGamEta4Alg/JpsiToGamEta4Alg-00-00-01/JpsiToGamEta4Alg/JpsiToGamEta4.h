//change all "JpsiToGamEta4" to "YourScriptName"
#ifndef Physics_Analysis_JpsiToGamEta4_H
#define Physics_Analysis_JpsiToGamEta4_H 

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Matrix/SymMatrix.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Geometry/Point3D.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"
#include "DstEvent/TofHitStatus.h"
#include "EmcRawEvent/EmcDigi.h"
#include "EmcRecEventModel/RecEmcHit.h"
#include "EmcRecEventModel/RecEmcShower.h"
#include "EventModel/EventHeader.h"
#include "EventModel/EventModel.h"
#include "EventModel/Event.h"
#include "EventNavigator/EventNavigator.h"
#include "EvTimeEvent/RecEsTime.h"
#include "EvtRecEvent/EvtRecEvent.h"
#include "EvtRecEvent/EvtRecTrack.h"
#include "GaudiKernel/AlgFactory.h"
#include "GaudiKernel/Algorithm.h"
#include "GaudiKernel/Bootstrap.h"
#include "GaudiKernel/IDataProviderSvc.h"
#include "GaudiKernel/IHistogramSvc.h"
#include "GaudiKernel/IJobOptionsSvc.h"
#include "GaudiKernel/INTupleSvc.h"
#include "GaudiKernel/ISvcLocator.h"
#include "GaudiKernel/MsgStream.h"
//#include "GaudiKernel/NTuple.h"    //No NTuple!
#include "GaudiKernel/PropertyMgr.h"
#include "GaudiKernel/SmartDataPtr.h"
#include "Identifier/Identifier.h"
#include "MdcRecEvent/RecMdcKalTrack.h"
#include "MdcRecEvent/RecMdcTrack.h"
#include "MdcRecEvent/RecMdcHit.h"
#include "ParticleID/ParticleID.h"
//#include "PartPropSvc/PartPropSvc.h"
#include "RootCnvSvc/RootCnvSvc.h"
#include "RootCnvSvc/RootInterface.h"
#include "TBenchmark.h"
#include "VertexFit/Helix.h"
#include "VertexFit/IVertexDbSvc.h"
#include "VertexFit/KalmanKinematicFit.h"
#include "VertexFit/KinematicFit.h"
//#include "VertexFit/ReadBeamParFromDb.h"
#include "VertexFit/SecondVertexFit.h"
#include "VertexFit/VertexFit.h"
#include "VertexFit/WTrackParameter.h"
#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <TFile.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TTree.h>
#include <vector>

class JpsiToGamEta4 : public Algorithm {

        public:
                JpsiToGamEta4(const std::string& name, ISvcLocator* pSvcLocator);
                StatusCode initialize();
                StatusCode execute();
                StatusCode finalize();

        private:
                bool getPrimaryVertex(Hep3Vector &, HepSymMatrix &);
                int selectGoodChargedTrack(Hep3Vector, vector<int> &, vector<double> *);
                int selectGoodPhoton(Hep3Vector, vector<int> &, vector<double> *, vector<double> *);
                int selectGoodPhoton2(Hep3Vector, vector<int> &, vector<double> *, vector<double> *);
//                int selectGoodPhoton(Hep3Vector, vector<int> *, vector<double> *, vector<double> *, TLorentzVector *);
                int identifyPID(vector<int>, vector<int> iTen[], vector<double> *, int PID_0=-1, int PID_1=-1, int PID_2=-1, int PID_3=-1, int PID_4=-1);
                bool doVertexFit(VertexFit*, int, vector<int>, vector<int>);
                bool doVertexFitC(VertexFit*, int, vector<int>, vector<int>);
                bool doSecondVertexFit(SecondVertexFit*, VertexFit*, Hep3Vector, HepSymMatrix);
                void assignMomentumToPhoton(Hep3Vector, vector<TLorentzVector> &, vector<int>, Int_t);
                void assignMomentumToCharged(vector<TLorentzVector> &, vector<int>, Int_t, int);
                void assignMomentumToCharged(vector<TLorentzVector> &, vector<int>, Int_t, int, vector<double> *emcE);
                bool buildPi0(double, double, vector<int> &, vector<TLorentzVector> &);
                bool buildPi0(double, double, int, vector<int> &, vector<TLorentzVector> &);
                bool buildPi0(double, double, vector<int> &, vector<TLorentzVector> &, vector<int>);
                bool buildPi0(double, double, int, vector<int> &, vector<TLorentzVector> &, vector<int>);
                void corgen(HepMatrix &, HepVector &, int );
                void corset(HepSymMatrix &, HepMatrix &, int );
                void calibration(RecMdcKalTrack * , HepVector &, int );

                //IPartPropSvc *p_PartPropSvc;
                //HepPDT::ParticleDataTable* m_particleTable;
	        SmartDataPtr<Event::EventHeader> *_eventHeader;
                SmartDataPtr<EvtRecEvent> *_evtRecEvent;
                SmartDataPtr<EvtRecTrackCol> *_evtRecTrkCol;

                bool m_MCTruth;
                //Declare r0, z0 and cos cut for charged tracks
                Double_t m_vr0cut;
                Double_t m_vz0cut;
                Double_t m_ccoscut;
                Double_t m_ptcut;
                Double_t m_pcut;
                //Declare energy, dphi, dthe cuts for fake gamma's
                Double_t m_isoAngleCut;
                //Declare type of psi. For Jpsi psitype=1, for psip psitype=0
                Int_t m_psiType;

                //Declare name of output file
                std::string m_OutputFileName;
                TFile *saveFile;

                //Define TreeAna here
                TTree *TreeAna;
};

#endif

