//change all "JpsiToGamEta" to "YourScriptName"
#include "JpsiToGamEtaAlg/JpsiToGamEta.h"
//#include "GammaConv/GammaConv.h"
#include "GammaConv/GammaConv.h"
#ifndef ENABLE_BACKWARDS_COMPATIBILITY
typedef HepGeom::Point3D<double> HepPoint3D;
#endif
typedef std::vector<int> Vint;
typedef std::vector<TLorentzVector> Vp4;
using CLHEP::Hep3Vector;
using CLHEP::Hep2Vector;
using CLHEP::HepLorentzVector;
using namespace std;

const double me = 0.00051099891, mmu = 0.105658367, mpi = 0.13957, mpi0 = 0.1349766, mk = 0.493677, mp = 0.93827203, mk0 = 0.497611, meta = 0.547853, momega = 0.78265;
const double xmass[5] = {0.000511, 0.105658, 0.139570, 0.493677, 0.938272};
const double Ejpsi=3.097, Epsip=3.686, Econt=3.650;

Int_t dataCutFlow2Trk, dataCutFlowPID, dataCutFlowGam, dataCutFlowVtx, dataCutFlow2Vtx, dataCutFlowLenLenerr, dataCutFlowPi0, dataCutFlow4C;
Int_t runid, evtid;
double myHelix[6];
//--------------MC Truth-------------------------
Int_t indexmc, pdgid[100], motheridx[100];
TClonesArray *mcParticle;
//--------------Charged Tracks-------------------
Int_t nCharge, nGoodCharge, netGoodCharge, nEMuPiKP[10];//nPip, nPim, nKp, nKm, nProtonp, nProtonm, nElectronp, nElectronm, nMuonp, nMuonm;
vector<double> *chargeCut;
//--------------PID------------------------------
vector<double> *pidCut;
//--------------EMC Showers----------------------
Int_t nShower, nGamma;
vector<double> *showerCut, *showerPos;
//--------------Before 4C------------------------
vector<int> iGamma, iGood, iEMuPiKP[10];
vector<TLorentzVector> pGamma, pEMuPiKP[10];
TClonesArray *ClonesArray[11];
Double_t chisqVtx, chisqKmfit;
TClonesArray *recoilEE;
TClonesArray *cGamma;
TClonesArray *allCharged;
TLorentzVector *Total;
vector<double> *deltaXY, *deltaZ, *angleEE, *psiPair, *cosEG, *eEP, *Rxy, *Rx, *Ry;
vector<int> *idxEE;

//enum DstMdcKalTrack::PidType  null = -1, electron = 0, muon = 1, pion = 2, kaon = 3, proton = 4
TClonesArray *&Electronp = ClonesArray[0], *&Electronm = ClonesArray[1],
             *&Muonp = ClonesArray[2],     *&Muonm = ClonesArray[3],
             *&Pip = ClonesArray[4],       *&Pim = ClonesArray[5],
             *&Kp = ClonesArray[6],        *&Km = ClonesArray[7],
             *&Protonp = ClonesArray[8],   *&Protonm = ClonesArray[9],
             *&Gamma = ClonesArray[10];

Int_t &nelectronp = nEMuPiKP[0], &nelectronm = nEMuPiKP[1],
      &nmuonp = nEMuPiKP[2],     &nmuonm = nEMuPiKP[3],
      &npip = nEMuPiKP[4],       &npim = nEMuPiKP[5],
      &nkp = nEMuPiKP[6],        &nkm = nEMuPiKP[7],
      &nprotonp = nEMuPiKP[8],   &nprotonm = nEMuPiKP[9];

Vint  &ielectronp = iEMuPiKP[0], &ielectronm = iEMuPiKP[1],
      &imuonp = iEMuPiKP[2],     &imuonm = iEMuPiKP[3],
      &ipip = iEMuPiKP[4],       &ipim = iEMuPiKP[5],
      &ikp = iEMuPiKP[6],        &ikm = iEMuPiKP[7],
      &iprotonp = iEMuPiKP[8],   &iprotonm = iEMuPiKP[9];

Vp4   &pelectronp = pEMuPiKP[0], &pelectronm = pEMuPiKP[1],
      &pmuonp = pEMuPiKP[2],     &pmuonm = pEMuPiKP[3],
      &ppip = pEMuPiKP[4],       &ppim = pEMuPiKP[5],
      &pkp = pEMuPiKP[6],        &pkm = pEMuPiKP[7],
      &pprotonp = pEMuPiKP[8],   &pprotonm = pEMuPiKP[9];
/////////////////////////////////////////////////////////////////////////////

JpsiToGamEta::JpsiToGamEta(const std::string& name, ISvcLocator* pSvcLocator) : Algorithm(name, pSvcLocator){
   //Declare the properties
   declareProperty("OutputFileName",  m_OutputFileName = "JpsiToGamEtatest.root");
   declareProperty("PsiType",m_psiType = 1);
   declareProperty("Vr0cut", m_vr0cut=2000.0);
   declareProperty("Vz0cut", m_vz0cut=30.0);
   declareProperty("Ccoscut", m_ccoscut=0.93);
   declareProperty("IsoAngleCut", m_isoAngleCut=10.0);
   declareProperty("MCTruth", m_MCTruth=true);
}


StatusCode JpsiToGamEta::initialize(){
   MsgStream log(msgSvc(), name());
   log << MSG::INFO << "in initialize()" << endmsg;
   StatusCode status;

   dataCutFlow2Trk=0; dataCutFlowPID=0; dataCutFlowGam=0; dataCutFlowVtx=0; dataCutFlow2Vtx=0; dataCutFlowLenLenerr=0; dataCutFlowPi0=0; dataCutFlow4C=0;

   //-***********Initialize the output structure**************
   TString s_OutputFileName(m_OutputFileName);
   s_OutputFileName.ReplaceAll("[\"","");
   s_OutputFileName.ReplaceAll("\"]","");
   saveFile = new TFile(s_OutputFileName, "recreate");

   //-***********Initialize the Analysis Tree*******************
   deltaXY=new vector<double>();
   deltaZ=new vector<double>();
   angleEE=new vector<double>();
   psiPair=new vector<double>();
   cosEG=new vector<double>();
   eEP=new vector<double>();
   Rxy=new vector<double>();
   Rx=new vector<double>();
   Ry=new vector<double>();
   idxEE=new vector<int>();

   recoilEE   = new TClonesArray("TLorentzVector");
   cGamma   = new TClonesArray("TLorentzVector");
   allCharged   = new TClonesArray("TLorentzVector");
   for(int i=0;i<11;i++) ClonesArray[i] = new TClonesArray("TLorentzVector");
   Total   = new TLorentzVector();

   TreeAna = new TTree("TreeAna", "analysis");


   TreeAna->Branch("runid", &runid, "runid/I");
   TreeAna->Branch("evtid", &evtid, "evtid/I");
   TreeAna->Branch("nGoodCharge", &nGoodCharge, "nGoodCharge/I");
   TreeAna->Branch("nGamma", &nGamma, "nGamma/I");
   TreeAna->Branch("nEMuPiKP", &nEMuPiKP, "nEMuPiKP[10]/I");

//   TreeAna->Branch("chargeCut", "vector<double>", &chargeCut);
//   TreeAna->Branch("pidCut", "vector<double>", &pidCut);
//   TreeAna->Branch("showerCut", "vector<double>", &showerCut);
//   TreeAna->Branch("showerPos", "vector<double>", &showerPos);

   TreeAna->Branch("Gamma", "TClonesArray", &Gamma, 256000, 0);
   //TreeAna->Branch("Ep", "TClonesArray", &Electronp, 256000, 0);
   //TreeAna->Branch("Em", "TClonesArray", &Electronm, 256000, 0);
   TreeAna->Branch("recoilEE","TClonesArray",&recoilEE,256000,0);
   TreeAna->Branch("cGamma","TClonesArray",&cGamma,256000,0);
   TreeAna->Branch("eEP", "vector<double>", &eEP);
   TreeAna->Branch("allCharged","TClonesArray",&allCharged,256000,0);
   TreeAna->Branch("Total", &Total, 32000,0);
   TreeAna->Branch("deltaXY", "vector<double>", &deltaXY);
   TreeAna->Branch("deltaZ", "vector<double>", &deltaZ);
   TreeAna->Branch("angleEE", "vector<double>", &angleEE);
   TreeAna->Branch("psiPair", "vector<double>", &psiPair);
   TreeAna->Branch("cosEG", "vector<double>", &cosEG);
   TreeAna->Branch("Rxy", "vector<double>", &Rxy);
   TreeAna->Branch("Rx", "vector<double>", &Rx);
   TreeAna->Branch("Ry", "vector<double>", &Ry);
   TreeAna->Branch("idxEE", "vector<int>", &idxEE);
   TreeAna->Branch("myHelix", myHelix,"myHelix[6]/D");

   for(int i=0;i<11;i++) ClonesArray[i]->BypassStreamer();
   recoilEE->BypassStreamer();
   cGamma->BypassStreamer();
   allCharged->BypassStreamer();

   if(m_MCTruth){
      TreeAna->Branch("indexmc",&indexmc,"indexmc/I");
      TreeAna->Branch("pdgid", pdgid,"pdgID[indexmc]/I");
      TreeAna->Branch("motheridx",motheridx,"motheridx[indexmc]/I");

      mcParticle   = new TClonesArray("TLorentzVector");
      TreeAna->Branch("mcParticle","TClonesArray",&mcParticle,256000,0);
      mcParticle->BypassStreamer();
   }
   //-********************************************************

   //--------end of book--------
   log << MSG::INFO << "successfully return from initialize()" <<endmsg;
   return StatusCode::SUCCESS;

}



StatusCode JpsiToGamEta::execute() {

   //-************Initialize Global Variables**********************
   runid = -1; evtid = -1;
   nCharge = 0; nGoodCharge = 0; netGoodCharge = 0; nShower = 0; nGamma = 0;
   for(int i=0;i<10;i++) nEMuPiKP[i]=0;
   iGood.clear(); iGamma.clear();
   for(int i=0;i<10;i++) (iEMuPiKP[i]).clear();
   pGamma.clear();
   for(int i=0;i<10;i++) (pEMuPiKP[i]).clear();

//   chargeCut->clear(); pidCut->clear(); showerCut->clear(); showerPos->clear();
   for(int i=0;i<11;i++) ClonesArray[i]->Clear();
   recoilEE->Clear();
   cGamma->Clear();
   allCharged->Clear();
   Total->SetPxPyPzE(0,0,0,0);
   deltaXY->clear(); deltaZ->clear(); angleEE->clear(); psiPair->clear(); cosEG->clear(); eEP->clear(); Rxy->clear(); Rx->clear(); Ry->clear(); idxEE->clear();
   //-************************************************

   //-******** Do not change this section ***************************************************-//--
   MsgStream log(msgSvc(), name());                                                         //--
   log << MSG::INFO << "in execute()" << endreq;                                            //--
                                                                                            //--
   SmartDataPtr<Event::EventHeader> eventHeader(eventSvc(),"/Event/EventHeader");           //--
   _eventHeader = &eventHeader;                                                             //--
   int runNo=eventHeader->runNumber(), event=eventHeader->eventNumber();                    //--
   log << MSG::DEBUG <<"run, evtnum = "<< runNo << " , " << event <<endreq;                 //--
                                                                                            //--
   SmartDataPtr<EvtRecEvent> evtRecEvent(eventSvc(),EventModel::EvtRec::EvtRecEvent);       //--
   _evtRecEvent = &evtRecEvent;                                                             //--
   log << MSG::DEBUG <<"ncharg, nneu, tottks = " << evtRecEvent->totalCharged() << " , "    //--
       << evtRecEvent->totalNeutral() << " , " << evtRecEvent->totalTracks() <<endreq;      //--
                                                                                            //--
   SmartDataPtr<EvtRecTrackCol> evtRecTrkCol(eventSvc(),EventModel::EvtRec::EvtRecTrackCol);//--
   _evtRecTrkCol = &evtRecTrkCol;                                                           //--
                                                                                            //--
   if(event%1000 == 0) cout<<"Processing "<<event<<"th event..."<<endl;                     //--
                                                                                            //--
   //Global Event Parameters(do not change below)                                           //--
   double ecms, ESpread;                                                                    //--
   int psi;                                                                                 //--
   if(m_psiType == 0){ ecms=Epsip; psi=100443; ESpread=0.0013; }//psi                       //--
   if(m_psiType == 1){ ecms=Ejpsi; psi=443; ESpread=0.0008; }   //Jpsi                      //--
   if(m_psiType == 2){ ecms=Econt; psi=100443; ESpread=0.0013; }//continume                 //--
                                                                                            //--
   HepLorentzVector cms(0.011*ecms, 0., 0., ecms);    //for kinematic fit                   //--
   TLorentzVector p_cms(0.011*ecms, 0., 0., ecms);                                          //--
   //-***************************************************************************************-//--

   //-*****************Primary Vertex*****************
   Hep3Vector xorigin(0,0,0);
   HepSymMatrix VtxErr(3,0);
   if(!getPrimaryVertex(xorigin,VtxErr)) return StatusCode::SUCCESS; // get information of average IP in each run
   //-*****************************************************

   //-***********Good Charged Track Selection**************
   nCharge = evtRecEvent->totalCharged();
   netGoodCharge = selectGoodChargedTrack(xorigin,iGood,chargeCut);
   nGoodCharge = iGood.size();
   if(nGoodCharge<2) return StatusCode::SUCCESS;
dataCutFlow2Trk++;
   //-*********** Finish Good Charged Track Selection ***********

   //-********************PID identification*****************
   int netChargePID = identifyPID(iGood,iEMuPiKP,pidCut);
   for(int i=0;i<10;i++) nEMuPiKP[i] = (iEMuPiKP[i]).size();
   if(nelectronp<1||nelectronm<1) return StatusCode::SUCCESS;
dataCutFlowPID++;
   assignMomentumToCharged(pEMuPiKP[0], iEMuPiKP[0], nEMuPiKP[0], 0, eEP);
   assignMomentumToCharged(pEMuPiKP[1], iEMuPiKP[1], nEMuPiKP[1], 0, eEP);
   //for(int i=0;i<pelectronp.size();i++) new ((*Electronp)[i]) TLorentzVector(pelectronp[i]);
   //for(int i=0;i<pelectronm.size();i++) new ((*Electronm)[i]) TLorentzVector(pelectronm[i]);
   for(int i=2;i<10;i++) assignMomentumToCharged(pEMuPiKP[i], iEMuPiKP[i], nEMuPiKP[i], i/2);
   //-********************End PID identification*****************

   //-********************Gamma Conversion*****************
   int iPair=0;
   HepPoint3D IPG(xorigin[0],xorigin[1],xorigin[2]);
   for(int iep=0;iep<nelectronp;iep++){
      for(int iem=0;iem<nelectronm;iem++){
double mass=(p_cms-pelectronp[iep]-pelectronm[iem]).M();
//if(mass<0||mass>1.5||(*eEP)[iep]<0.8||(*eEP)[iem+nelectronp]<0.8) continue;
if(mass<0||mass>1.5) continue;
         if((pelectronp[iep]+pelectronm[iem]).M()>0.15) continue;
         RecMdcKalTrack *epTrk =(*((*_evtRecTrkCol)->begin()+ielectronp[iep]))->mdcKalTrack();
         RecMdcKalTrack *emTrk =(*((*_evtRecTrkCol)->begin()+ielectronm[iem]))->mdcKalTrack();
         myHelix[0]=(epTrk->helix())[1];
         myHelix[1]=(epTrk->helix())[2];
         myHelix[2]=(epTrk->helix())[4];
         myHelix[3]=(emTrk->helix())[1];
         myHelix[4]=(emTrk->helix())[2];
         myHelix[5]=(emTrk->helix())[4];
         RecMdcKalTrack::setPidType(RecMdcKalTrack::electron);
         GammaConv gconv = GammaConv(epTrk->helix(),emTrk->helix(),IPG);
         if(gconv.getRXY()>15) continue;
         //if(gconv.getRXY()<1||gconv.getRXY()>10) continue;
         if(gconv.getDeltaXY()<-5||gconv.getDeltaXY()>3) continue;
         if(fabs(gconv.getPsipair())>1.5) continue;
         new ((*cGamma)[iPair]) TLorentzVector(pelectronp[iep]+pelectronm[iem]);
         new ((*recoilEE)[iPair]) TLorentzVector(p_cms-pelectronp[iep]-pelectronm[iem]);

         deltaXY->push_back(gconv.getDeltaXY());
         deltaZ->push_back(gconv.getDeltaZ());
         angleEE->push_back(gconv.getXiep());
         psiPair->push_back(gconv.getPsipair());
         cosEG->push_back(gconv.getCthe());
         Rxy->push_back(gconv.getRXY());
         Rx->push_back(gconv.getRX());
         Ry->push_back(gconv.getRY());
         idxEE->push_back(iep); idxEE->push_back(nelectronp+iem);//index in "allCharged" of positron and electron

         iPair++;
      }
   }
   if(iPair==0) return StatusCode::SUCCESS;
   //-********************end Gamma Conversion*****************

   //-********************Good Photon selection*****************
   nGamma = selectGoodPhoton(xorigin,iGamma,showerCut,showerPos);
dataCutFlowGam++;
   assignMomentumToPhoton(xorigin, pGamma, iGamma, nGamma);
   for(int i=0;i<pGamma.size();i++){
//if(pGamma[i].E()>1.2) return StatusCode::SUCCESS;
new ((*Gamma)[i]) TLorentzVector(pGamma[i]);
}
   //-********************End Good Photon selection*****************

   Total->SetPxPyPzE(0,0,0,0);
   for(int i=0;i<10;i++) for(int j=0;j<(pEMuPiKP[i]).size();j++) *Total=*Total+(pEMuPiKP[i])[j];
   for(int i=0;i<pGamma.size();i++) *Total=*Total+pGamma[i];

   // ******************** MC Truth *****************
   if(m_MCTruth){

      mcParticle->Clear();
      SmartDataPtr<Event::McParticleCol> mcParticleCol(eventSvc(),"/Event/MC/McParticleCol");
      if (!mcParticleCol) return  (StatusCode::FAILURE);
      else{
         int m_numParticle = 0;
         bool jpsiDecay = false;
         bool m_strange = false;
         int  jpsiIndex = -1;

         Event::McParticleCol::iterator iter_mc = mcParticleCol->begin();
         for(; iter_mc != mcParticleCol->end(); iter_mc++){

            if ((*iter_mc)->primaryParticle()&&(*iter_mc)->particleProperty()==11&&((*iter_mc)->mother()).particleProperty()== 11) m_strange=true;
            if ((*iter_mc)->primaryParticle()) continue;//e+e-
            if (!(*iter_mc)->decayFromGenerator()) continue;
            if ((*iter_mc)->particleProperty()== 443){
              jpsiDecay = true;
              jpsiIndex = (*iter_mc)->trackIndex();
            }
            if (!jpsiDecay) continue;

            int m_mcidx = ((*iter_mc)->mother()).trackIndex() - jpsiIndex;
            int m_pdgid = (*iter_mc)->particleProperty();
            if(m_strange&&((*iter_mc)->mother()).particleProperty()!= 443) m_mcidx--;

            pdgid[m_numParticle] = m_pdgid;
            motheridx[m_numParticle] = m_mcidx;
            HepLorentzVector mcParticleH = (*iter_mc)->initialFourMomentum();
            TLorentzVector mcParticleT(mcParticleH.px(),mcParticleH.py(),mcParticleH.pz(),mcParticleH.e());
            new ((*mcParticle)[m_numParticle]) TLorentzVector(mcParticleT);
            m_numParticle ++;
         }//end for mcParticleCol
         indexmc= m_numParticle;
      }//end else

   }//end if m_MCTruth
   // ******************** End MC Truth *****************

   int num=0;
   for(int i=0;i<10;i++)
      for(int j=0;j<(pEMuPiKP[i]).size();j++){
         new ((*allCharged)[num]) TLorentzVector((pEMuPiKP[i])[j]);
         num++;
      }

   runid = runNo; evtid = event;
//if((nShower*6)!=showerCut->size()) cout<<runNo<<","<<event<<": nShower*6!=showerCut.size()"<<endl;
   TreeAna->Fill();
   return StatusCode::SUCCESS;
}


StatusCode JpsiToGamEta::finalize() {
   cout << "In finalize()..." << endl;

   cout<<"dataCutFlowwwww:"<<dataCutFlow2Trk<<", "<<dataCutFlowPID<<", "<<dataCutFlowGam<<", "<<dataCutFlowVtx<<", "<<dataCutFlow2Vtx<<", "<<dataCutFlowLenLenerr<<","<<dataCutFlowPi0<<", "<<dataCutFlow4C<<endl;

   saveFile->cd();
   TreeAna->Write();
   saveFile->Close();

   return StatusCode::SUCCESS;
}


////////////////////////////////////////////////////////////////////
///////////////////// MEMBER FUNCTIONS /////////////////////////////
////////////////////////////////////////////////////////////////////

//-----------Get Primary Vertex-------------
bool JpsiToGamEta::getPrimaryVertex(Hep3Vector &xorigin, HepSymMatrix &VtxErr){
   IVertexDbSvc*  vtxsvc;
   Gaudi::svcLocator()->service("VertexDbSvc", vtxsvc);
   if(vtxsvc->isVertexValid()){
      double* dbv = vtxsvc->PrimaryVertex(); 
      double*  vv = vtxsvc->SigmaPrimaryVertex();  
      xorigin.setX(dbv[0]);
      xorigin.setY(dbv[1]);
      xorigin.setZ(dbv[2]);
      VtxErr[0][0] = vv[0]*vv[0];
      VtxErr[1][1] = vv[1]*vv[1];
      VtxErr[2][2] = vv[2]*vv[2];
      return true;
   }
   else return false;// if cannot load vertex information, will go to another event
}
//-----------End Get Primary Vertex-------------

//-----------Good Charged Track Selection-------------
int JpsiToGamEta::selectGoodChargedTrack(Hep3Vector xorigin, vector<int> &iGood, vector<double> *chargeCut){
   int netCharge = 0;
   for(int i = 0; i < (*_evtRecEvent)->totalCharged(); i++){
      if(i >= (*_evtRecTrkCol)->size()) break;

      EvtRecTrackIterator itTrk=(*_evtRecTrkCol)->begin() + i;
      if(!(*itTrk)->isMdcTrackValid()) continue;
      if(!(*itTrk)->isMdcKalTrackValid()) continue;

      RecMdcTrack *mdcTrk = (*itTrk)->mdcTrack();
      double pch=mdcTrk->p();
      double ptch=mdcTrk->pxy();
      double x0=mdcTrk->x();
      double y0=mdcTrk->y();
      double z0=mdcTrk->z();
      double phi0=mdcTrk->helix(1);
      double xv=xorigin.x();
      double yv=xorigin.y();
      double Rxy=(x0-xv)*cos(phi0)+(y0-yv)*sin(phi0);

      HepVector a = mdcTrk->helix();
      HepSymMatrix Ea = mdcTrk->err();
      HepPoint3D point0(0.,0.,0.);   // the initial point for MDC recosntruction
      HepPoint3D IP(xorigin[0],xorigin[1],xorigin[2]); 
      VFHelix helixip(point0,a,Ea); 
      helixip.pivot(IP);
      HepVector vecipa = helixip.a();
      double  Rvxy0=fabs(vecipa[0]);  //the nearest distance to IP in xy plane
      double  Rvz0=vecipa[3];         //the nearest distance to IP in z direction
      double  Rvphi0=vecipa[1];
      double ccos=cos(mdcTrk->theta());

      RecMdcKalTrack *mdcKalTrk = (*itTrk)->mdcKalTrack();
//      chargeCut->push_back(mdcKalTrk->charge());
//      chargeCut->push_back(ccos);
//      chargeCut->push_back(Rvz0);
//      chargeCut->push_back(Rvxy0);

      if(fabs(ccos) > m_ccoscut) continue;
      if(fabs(Rvz0) >= m_vz0cut ) continue;
      if(fabs(Rvxy0) >= m_vr0cut) continue;

      if(mdcKalTrk->charge() == 0) continue;

      iGood.push_back((*itTrk)->trackId());
      netCharge += mdcKalTrk->charge();
   }
   return netCharge;
}
//-----------End Good Charged Track Selection-------------

//-----------Good Photon Selection-------------
int JpsiToGamEta::selectGoodPhoton(Hep3Vector xorigin, vector<int> &iGamma, vector<double> *showerCut, vector<double> *showerPos){
   for(int i = (*_evtRecEvent)->totalCharged(); i< (*_evtRecEvent)->totalTracks(); i++){
      if(i >= (*_evtRecTrkCol)->size()) break;

      EvtRecTrackIterator itTrk=(*_evtRecTrkCol)->begin() + i;
      if(!(*itTrk)->isEmcShowerValid()) continue;
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z());
      //emcpos = emcpos - xorigin;//Added by XQYuan***********************************

      // find the nearest charged track  
      double dthe = 200.;
      double dphi = 200.;
      double dang = 200.; 
      for(int j = 0; j < (*_evtRecEvent)->totalCharged(); j++){
         if(j >= (*_evtRecTrkCol)->size()) break;

         EvtRecTrackIterator jtTrk = (*_evtRecTrkCol)->begin() + j;
         if(!(*jtTrk)->isExtTrackValid()) continue;
         RecExtTrack *extTrk = (*jtTrk)->extTrack();
         if(extTrk->emcVolumeNumber() == -1) continue;
         Hep3Vector extpos = extTrk->emcPosition();// - xorigin;

         double angd = extpos.angle(emcpos);
         double thed = extpos.theta() - emcpos.theta();
         double phid = extpos.deltaPhi(emcpos);
         thed = fmod(thed+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;
         phid = fmod(phid+CLHEP::twopi+CLHEP::twopi+pi, CLHEP::twopi) - CLHEP::pi;

         if(angd < dang){
            dang = angd;
            dthe = thed;
            dphi = phid;
         }
      }
//      if(dang >= 200) continue;
      double eraw = emcTrk->energy();
      dthe = dthe * 180 / (CLHEP::pi);
      dphi = dphi * 180 / (CLHEP::pi);
      dang = dang * 180 / (CLHEP::pi);

      // good photon cut will be set here
      double the = emcpos.theta();
      double e_threshold = 10.0;
      if(fabs(cos(the)) < 0.8)   e_threshold = 0.025;
      else if((fabs(cos(the)) > 0.86) && (fabs(cos(the)) < 0.92)) e_threshold = 0.050;

//      showerCut->push_back(eraw); showerCut->push_back(dAngKlong); showerCut->push_back(dang); showerCut->push_back(the); showerCut->push_back(emcpos.phi()); showerCut->push_back(emcTrk->time());
//      showerPos->push_back(emcTrk->x()); showerPos->push_back(emcTrk->y()); showerPos->push_back(emcTrk->z());
      if(eraw < e_threshold) continue;
      if(fabs(dang) < m_isoAngleCut) continue;
      if(emcTrk->time()>14  || emcTrk->time()<0) continue;
      iGamma.push_back((*itTrk)->trackId());
   }
   return iGamma.size();
}
//-----------End Good Photon Selection-------------

//-----------PID identification-------------
int JpsiToGamEta::identifyPID(vector<int> iGood, vector<int> iEMuPiKP[], vector<double> *pidCut, int PID_0, int PID_1, int PID_2, int PID_3, int PID_4){
   int netChargePID = 0;
   ParticleID *pid = ParticleID::instance();
   for(int i = 0; i < nGoodCharge; i++) {
      EvtRecTrackIterator itTrk = (*_evtRecTrkCol)->begin() + iGood[i];
      //if(pid) delete pid;
      pid->init();
      pid->setMethod(pid->methodProbability());
      //pid->setMethod(pid->methodLikelihood());  //for Likelihood Method  

      pid->setChiMinCut(4);//????????????????????
      pid->setRecTrack(*itTrk);
      pid->usePidSys(pid->useDedx() | pid->useTof1() | pid->useTof2() | pid->useTofE() | pid->useEmc() | pid->useMuc()); // use PID sub-system
      pid->identify(pid->all()|pid->onlyMuon());//pid->identify(pid->onlyPion() | pid->onlyKaon() | pid->onlyMuon() | pid->onlyProton() | pid->onlyElectron());

      pid->calculate();
      if(!(pid->IsPidInfoValid())) continue;

      Double_t prob_pid[5];
      prob_pid[0] = pid->probElectron(); //prob_e
      prob_pid[1] = pid->probMuon(); //prob_mu
      prob_pid[2] = pid->probPion(); //prob_pi
      prob_pid[3] = pid->probKaon(); //prob_k
      prob_pid[4] = pid->probProton(); //prob_p
      int maxPID = 0;
      for(int iPID=1;iPID<5;iPID++){ if(prob_pid[maxPID]<=prob_pid[iPID]) maxPID = iPID;}
//      pidCut->push_back(maxPID); pidCut->push_back(prob_pid[maxPID]);

      if((maxPID==PID_0||maxPID==PID_1||maxPID==PID_2||maxPID==PID_3||maxPID==PID_4||PID_0+PID_1+PID_2+PID_3+PID_4==-5) && prob_pid[maxPID]>0.001){
         RecMdcKalTrack* mdcKalTrk = (*itTrk)->mdcKalTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
         RecMdcKalTrack::setPidType  (RecMdcKalTrack::PidType(maxPID));//enum DstMdcKalTrack::PidType  null = -1, electron = 0, muon = 1, pion = 2, kaon = 3, proton = 4
         if(mdcKalTrk->charge() > 0) {
            (iEMuPiKP[maxPID*2]).push_back(iGood[i]);
            netChargePID ++;
         }
         else {
            (iEMuPiKP[maxPID*2+1]).push_back(iGood[i]);
            netChargePID --;
         }
      }
   }
   return netChargePID;
}
//-----------End PID identification-------------

//-----------Vertex Fit-------------
bool JpsiToGamEta::doVertexFit(VertexFit* vtxfit, int nTrk, vector<int> trkID, vector<int> trkType){

   HepPoint3D vx(0., 0., 0.);
   HepSymMatrix Evx(3, 0);
   double bx = 1E+6;
   double by = 1E+6;
   double bz = 1E+6;
   Evx[0][0] = bx*bx;
   Evx[1][1] = by*by;
   Evx[2][2] = bz*bz;
   VertexParameter vxpar;
   vxpar.setVx(vx);
   vxpar.setEvx(Evx);

   vtxfit->init();
   WTrackParameter *wvTrk = new WTrackParameter[nTrk];
   RecMdcKalTrack **RTrk = new RecMdcKalTrack*[nTrk];
   vector<int> list;
   for(int i_trk=0;i_trk<nTrk;i_trk++){
      RTrk[i_trk] = (*((*_evtRecTrkCol)->begin()+trkID[i_trk]))->mdcKalTrack();
      RTrk[i_trk]->setPidType(RecMdcKalTrack::PidType(trkType[i_trk]));
      wvTrk[i_trk] = WTrackParameter(xmass[trkType[i_trk]], RTrk[i_trk]->getZHelix(), RTrk[i_trk]->getZError());
      vtxfit->AddTrack(i_trk, wvTrk[i_trk]);
      list.push_back(i_trk);
   }
   vtxfit->AddVertex(0, vxpar, list);
   if(vtxfit->Fit()){
      vtxfit->Swim(0);
      delete []wvTrk;
      delete []RTrk;
      return true;
   }
   else{
      delete []wvTrk;
      delete []RTrk;
      return false;
   }
}
//-----------End Vertex Fit-------------

//-----------Second Vertex Fit-------------
bool JpsiToGamEta::doSecondVertexFit(SecondVertexFit *svtxfit, VertexFit* vtxfit, Hep3Vector xorigin, HepSymMatrix VtxErr){
   VertexParameter vxrawpar;
   vxrawpar.setVx(xorigin);
   vxrawpar.setEvx(VtxErr); 
   svtxfit->init();
   svtxfit->setPrimaryVertex(vxrawpar);
   svtxfit->AddTrack(0, vtxfit->wVirtualTrack(0));
   svtxfit->setVpar(vtxfit->vpar(0));
   if(svtxfit->Fit()) return true;
   else return false;
}
//-----------End Second Vertex Fit-------------

//-----------Assign Momentum to Photons-------------
void JpsiToGamEta::assignMomentumToPhoton(Hep3Vector xorigin, Vp4 &pGamma, Vint iGamma, Int_t nGamma){
   for(int i = 0; i < nGamma; i++){
      RecEmcShower* emcTrk = (*((*_evtRecTrkCol)->begin() + iGamma[i]))->emcShower();
      double eraw = emcTrk->energy();

      Hep3Vector emcpos(emcTrk->x(), emcTrk->y(), emcTrk->z()); 
      Hep3Vector gammaDirection = emcpos;// - xorigin;
      double phi = gammaDirection.phi();
      double the = gammaDirection.theta();
      //if(fabs(gammaDirection.cosTheta())>0.93)continue;
/*
      double phi = emcTrk->phi();
      double the = emcTrk->theta();
*/
      TLorentzVector ptrk;
      ptrk.SetPx(eraw*sin(the)*cos(phi));
      ptrk.SetPy(eraw*sin(the)*sin(phi));
      ptrk.SetPz(eraw*cos(the));
      ptrk.SetE(eraw);
      pGamma.push_back(ptrk);
   }
}
//-----------End Assign Momentum to Photons-------------

//-----------Assign Momentum to Charged Particles-------------
void JpsiToGamEta::assignMomentumToCharged(Vp4 &pCharged, Vint iCharged, Int_t nGoodCharge, int pidType){
   for(int i = 0; i < nGoodCharge; i++){
      RecMdcKalTrack* mdcKalTrk = (*((*_evtRecTrkCol)->begin() + iCharged[i]))->mdcKalTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
      RecMdcKalTrack::setPidType(RecMdcKalTrack::PidType(pidType));//enum DstMdcKalTrack::PidType  null = -1, electron = 0, muon = 1, pion = 2, kaon = 3, proton = 4

      TLorentzVector ptrk;
      ptrk.SetPx(mdcKalTrk->px());
      ptrk.SetPy(mdcKalTrk->py());
      ptrk.SetPz(mdcKalTrk->pz());
      double p3 = ptrk.P();
      ptrk.SetE(sqrt(p3*p3+xmass[pidType]*xmass[pidType]));
      pCharged.push_back(ptrk);
   }
}
//-----------End Assign Momentum to Charged Particles-------------

//-----------Assign Momentum to Charged Particles-------------
void JpsiToGamEta::assignMomentumToCharged(Vp4 &pCharged, Vint iCharged, Int_t nGoodCharge, int pidType, vector<double> *emcE){
   for(int i = 0; i < nGoodCharge; i++){
      RecMdcKalTrack* mdcKalTrk = (*((*_evtRecTrkCol)->begin() + iCharged[i]))->mdcKalTrack();//After ParticleID, use RecMdcKalTrack substitute RecMdcTrack
      RecMdcKalTrack::setPidType(RecMdcKalTrack::PidType(pidType));//enum DstMdcKalTrack::PidType  null = -1, electron = 0, muon = 1, pion = 2, kaon = 3, proton = 4

      TLorentzVector ptrk;
      ptrk.SetPx(mdcKalTrk->px());
      ptrk.SetPy(mdcKalTrk->py());
      ptrk.SetPz(mdcKalTrk->pz());
      double p3 = ptrk.P();
      ptrk.SetE(sqrt(p3*p3+xmass[pidType]*xmass[pidType]));
      pCharged.push_back(ptrk);

      if(!(*((*_evtRecTrkCol)->begin()+iCharged[i]))->isEmcShowerValid()) emcE->push_back(-1);
      else emcE->push_back(((*((*_evtRecTrkCol)->begin()+iCharged[i]))->emcShower())->energy()/p3);
   }
}
//-----------End Assign Momentum to Charged Particles-------------

//-----------Build Particle with Mass Window-------------
bool JpsiToGamEta::buildPi0(double massWin_lower, double massWin_upper, vector<int> &gammaID, vector<TLorentzVector> &pPi0){
//gammaID:gamma在pGamma中的排序，从0开始
   if(pGamma.size() < 2){
      cout<<"ERROR! There're less than 2 gammas in pGamma"<<endl;
      return false;
   }

   bool pi0Built = false;
   TLorentzVector pion0(0,0,0,0);
   for(int i=0;i<pGamma.size()-1;i++){
      for(int i2=i+1;i2<pGamma.size();i2++){
         pion0 = pGamma[i]+pGamma[i2];//pGamma是否加*?
         if(pion0.M() < massWin_lower || pion0.M() > massWin_upper) continue;

         pPi0.push_back(pion0);
         gammaID.push_back(i);
         gammaID.push_back(i2);
         pi0Built = true;
      }
   }
   return pi0Built;
}
//-----------End Build Particle with Mass Window-------------

//-----------Build Particle with Mass Window-------------
bool JpsiToGamEta::buildPi0(double massWin_lower, double massWin_upper, vector<int> &gammaID, vector<TLorentzVector> &pPi0, vector<int> inputGamID){
//gammaID:gamma在pGamma中的排序，从0开始
   if(inputGamID.size() < 2){
      cout<<"ERROR! There're less than 2 gammas in inputGamma"<<endl;
      return false;
   }

   bool pi0Built = false;
   TLorentzVector pion0(0,0,0,0);
   for(int i=0;i<inputGamID.size()-1;i++){
      for(int i2=i+1;i2<inputGamID.size();i2++){
         pion0 = pGamma[inputGamID[i]]+pGamma[inputGamID[i2]];//pGamma是否加*?
         if(pion0.M() < massWin_lower || pion0.M() > massWin_upper) continue;

         pPi0.push_back(pion0);
         gammaID.push_back(inputGamID[i]);
         gammaID.push_back(inputGamID[i2]);
         pi0Built = true;
      }
   }
   return pi0Built;
}
//-----------End Build Particle with Mass Window-------------

//-----------Build Particles with Mass Window-------------
bool JpsiToGamEta::buildPi0(double massWin_lower, double massWin_upper, int nPi0ToBuild, vector<int> &gammaID, vector<TLorentzVector> &pPi0, vector<int> inputGamID){
//gammaID:gamma在pGamma中的排序，从0开始
   if(inputGamID.size() < 2*nPi0ToBuild){
      cout<<"ERROR! There're less than "<<2*nPi0ToBuild<<" gammas in pGamma to build "<<nPi0ToBuild<<" #pi0"<<endl;
      return false;
   }
   if(nPi0ToBuild<=0){
      cout<<"ERROR! Can't build "<<nPi0ToBuild<<" #pi0"<<endl;
      return false;
   }

   if(nPi0ToBuild == 1) return buildPi0(massWin_lower,massWin_upper,gammaID,pPi0,inputGamID);
   bool pi0Built = false;
   TLorentzVector pion0(0,0,0,0);
   for(int i=0;i<inputGamID.size()-nPi0ToBuild*2+1;i++){
      for(int i2=i+1;i2<inputGamID.size();i2++){
         pion0 = pGamma[inputGamID[i]]+pGamma[inputGamID[i2]];
         if(pion0.M() < massWin_lower || pion0.M() > massWin_upper) continue;
         vector<int> subGammaID;
         vector<TLorentzVector> subPi0;
         vector<int> subInputGamID;
         for(int i3=i+1;i3<inputGamID.size();i3++) if(i3!=i&&i3!=i2) subInputGamID.push_back(inputGamID[i3]);
         if(!buildPi0(massWin_lower,massWin_upper,nPi0ToBuild-1,subGammaID,subPi0,subInputGamID)) continue;
         pi0Built = true;
         for(int i3=0;i3<(subPi0.size()/(nPi0ToBuild-1));i3++){
            gammaID.push_back(inputGamID[i]);
            gammaID.push_back(inputGamID[i2]);
            for(int i4=0;i4<2*(nPi0ToBuild-1);i4++) gammaID.push_back(subGammaID[i4+i3*2*(nPi0ToBuild-1)]);
            pPi0.push_back(pion0);
            for(int i4=0;i4<(nPi0ToBuild-1);i4++) pPi0.push_back(subPi0[i4+i3*(nPi0ToBuild-1)]);
         }
      }
   }

   return pi0Built;
}
//-----------End Build Particles with Mass Window-------------

//-----------Build Particles with Mass Window-------------
bool JpsiToGamEta::buildPi0(double massWin_lower, double massWin_upper, int nPi0ToBuild, vector<int> &gammaID, vector<TLorentzVector> &pPi0){
//gammaID:gamma在pGamma中的排序，从0开始
   if(pGamma.size() < 2*nPi0ToBuild){
      cout<<"ERROR! There're less than "<<2*nPi0ToBuild<<" gammas in pGamma to build "<<nPi0ToBuild<<" #pi0"<<endl;
      return false;
   }
   if(nPi0ToBuild<=0){
      cout<<"ERROR! Can't build "<<nPi0ToBuild<<" #pi0"<<endl;
      return false;
   }

   if(nPi0ToBuild == 1) return buildPi0(massWin_lower,massWin_upper,gammaID,pPi0);
   bool pi0Built = false;
   TLorentzVector pion0(0,0,0,0);
   for(int i=0;i<pGamma.size()-2*nPi0ToBuild+1;i++){
      for(int ii=i+1;ii<pGamma.size();ii++){
         pion0 = pGamma[i]+pGamma[ii];
         if(pion0.M() < massWin_lower || pion0.M() > massWin_upper) continue;
         vector<int> subGammaID;
         vector<TLorentzVector> subPi0;
         vector<int> subInputGamID;
         for(int i3=i+1;i3<pGamma.size();i3++) if(i3!=i&&i3!=ii) subInputGamID.push_back(i3);
         if(!buildPi0(massWin_lower,massWin_upper,nPi0ToBuild-1,subGammaID,subPi0,subInputGamID)) continue;
         for(int i3=0;i3<(subPi0.size()/(nPi0ToBuild-1));i3++){
            gammaID.push_back(i);
            gammaID.push_back(ii);
            for(int i4=0;i4<2*(nPi0ToBuild-1);i4++) gammaID.push_back(subGammaID[i4+i3*2*(nPi0ToBuild-1)]);
            pPi0.push_back(pion0);
            for(int i4=0;i4<(nPi0ToBuild-1);i4++) pPi0.push_back(subPi0[i4+i3*(nPi0ToBuild-1)]);
            pi0Built = true;
         }
      }
   }

   return pi0Built;
}
//-----------End Build Particles with Mass Window-------------

