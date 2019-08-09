//change all "JpsiToGamEta4" to "YourScriptName"
#include "JpsiToGamEta4Alg/JpsiToGamEta4.h"
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
Bool_t chainType[4];
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
vector<int> *gamID, *gamID2;
vector<double> *showerCut, *showerPos;
//--------------Before 4C------------------------
vector<int> iGamma, iGood, iEMuPiKP[10];
vector<TLorentzVector> pGamma, pEMuPiKP[10];
TClonesArray *ClonesArray[11];
Double_t chisqVtx, chisqKmfit[2];
TClonesArray *GamAf,*GamAf2,*Pi0;
TLorentzVector *Eta, *Eta2, *Pip_b4, *Pim_b4, *Pip_af, *Pim_af, *Pip_af2, *Pim_af2, *cGam_af;

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

JpsiToGamEta4::JpsiToGamEta4(const std::string& name, ISvcLocator* pSvcLocator) : Algorithm(name, pSvcLocator){
   //Declare the properties
   declareProperty("OutputFileName",  m_OutputFileName = "JpsiToGamEta4test.root");
   declareProperty("PsiType",m_psiType = 1);
   declareProperty("Vr0cut", m_vr0cut=1.0);
   declareProperty("Vz0cut", m_vz0cut=10.0);
   declareProperty("Ccoscut", m_ccoscut=0.93);
   declareProperty("IsoAngleCut", m_isoAngleCut=10.0);
   declareProperty("MCTruth", m_MCTruth=false);
}


StatusCode JpsiToGamEta4::initialize(){
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

   for(int i=0;i<11;i++) ClonesArray[i] = new TClonesArray("TLorentzVector");
   GamAf   = new TClonesArray("TLorentzVector");
   GamAf2   = new TClonesArray("TLorentzVector");
   Pi0   = new TClonesArray("TLorentzVector");
   Eta = new TLorentzVector();
   Eta2 = new TLorentzVector();
   Pip_b4 = new TLorentzVector();
   Pim_b4 = new TLorentzVector();
   Pip_af = new TLorentzVector();
   Pim_af = new TLorentzVector();
   Pip_af2 = new TLorentzVector();
   Pim_af2 = new TLorentzVector();
   cGam_af = new TLorentzVector();

   TreeAna = new TTree("TreeAna", "analysis");


   TreeAna->Branch("runid", &runid, "runid/I");
   TreeAna->Branch("evtid", &evtid, "evtid/I");
   TreeAna->Branch("nGoodCharge", &nGoodCharge, "nGoodCharge/I");
   TreeAna->Branch("nGamma", &nGamma, "nGamma/I");
   TreeAna->Branch("chainType", &chainType, "chainType[4]/O");
   TreeAna->Branch("chisqVtx", &chisqVtx, "chisqVtx/D");
   TreeAna->Branch("chisqKmfit", &chisqKmfit, "chisqKmfit[2]/D");
   TreeAna->Branch("gamID", "vector<int>", &gamID);
   TreeAna->Branch("gamID2", "vector<int>", &gamID2);

//   TreeAna->Branch("chargeCut", "vector<double>", &chargeCut);
   TreeAna->Branch("pidCut", "vector<double>", &pidCut);
   TreeAna->Branch("showerCut", "vector<double>", &showerCut);
//   TreeAna->Branch("showerPos", "vector<double>", &showerPos);

   TreeAna->Branch("Gamma", "TClonesArray", &Gamma, 256000, 0);
   TreeAna->Branch("GamAf","TClonesArray",&GamAf,256000,0);
   TreeAna->Branch("GamAf2","TClonesArray",&GamAf2,256000,0);
   TreeAna->Branch("Pi0","TClonesArray",&Pi0,256000,0);
   TreeAna->Branch("Eta", &Eta, 32000,0);
   TreeAna->Branch("Eta2", &Eta2, 32000,0);
   TreeAna->Branch("Pip_b4", &Pip_b4, 32000,0);
   TreeAna->Branch("Pim_b4", &Pim_b4, 32000,0);
   TreeAna->Branch("Pip_af", &Pip_af, 32000,0);
   TreeAna->Branch("Pim_af", &Pim_af, 32000,0);
   TreeAna->Branch("Pip_af2", &Pip_af2, 32000,0);
   TreeAna->Branch("Pim_af2", &Pim_af2, 32000,0);
   //TreeAna->Branch("cGam_af", &cGam_af, 32000,0);

   for(int i=0;i<11;i++) ClonesArray[i]->BypassStreamer();
   GamAf->BypassStreamer();
   GamAf2->BypassStreamer();
   Pi0->BypassStreamer();

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



StatusCode JpsiToGamEta4::execute() {

   //-************Initialize Global Variables**********************
   runid = -1; evtid = -1; chisqVtx = 1000; chisqKmfit[0]=chisqKmfit[1]=1000;
   chainType[0]=chainType[1]=chainType[2]=chainType[3]=false;
   nCharge = 0; nGoodCharge = 0; netGoodCharge = 0; nShower = 0; nGamma = 0;
   for(int i=0;i<10;i++) nEMuPiKP[i]=0;
   iGood.clear(); iGamma.clear();
   for(int i=0;i<10;i++) (iEMuPiKP[i]).clear();
   pGamma.clear();
   for(int i=0;i<10;i++) (pEMuPiKP[i]).clear();

   gamID->clear(); gamID2->clear();
//   chargeCut->clear(); pidCut->clear(); showerCut->clear(); showerPos->clear();
   pidCut->clear();
   showerCut->clear(); 
   for(int i=0;i<11;i++) ClonesArray[i]->Clear();
   GamAf->Clear();
   GamAf2->Clear();
   Pi0->Clear();
   Eta->SetPxPyPzE(0,0,0,0);
   Eta2->SetPxPyPzE(0,0,0,0);
   Pip_b4->SetPxPyPzE(0,0,0,0);
   Pim_b4->SetPxPyPzE(0,0,0,0);
   Pip_af->SetPxPyPzE(0,0,0,0);
   Pim_af->SetPxPyPzE(0,0,0,0);
   Pip_af2->SetPxPyPzE(0,0,0,0);
   Pim_af2->SetPxPyPzE(0,0,0,0);
   cGam_af->SetPxPyPzE(0,0,0,0);
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
   if((nGoodCharge!=0&&nGoodCharge!=2)||netGoodCharge!=0) return StatusCode::SUCCESS;
dataCutFlow2Trk++;
   //-*********** Finish Good Charged Track Selection ***********

if(nGoodCharge==2){
   //-********************PID identification*****************
   int netChargePID = identifyPID(iGood,iEMuPiKP,pidCut);
   for(int i=0;i<10;i++) nEMuPiKP[i] = (iEMuPiKP[i]).size();
   if(npip!=1||npim!=1) return StatusCode::SUCCESS;
dataCutFlowPID++;
   for(int i=0;i<10;i++) assignMomentumToCharged(pEMuPiKP[i], iEMuPiKP[i], nEMuPiKP[i], i/2);
   //-********************End PID identification*****************

   //-********************Good Photon selection*****************
   nGamma = selectGoodPhoton(xorigin,iGamma,showerCut,showerPos);
   if(nGamma<2) return StatusCode::SUCCESS;
dataCutFlowGam++;
   assignMomentumToPhoton(xorigin, pGamma, iGamma, nGamma);
   //-********************End Good Photon selection*****************
}//end if nGoodCharge==2
else if(nGoodCharge==0){
   //-********************Good Photon selection*****************
   nGamma = selectGoodPhoton2(xorigin,iGamma,showerCut,showerPos);
   if(nGamma<3) return StatusCode::SUCCESS;
dataCutFlowGam++;
   assignMomentumToPhoton(xorigin, pGamma, iGamma, nGamma);
   //-********************End Good Photon selection*****************
}

   bool okFit = false; 
   KalmanKinematicFit * kmfit = KalmanKinematicFit::instance();
   RecEmcShower* gTrkForKmfit[7];
   //-********************eta selection*****************

   vector<int> bestID;
   //---------------- 1st chain --------------------------
   if(nGoodCharge==0&&nGamma>=3&&nGamma<20){ //Jpsi->GamEta->Gam GamGam
      chisqKmfit[0]=chisqKmfit[1]=200;
      int bestGam=0;
      for(int igamC=1;igamC<nGamma;igamC++)
         if(pGamma[igamC].E()>pGamma[bestGam].E()) bestGam=igamC;

      for(int igam1=0;igam1<nGamma-1;igam1++){
         if(igam1==bestGam) continue;
         for(int igam2=igam1+1;igam2<nGamma;igam2++){
            if(igam2==bestGam) continue;
               TLorentzVector myEta=pGamma[igam1]+pGamma[igam2];
               if(fabs(myEta.M()-0.55)>0.35) continue;

               gTrkForKmfit[0] = (*(evtRecTrkCol->begin()+iGamma[bestGam]))->emcShower();
               gTrkForKmfit[1] = (*(evtRecTrkCol->begin()+iGamma[igam1]))->emcShower();
               gTrkForKmfit[2] = (*(evtRecTrkCol->begin()+iGamma[igam2]))->emcShower();

               kmfit->init();
               for(int i_add=0;i_add<3;i_add++) kmfit->AddTrack(i_add,0.0,gTrkForKmfit[i_add]);
               kmfit->AddFourMomentum(0, cms);
               if(!kmfit->Fit()) continue;
               if(kmfit->chisq()>=chisqKmfit[0]) continue;
               okFit = true;
               chisqKmfit[0] = kmfit->chisq();
               bestID.clear();
               bestID.push_back(igam1);
               bestID.push_back(igam2);
         }
      }
      if(okFit){
         chainType[0]=true;

         gTrkForKmfit[0] = (*(evtRecTrkCol->begin()+iGamma[bestGam]))->emcShower();
         gTrkForKmfit[1] = (*(evtRecTrkCol->begin()+iGamma[bestID[0]]))->emcShower();
         gTrkForKmfit[2] = (*(evtRecTrkCol->begin()+iGamma[bestID[1]]))->emcShower();
         kmfit->init();
         for(int i_add=0;i_add<3;i_add++) kmfit->AddTrack(i_add,0.0,gTrkForKmfit[i_add]);
         kmfit->AddFourMomentum(0, cms);
         kmfit->Fit();
         chisqKmfit[0] = kmfit->chisq();
         gamID->clear();
         gamID->push_back(bestGam);
         gamID->push_back(bestID[0]);
         gamID->push_back(bestID[1]);
         GamAf->Clear();
         cGam_af->SetPxPyPzE((kmfit->pfit(0)).px(),(kmfit->pfit(0)).py(),(kmfit->pfit(0)).pz(),(kmfit->pfit(0)).e());
         new ((*GamAf)[0]) TLorentzVector(*cGam_af);
         cGam_af->SetPxPyPzE((kmfit->pfit(1)).px(),(kmfit->pfit(1)).py(),(kmfit->pfit(1)).pz(),(kmfit->pfit(1)).e());
         new ((*GamAf)[1]) TLorentzVector(*cGam_af);
         *Eta=*cGam_af;
         cGam_af->SetPxPyPzE((kmfit->pfit(2)).px(),(kmfit->pfit(2)).py(),(kmfit->pfit(2)).pz(),(kmfit->pfit(2)).e());
         new ((*GamAf)[2]) TLorentzVector(*cGam_af);
         *Eta+=*cGam_af;
      }
      //---------------- end 1st chain --------------------------

      //---------------- second chain ---------------------------
      if(nGamma>=7&&nGamma<=17){ //Jpsi->GamEta->Gam 3Pi0
//cout<<event<<": "<<nGamma<<","<<nGoodCharge<<endl;
         okFit=false;

         for(int igam1=0;igam1<nGamma-5;igam1++){
            if(igam1==bestGam) continue;
            for(int igam2=igam1+1;igam2<nGamma-4;igam2++){
               if(igam2==bestGam) continue;
               for(int igam3=igam2+1;igam3<nGamma-3;igam3++){
                  if(igam3==bestGam) continue;
                  for(int igam4=igam3+1;igam4<nGamma-2;igam4++){
                     if(igam4==bestGam) continue;
                     for(int igam5=igam4+1;igam5<nGamma-1;igam5++){
                        if(igam5==bestGam) continue;
                        for(int igam6=igam5+1;igam6<nGamma;igam6++){
                           if(igam6==bestGam) continue;
                           TLorentzVector myEta=pGamma[igam1]+pGamma[igam2]+pGamma[igam3]+pGamma[igam4]+pGamma[igam5]+pGamma[igam6];
                           if(fabs(myEta.M()-0.55)>0.35) continue;

                           gTrkForKmfit[0] = (*(evtRecTrkCol->begin()+iGamma[bestGam]))->emcShower();
                           gTrkForKmfit[1] = (*(evtRecTrkCol->begin()+iGamma[igam1]))->emcShower();
                           gTrkForKmfit[2] = (*(evtRecTrkCol->begin()+iGamma[igam2]))->emcShower();
                           gTrkForKmfit[3] = (*(evtRecTrkCol->begin()+iGamma[igam3]))->emcShower();
                           gTrkForKmfit[4] = (*(evtRecTrkCol->begin()+iGamma[igam4]))->emcShower();
                           gTrkForKmfit[5] = (*(evtRecTrkCol->begin()+iGamma[igam5]))->emcShower();
                           gTrkForKmfit[6] = (*(evtRecTrkCol->begin()+iGamma[igam6]))->emcShower();
                           kmfit->init();
                           for(int i_add=0;i_add<7;i_add++) kmfit->AddTrack(i_add,0.0,gTrkForKmfit[i_add]);
                           //kmfit->AddResonance(0,0.547862,1,2,3,4,5,6);
                           //kmfit->AddFourMomentum(1, cms);
                           kmfit->AddFourMomentum(0, cms);
                           if(!kmfit->Fit()) continue;
                           if(kmfit->chisq()>=chisqKmfit[1]) continue;
                           okFit = true;
                           chisqKmfit[1] = kmfit->chisq();
                           bestID.clear();
                           bestID.push_back(igam1);
                           bestID.push_back(igam2);
                           bestID.push_back(igam3);
                           bestID.push_back(igam4);
                           bestID.push_back(igam5);
                           bestID.push_back(igam6);
                        }
                     }
                  }
               }
            }
         }
         if(okFit){
            okFit=false; chisqKmfit[1]=200;
            vector<int> iPi0;//index of gammas of pi0 (serial number in pGamma)
            vector<TLorentzVector> pPi0;
            if(buildPi0(0.05,0.22,3,iPi0,pPi0,bestID)){
               for(int ipi01=0;ipi01<pPi0.size()-2;ipi01++){
                  for(int ipi02=ipi01+1;ipi02<pPi0.size()-1;ipi02++){
                     if(iPi0[2*ipi02]==iPi0[2*ipi01]||iPi0[2*ipi02]==iPi0[2*ipi01+1]||iPi0[2*ipi02+1]==iPi0[2*ipi01]||iPi0[2*ipi02+1]==iPi0[2*ipi01+1]) continue;
                     for(int ipi03=ipi02+1;ipi03<pPi0.size();ipi03++){
                        if(iPi0[2*ipi02]==iPi0[2*ipi03]||iPi0[2*ipi02]==iPi0[2*ipi03+1]||iPi0[2*ipi02+1]==iPi0[2*ipi03]||iPi0[2*ipi02+1]==iPi0[2*ipi03+1]) continue;
                        if(iPi0[2*ipi01]==iPi0[2*ipi03]||iPi0[2*ipi01]==iPi0[2*ipi03+1]||iPi0[2*ipi01+1]==iPi0[2*ipi03]||iPi0[2*ipi01+1]==iPi0[2*ipi03+1]) continue;

                        gTrkForKmfit[0] = (*(evtRecTrkCol->begin()+iGamma[bestGam]))->emcShower();
                        gTrkForKmfit[1] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*ipi01]]))->emcShower();
                        gTrkForKmfit[2] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*ipi01+1]]))->emcShower();
                        gTrkForKmfit[3] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*ipi02]]))->emcShower();
                        gTrkForKmfit[4] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*ipi02+1]]))->emcShower();
                        gTrkForKmfit[5] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*ipi03]]))->emcShower();
                        gTrkForKmfit[6] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*ipi03+1]]))->emcShower();

                        kmfit->init();
                        for(int i_add=0;i_add<7;i_add++) kmfit->AddTrack(i_add,0.0,gTrkForKmfit[i_add]);
                        kmfit->AddResonance(0,0.1349766,1,2);
                        kmfit->AddResonance(1,0.1349766,3,4);
                        kmfit->AddResonance(2,0.1349766,5,6);
                        kmfit->AddFourMomentum(3, cms);
                        if(!kmfit->Fit()) continue;
                        if(kmfit->chisq()>=chisqKmfit[1]) continue;
                        okFit = true;
                        chisqKmfit[1] = kmfit->chisq();
                        bestID.clear();
                        bestID.push_back(ipi01);
                        bestID.push_back(ipi02);
                        bestID.push_back(ipi03);
                     }
                  }
               }

               if(okFit){
                  chainType[1]=true;

                  gTrkForKmfit[0] = (*(evtRecTrkCol->begin()+iGamma[bestGam]))->emcShower();
                  gTrkForKmfit[1] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*bestID[0]]]))->emcShower();
                  gTrkForKmfit[2] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*bestID[0]+1]]))->emcShower();
                  gTrkForKmfit[3] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*bestID[1]]]))->emcShower();
                  gTrkForKmfit[4] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*bestID[1]+1]]))->emcShower();
                  gTrkForKmfit[5] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*bestID[2]]]))->emcShower();
                  gTrkForKmfit[6] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*bestID[2]+1]]))->emcShower();

                  kmfit->init();
                  for(int i_add=0;i_add<7;i_add++) kmfit->AddTrack(i_add,0.0,gTrkForKmfit[i_add]);
                  kmfit->AddResonance(0,0.1349766,1,2);
                  kmfit->AddResonance(1,0.1349766,3,4);
                  kmfit->AddResonance(2,0.1349766,5,6);
                  kmfit->AddFourMomentum(3, cms);
                  kmfit->Fit();
                  chisqKmfit[1] = kmfit->chisq();
                  gamID2->clear();
                  gamID2->push_back(bestGam);
                  gamID2->push_back(iPi0[2*bestID[0]]);
                  gamID2->push_back(iPi0[2*bestID[0]+1]);
                  gamID2->push_back(iPi0[2*bestID[1]]);
                  gamID2->push_back(iPi0[2*bestID[1]+1]);
                  gamID2->push_back(iPi0[2*bestID[2]]);
                  gamID2->push_back(iPi0[2*bestID[2]+1]);
                  GamAf2->Clear(); Pi0->Clear();
                  TLorentzVector temp;
                  cGam_af->SetPxPyPzE((kmfit->pfit(0)).px(),(kmfit->pfit(0)).py(),(kmfit->pfit(0)).pz(),(kmfit->pfit(0)).e());
                  new ((*GamAf2)[0]) TLorentzVector(*cGam_af);
                  cGam_af->SetPxPyPzE((kmfit->pfit(1)).px(),(kmfit->pfit(1)).py(),(kmfit->pfit(1)).pz(),(kmfit->pfit(1)).e());
                  new ((*GamAf2)[1]) TLorentzVector(*cGam_af);
                  *Eta2=*cGam_af;
                  cGam_af->SetPxPyPzE((kmfit->pfit(2)).px(),(kmfit->pfit(2)).py(),(kmfit->pfit(2)).pz(),(kmfit->pfit(2)).e());
                  new ((*GamAf2)[2]) TLorentzVector(*cGam_af);
                  *Eta2+=*cGam_af; temp=*Eta2;
                  new ((*Pi0)[0]) TLorentzVector(*Eta2);
                  cGam_af->SetPxPyPzE((kmfit->pfit(3)).px(),(kmfit->pfit(3)).py(),(kmfit->pfit(3)).pz(),(kmfit->pfit(3)).e());
                  new ((*GamAf2)[3]) TLorentzVector(*cGam_af);
                  *Eta2=*cGam_af;
                  cGam_af->SetPxPyPzE((kmfit->pfit(4)).px(),(kmfit->pfit(4)).py(),(kmfit->pfit(4)).pz(),(kmfit->pfit(4)).e());
                  new ((*GamAf2)[4]) TLorentzVector(*cGam_af);
                  *Eta2+=*cGam_af; temp+=*Eta2;
                  new ((*Pi0)[1]) TLorentzVector(*Eta2);
                  cGam_af->SetPxPyPzE((kmfit->pfit(5)).px(),(kmfit->pfit(5)).py(),(kmfit->pfit(5)).pz(),(kmfit->pfit(5)).e());
                  new ((*GamAf2)[5]) TLorentzVector(*cGam_af);
                  *Eta2=*cGam_af;
                  cGam_af->SetPxPyPzE((kmfit->pfit(6)).px(),(kmfit->pfit(6)).py(),(kmfit->pfit(6)).pz(),(kmfit->pfit(6)).e());
                  new ((*GamAf2)[6]) TLorentzVector(*cGam_af);
                  *Eta2+=*cGam_af; temp+=*Eta2;
                  new ((*Pi0)[2]) TLorentzVector(*Eta2);
                  *Eta2=temp;
               }//end if okFit 3 pi0
            }//end if buildPi0

         }//end if okfit 6 Gamma

      }
      //------------------- end second chain -------------------------------
/*
      //---------------- second chain - another choice ------------------
      if(nGamma>=7&&nGamma<=17){ //Jpsi->GamEta->Gam 3Pi0
//cout<<event<<": "<<nGamma<<","<<nGoodCharge<<endl;

         vector<int> iPi0;//index of gammas of pi0 (serial number in pGamma)
         vector<TLorentzVector> pPi0;
         for(int igam1=0;igam1<nGamma-1;igam1++){
            if(igam1==bestGam) continue;
            for(int igam2=igam1+1;igam2<nGamma;igam2++){
               if(igam2=bestGam) continue;
               double mass=(pGamma[igam1]+pGamma[igam2]).M();
               if(mass<0.05||mass>0.22) continue;
               gTrkForKmfit[0] = (*(evtRecTrkCol->begin()+iGamma[igam1]))->emcShower();
               gTrkForKmfit[1] = (*(evtRecTrkCol->begin()+iGamma[igam2]))->emcShower();

               kmfit->init();
               for(int i_add=0;i_add<2;i_add++) kmfit->AddTrack(i_add,0.0,gTrkForKmfit[i_add]);
               kmfit->AddResonance(0,0.1349766,0,1);
               if(!kmfit->Fit()) continue;
               if(kmfit->chisq()>25) continue;
               iPi0.push_back(igam1); iPi0.push_back(igam2);
               Eta2->SetPxPyPzE((kmfit->pfit(1)).px()+(kmfit->pfit(0)).px(),(kmfit->pfit(1)).py()+(kmfit->pfit(0)).py(),(kmfit->pfit(1)).pz()+(kmfit->pfit(0)).pz(),(kmfit->pfit(1)).e()+(kmfit->pfit(0)).e());
               pPi0.push_back(*Eta2);
            }
         }

         if(pPi0.size()>=3){
            okFit=false;
            for(int ipi01=0;ipi01<pPi0.size()-2;ipi01++){
               for(int ipi02=ipi01+1;ipi02<pPi0.size()-1;ipi02++){
                  if(iPi0[2*ipi02]==iPi0[2*ipi01]||iPi0[2*ipi02]==iPi0[2*ipi01+1]||iPi0[2*ipi02+1]==iPi0[2*ipi01]||iPi0[2*ipi02+1]==iPi0[2*ipi01+1]) continue;
                  for(int ipi03=ipi02+1;ipi03<pPi0.size();ipi03++){
                     if(iPi0[2*ipi02]==iPi0[2*ipi03]||iPi0[2*ipi02]==iPi0[2*ipi03+1]||iPi0[2*ipi02+1]==iPi0[2*ipi03]||iPi0[2*ipi02+1]==iPi0[2*ipi03+1]) continue;
                     if(iPi0[2*ipi01]==iPi0[2*ipi03]||iPi0[2*ipi01]==iPi0[2*ipi03+1]||iPi0[2*ipi01+1]==iPi0[2*ipi03]||iPi0[2*ipi01+1]==iPi0[2*ipi03+1]) continue;
//cout<<"\t"<<iPi0.size()<<endl;
                     double mass=(pPi0[ipi01]+pPi0[ipi02]+pPi0[ipi03]).M();
                     if(mass<0||mass>1.5) continue;

                     gTrkForKmfit[0] = (*(evtRecTrkCol->begin()+iGamma[bestGam]))->emcShower();
                     gTrkForKmfit[1] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*ipi01]]))->emcShower();
                     gTrkForKmfit[2] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*ipi01+1]]))->emcShower();
                     gTrkForKmfit[3] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*ipi02]]))->emcShower();
                     gTrkForKmfit[4] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*ipi02+1]]))->emcShower();
                     gTrkForKmfit[5] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*ipi03]]))->emcShower();
                     gTrkForKmfit[6] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*ipi03+1]]))->emcShower();

                     kmfit->init();
                     for(int i_add=0;i_add<7;i_add++) kmfit->AddTrack(i_add,0.0,gTrkForKmfit[i_add]);
                     kmfit->AddResonance(0,0.1349766,1,2);
                     kmfit->AddResonance(1,0.1349766,3,4);
                     kmfit->AddResonance(2,0.1349766,5,6);
                     kmfit->AddFourMomentum(3, cms);
                     if(!kmfit->Fit()) continue;
                     if(kmfit->chisq()>=chisqKmfit[1]) continue;
                     okFit = true;
                     chisqKmfit[1] = kmfit->chisq();
                     bestID.clear();
                     bestID.push_back(ipi01);
                     bestID.push_back(ipi02);
                     bestID.push_back(ipi03);
                  }
               }
            }
            if(okFit){
               chainType[1]=true;

               gTrkForKmfit[0] = (*(evtRecTrkCol->begin()+iGamma[bestGam]))->emcShower();
               gTrkForKmfit[1] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*bestID[0]]]))->emcShower();
               gTrkForKmfit[2] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*bestID[0]+1]]))->emcShower();
               gTrkForKmfit[3] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*bestID[1]]]))->emcShower();
               gTrkForKmfit[4] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*bestID[1]+1]]))->emcShower();
               gTrkForKmfit[5] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*bestID[2]]]))->emcShower();
               gTrkForKmfit[6] = (*(evtRecTrkCol->begin()+iGamma[iPi0[2*bestID[2]+1]]))->emcShower();

               kmfit->init();
               for(int i_add=0;i_add<7;i_add++) kmfit->AddTrack(i_add,0.0,gTrkForKmfit[i_add]);
               kmfit->AddResonance(0,0.1349766,1,2);
               kmfit->AddResonance(1,0.1349766,3,4);
               kmfit->AddResonance(2,0.1349766,5,6);
               kmfit->AddFourMomentum(3, cms);
               kmfit->Fit();
               chisqKmfit[1] = kmfit->chisq();
               gamID2->clear();
               gamID2->push_back(bestGam);
               gamID2->push_back(iPi0[2*bestID[0]]);
               gamID2->push_back(iPi0[2*bestID[0]+1]);
               gamID2->push_back(iPi0[2*bestID[1]]);
               gamID2->push_back(iPi0[2*bestID[1]+1]);
               gamID2->push_back(iPi0[2*bestID[2]]);
               gamID2->push_back(iPi0[2*bestID[2]+1]);
               GamAf2->Clear(); Pi0->Clear();
               TLorentzVector temp;
               cGam_af->SetPxPyPzE((kmfit->pfit(0)).px(),(kmfit->pfit(0)).py(),(kmfit->pfit(0)).pz(),(kmfit->pfit(0)).e());
               new ((*GamAf2)[0]) TLorentzVector(*cGam_af);
               cGam_af->SetPxPyPzE((kmfit->pfit(1)).px(),(kmfit->pfit(1)).py(),(kmfit->pfit(1)).pz(),(kmfit->pfit(1)).e());
               new ((*GamAf2)[1]) TLorentzVector(*cGam_af);
               *Eta2=*cGam_af;
               cGam_af->SetPxPyPzE((kmfit->pfit(2)).px(),(kmfit->pfit(2)).py(),(kmfit->pfit(2)).pz(),(kmfit->pfit(2)).e());
               new ((*GamAf2)[2]) TLorentzVector(*cGam_af);
               *Eta2+=*cGam_af; temp=*Eta2;
               new ((*Pi0)[0]) TLorentzVector(*Eta2);
               cGam_af->SetPxPyPzE((kmfit->pfit(3)).px(),(kmfit->pfit(3)).py(),(kmfit->pfit(3)).pz(),(kmfit->pfit(3)).e());
               new ((*GamAf2)[3]) TLorentzVector(*cGam_af);
               *Eta2=*cGam_af;
               cGam_af->SetPxPyPzE((kmfit->pfit(4)).px(),(kmfit->pfit(4)).py(),(kmfit->pfit(4)).pz(),(kmfit->pfit(4)).e());
               new ((*GamAf2)[4]) TLorentzVector(*cGam_af);
               *Eta2+=*cGam_af; temp+=*Eta2;
               new ((*Pi0)[1]) TLorentzVector(*Eta2);
               cGam_af->SetPxPyPzE((kmfit->pfit(5)).px(),(kmfit->pfit(5)).py(),(kmfit->pfit(5)).pz(),(kmfit->pfit(5)).e());
               new ((*GamAf2)[5]) TLorentzVector(*cGam_af);
               *Eta2=*cGam_af;
               cGam_af->SetPxPyPzE((kmfit->pfit(6)).px(),(kmfit->pfit(6)).py(),(kmfit->pfit(6)).pz(),(kmfit->pfit(6)).e());
               new ((*GamAf2)[6]) TLorentzVector(*cGam_af);
               *Eta2+=*cGam_af; temp+=*Eta2;
               new ((*Pi0)[2]) TLorentzVector(*Eta2);
               *Eta2=temp;
            }
         }// end if pPi0.size()>3
      }
      //------------------- end second chain - another choice -----------------
*/
   }

   else if(nGamma>=2&&nGoodCharge==2&&nGamma<20){ //Jpsi->GamEta->Gam Pi0PipPim
      VertexFit* vtxfit = VertexFit::instance();
      vector<int> trkID, trkType;
      trkID.push_back(ipip[0]);trkType.push_back(2);
      trkID.push_back(ipim[0]);trkType.push_back(2);
      if(!doVertexFit(vtxfit,2,trkID,trkType)) return StatusCode::SUCCESS;
      chisqVtx=vtxfit->chisq();
      WTrackParameter wpip = vtxfit->wtrk(0);
      WTrackParameter wpim = vtxfit->wtrk(1);

      int bestGam=0;
      for(int igamC=1;igamC<nGamma;igamC++)
         if(pGamma[igamC].E()>pGamma[bestGam].E()) bestGam=igamC;
      //------------------------- third chain ---------------------------
      chisqKmfit[0]=chisqKmfit[1]=200; okFit=false;
      for(int igam1=0;igam1<nGamma;igam1++){
         if(igam1==bestGam) continue;
         double mass=(pGamma[igam1]+ppip[0]+ppim[0]).M();
         if(fabs(mass-0.55)>0.35) continue;
         gTrkForKmfit[0] = (*(evtRecTrkCol->begin()+iGamma[bestGam]))->emcShower();
         gTrkForKmfit[1] = (*(evtRecTrkCol->begin()+iGamma[igam1]))->emcShower();

         kmfit->init();
         for(int i_add=0;i_add<2;i_add++) kmfit->AddTrack(i_add,0.0,gTrkForKmfit[i_add]);
         kmfit->AddTrack(2,wpip);
         kmfit->AddTrack(3,wpim);
         kmfit->AddFourMomentum(0, cms);
         if(!kmfit->Fit()) continue;
         if(kmfit->chisq()>=chisqKmfit[0]) continue;
         okFit = true;
         chisqKmfit[0] = kmfit->chisq();
         bestID.clear(); bestID.push_back(igam1);
      }
      if(okFit){
         chainType[2]=true;

         gTrkForKmfit[0] = (*(evtRecTrkCol->begin()+iGamma[bestGam]))->emcShower();
         gTrkForKmfit[1] = (*(evtRecTrkCol->begin()+iGamma[bestID[0]]))->emcShower();
         kmfit->init();
         for(int i_add=0;i_add<2;i_add++) kmfit->AddTrack(i_add,0.0,gTrkForKmfit[i_add]);
         kmfit->AddTrack(2,wpip);
         kmfit->AddTrack(3,wpim);
         kmfit->AddFourMomentum(0, cms);
         kmfit->Fit();
         chisqKmfit[0] = kmfit->chisq();
         gamID->clear();
         gamID->push_back(bestGam);
         gamID->push_back(bestID[0]);
         GamAf->Clear();
         cGam_af->SetPxPyPzE((kmfit->pfit(0)).px(),(kmfit->pfit(0)).py(),(kmfit->pfit(0)).pz(),(kmfit->pfit(0)).e());
         new ((*GamAf)[0]) TLorentzVector(*cGam_af);
         cGam_af->SetPxPyPzE((kmfit->pfit(1)).px(),(kmfit->pfit(1)).py(),(kmfit->pfit(1)).pz(),(kmfit->pfit(1)).e());
         new ((*GamAf)[1]) TLorentzVector(*cGam_af);
         *Eta=*cGam_af;
         Pip_af->SetPxPyPzE((kmfit->pfit(2)).px(),(kmfit->pfit(2)).py(),(kmfit->pfit(2)).pz(),(kmfit->pfit(2)).e());
         Pim_af->SetPxPyPzE((kmfit->pfit(3)).px(),(kmfit->pfit(3)).py(),(kmfit->pfit(3)).pz(),(kmfit->pfit(3)).e());
         *Eta+=*Pip_af+*Pim_af;
      }
      //------------------------- end third chain ---------------------------

      //------------------------- fourth chain ---------------------------
      if(nGamma>=3){ //Jpsi->GamEta->Gam GamPipPim
         okFit=false;
         for(int igam1=0;igam1<nGamma-1;igam1++){
            if(igam1==bestGam) continue;
            for(int igam2=igam1+1;igam2<nGamma;igam2++){
               if(igam2==bestGam) continue;

               double mass=(pGamma[igam1]+pGamma[igam2]+ppip[0]+ppim[0]).M();
               //if(mass<0.06||mass>0.21) continue;
               if(fabs(mass-0.55)>0.35) continue;

               gTrkForKmfit[0] = (*(evtRecTrkCol->begin()+iGamma[bestGam]))->emcShower();
               gTrkForKmfit[1] = (*(evtRecTrkCol->begin()+iGamma[igam1]))->emcShower();
               gTrkForKmfit[2] = (*(evtRecTrkCol->begin()+iGamma[igam2]))->emcShower();

               kmfit->init();
               for(int i_add=0;i_add<3;i_add++) kmfit->AddTrack(i_add,0.0,gTrkForKmfit[i_add]);
               kmfit->AddTrack(3,wpip);
               kmfit->AddTrack(4,wpim);
               kmfit->AddResonance(0,0.1349766,1,2);
               kmfit->AddFourMomentum(1, cms);
               if(!kmfit->Fit()) continue;
               if(kmfit->chisq()>=chisqKmfit[1]) continue;
               okFit = true;
               chisqKmfit[1] = kmfit->chisq();
               bestID.clear();
               bestID.push_back(igam1);
               bestID.push_back(igam2);
            }
         }
         if(okFit){
            chainType[3]=true;

            gTrkForKmfit[0] = (*(evtRecTrkCol->begin()+iGamma[bestGam]))->emcShower();
            gTrkForKmfit[1] = (*(evtRecTrkCol->begin()+iGamma[bestID[0]]))->emcShower();
            gTrkForKmfit[2] = (*(evtRecTrkCol->begin()+iGamma[bestID[1]]))->emcShower();

            kmfit->init();
            for(int i_add=0;i_add<3;i_add++) kmfit->AddTrack(i_add,0.0,gTrkForKmfit[i_add]);
            kmfit->AddTrack(3,wpip);
            kmfit->AddTrack(4,wpim);
            kmfit->AddResonance(0,0.1349766,1,2);
            kmfit->AddFourMomentum(1, cms);
            kmfit->Fit();
            chisqKmfit[1] = kmfit->chisq();
            gamID2->clear();
            gamID2->push_back(bestGam);
            gamID2->push_back(bestID[0]);
            gamID2->push_back(bestID[1]);
            GamAf2->Clear();
            cGam_af->SetPxPyPzE((kmfit->pfit(0)).px(),(kmfit->pfit(0)).py(),(kmfit->pfit(0)).pz(),(kmfit->pfit(0)).e());
            new ((*GamAf2)[0]) TLorentzVector(*cGam_af);
            cGam_af->SetPxPyPzE((kmfit->pfit(1)).px(),(kmfit->pfit(1)).py(),(kmfit->pfit(1)).pz(),(kmfit->pfit(1)).e());
            new ((*GamAf2)[1]) TLorentzVector(*cGam_af);
            *Eta2=*cGam_af;
            cGam_af->SetPxPyPzE((kmfit->pfit(2)).px(),(kmfit->pfit(2)).py(),(kmfit->pfit(2)).pz(),(kmfit->pfit(2)).e());
            new ((*GamAf2)[2]) TLorentzVector(*cGam_af);
            *Eta2+=*cGam_af;
            Pip_af2->SetPxPyPzE((kmfit->pfit(3)).px(),(kmfit->pfit(3)).py(),(kmfit->pfit(3)).pz(),(kmfit->pfit(3)).e());
            Pim_af2->SetPxPyPzE((kmfit->pfit(4)).px(),(kmfit->pfit(4)).py(),(kmfit->pfit(4)).pz(),(kmfit->pfit(4)).e());
            Pi0->Clear(); new ((*Pi0)[0]) TLorentzVector(*Eta2);
            *Eta2+=*Pip_af2+*Pim_af2;
         }
      }
      //------------------------- end fourth chain ---------------------------

      *Pip_b4=ppip[0]; *Pim_b4=ppim[0];
   }
   if(!(chainType[0]||chainType[1]||chainType[2]||chainType[3])) return StatusCode::SUCCESS;
   //-********************end eta selection*****************

   for(int i=0;i<pGamma.size();i++) new ((*Gamma)[i]) TLorentzVector(pGamma[i]);

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

   runid = runNo; evtid = event;
//if((nShower*6)!=showerCut->size()) cout<<runNo<<","<<event<<": nShower*6!=showerCut.size()"<<endl;
   TreeAna->Fill();
   return StatusCode::SUCCESS;
}


StatusCode JpsiToGamEta4::finalize() {
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
bool JpsiToGamEta4::getPrimaryVertex(Hep3Vector &xorigin, HepSymMatrix &VtxErr){
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
int JpsiToGamEta4::selectGoodChargedTrack(Hep3Vector xorigin, vector<int> &iGood, vector<double> *chargeCut){
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
int JpsiToGamEta4::selectGoodPhoton(Hep3Vector xorigin, vector<int> &iGamma, vector<double> *showerCut, vector<double> *showerPos){
   for(int i = (*_evtRecEvent)->totalCharged(); i< (*_evtRecEvent)->totalTracks(); i++){
      if(i >= (*_evtRecTrkCol)->size()) break;

      EvtRecTrackIterator itTrk=(*_evtRecTrkCol)->begin() + i;
      if(!(*itTrk)->isEmcShowerValid()) continue;
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      if(emcTrk->time()>14 || emcTrk->time()<0) continue;
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

      showerCut->push_back(cos(the));
//      showerCut->push_back(eraw); showerCut->push_back(dAngKlong); showerCut->push_back(dang); showerCut->push_back(the); showerCut->push_back(emcpos.phi()); showerCut->push_back(emcTrk->time());
//      showerPos->push_back(emcTrk->x()); showerPos->push_back(emcTrk->y()); showerPos->push_back(emcTrk->z());
      if(eraw < e_threshold) continue;
      if(fabs(dang) < m_isoAngleCut) continue;
      iGamma.push_back((*itTrk)->trackId());
   }
   return iGamma.size();
}
//-----------End Good Photon Selection-------------

//-----------Good Photon Selection 2-------------
int JpsiToGamEta4::selectGoodPhoton2(Hep3Vector xorigin, vector<int> &iGamma, vector<double> *showerCut, vector<double> *showerPos){
   double bestShowerE=0, initialTime=-100;
   for(int i = (*_evtRecEvent)->totalCharged(); i< (*_evtRecEvent)->totalTracks(); i++){
      if(i >= (*_evtRecTrkCol)->size()) break;
      EvtRecTrackIterator itTrk=(*_evtRecTrkCol)->begin() + i;
      if(!(*itTrk)->isEmcShowerValid()) continue;
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      if(emcTrk->energy()>bestShowerE){ bestShowerE=emcTrk->energy(); initialTime=emcTrk->time(); }
   }
   if(bestShowerE<0.025||initialTime==-100) return 0;

   for(int i = (*_evtRecEvent)->totalCharged(); i< (*_evtRecEvent)->totalTracks(); i++){
      if(i >= (*_evtRecTrkCol)->size()) break;

      EvtRecTrackIterator itTrk=(*_evtRecTrkCol)->begin() + i;
      if(!(*itTrk)->isEmcShowerValid()) continue;
      RecEmcShower *emcTrk = (*itTrk)->emcShower();
      if(fabs(emcTrk->time()-initialTime)>10) continue;
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

      showerCut->push_back(cos(the));
//      showerCut->push_back(eraw); showerCut->push_back(dAngKlong); showerCut->push_back(dang); showerCut->push_back(the); showerCut->push_back(emcpos.phi()); showerCut->push_back(emcTrk->time());
//      showerPos->push_back(emcTrk->x()); showerPos->push_back(emcTrk->y()); showerPos->push_back(emcTrk->z());
      if(eraw < e_threshold) continue;
      if(fabs(dang) < m_isoAngleCut) continue;
      //if(emcTrk->time()>14  || emcTrk->time()<0) continue;
      iGamma.push_back((*itTrk)->trackId());
   }
   return iGamma.size();
}
//-----------End Good Photon Selection 2-------------

//-----------PID identification-------------
int JpsiToGamEta4::identifyPID(vector<int> iGood, vector<int> iEMuPiKP[], vector<double> *pidCut, int PID_0, int PID_1, int PID_2, int PID_3, int PID_4){
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

      pidCut->push_back(maxPID); pidCut->push_back(prob_pid[maxPID]);
      maxPID=2;
      if(true){
      //if((maxPID==PID_0||maxPID==PID_1||maxPID==PID_2||maxPID==PID_3||maxPID==PID_4||PID_0+PID_1+PID_2+PID_3+PID_4==-5) && prob_pid[maxPID]>0.001){
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
bool JpsiToGamEta4::doVertexFit(VertexFit* vtxfit, int nTrk, vector<int> trkID, vector<int> trkType){

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
bool JpsiToGamEta4::doSecondVertexFit(SecondVertexFit *svtxfit, VertexFit* vtxfit, Hep3Vector xorigin, HepSymMatrix VtxErr){
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
void JpsiToGamEta4::assignMomentumToPhoton(Hep3Vector xorigin, Vp4 &pGamma, Vint iGamma, Int_t nGamma){
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
void JpsiToGamEta4::assignMomentumToCharged(Vp4 &pCharged, Vint iCharged, Int_t nGoodCharge, int pidType){
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
void JpsiToGamEta4::assignMomentumToCharged(Vp4 &pCharged, Vint iCharged, Int_t nGoodCharge, int pidType, vector<double> *emcE){
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
bool JpsiToGamEta4::buildPi0(double massWin_lower, double massWin_upper, vector<int> &gammaID, vector<TLorentzVector> &pPi0){
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
bool JpsiToGamEta4::buildPi0(double massWin_lower, double massWin_upper, vector<int> &gammaID, vector<TLorentzVector> &pPi0, vector<int> inputGamID){
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
bool JpsiToGamEta4::buildPi0(double massWin_lower, double massWin_upper, int nPi0ToBuild, vector<int> &gammaID, vector<TLorentzVector> &pPi0, vector<int> inputGamID){
//gammaID:gamma在pGamma中的排序，从0开始
   if(inputGamID.size() < 2*nPi0ToBuild){
      cout<<"ERROR! There're less than "<<2*nPi0ToBuild<<" gammas in inputGamID to build "<<nPi0ToBuild<<" #pi0"<<endl;
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
bool JpsiToGamEta4::buildPi0(double massWin_lower, double massWin_upper, int nPi0ToBuild, vector<int> &gammaID, vector<TLorentzVector> &pPi0){
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

// Kinematic fit correction: thx Guoyp and Liaogy
// corset(): sets up the generation by calculating C from V.
void JpsiToGamEta4::corset(HepSymMatrix &V, HepMatrix &C, int n){
        double sum;
        // Compute square root of matrix sigma
        for (int j=0; j<n; j++) {
                sum = 0;
                for (int k=0; k<j; k++) {
                        sum = sum + C[j][k]*C[j][k];
                }
                C[j][j] = sqrt(abs(V[j][j] - sum));
                // Off Diagonal terms
                for (int i=j+1; i<n; i++) {
                        sum = 0;
                        for (int k=0; k<j; k++) {
                                sum = sum + C[i][k]*C[j][k];
                        }
                        C[i][j] = (V[i][j] - sum)/C[j][j];
                }
        }
}

//
void JpsiToGamEta4::corgen(HepMatrix &C, HepVector &x, int n)
{
        int i,j;
        int nmax = 100;

        if (n > nmax ) {
                printf("Error in corgen: array overflown");
        }

        double tmp[3];
        for(int p = 0 ; p < n; p ++){
                tmp[p] = gRandom->Gaus(0,1);
                //cout<<"tmp["<<p<<"]="<<tmp[p]<<endl;
        }
        for ( i=0; i<n; i++) {
                x[i] = 0.0;
                for (j=0; j<=i; j++) {
                        x[i] = x[i]+C[i][j]*tmp[j];
                }
        }

}

void JpsiToGamEta4::calibration(RecMdcKalTrack *trk , HepVector &wtrk_zHel, int n){
 
   HepVector pip_calerr_d2(5,0);
   HepVector pim_calerr_d2(5,0);
   HepVector kp_calerr_d2(5,0);
   HepVector km_calerr_d2(5,0);
   HepVector prtp_calerr_d2(5,0);
   HepVector prtm_calerr_d2(5,0);
   //wangzh 
   pip_calerr_d2[0] = 1.0;
   pip_calerr_d2[1] = 1.207;
   pip_calerr_d2[2] = 1.211;
   pip_calerr_d2[3] = 1.0;
   pip_calerr_d2[4] = 1.138;

   pim_calerr_d2[0] = 1.0;
   pim_calerr_d2[1] = 1.149;
   pim_calerr_d2[2] = 1.263;
   pim_calerr_d2[3] = 1.0;
   pim_calerr_d2[4] = 1.117;
 
   kp_calerr_d2[0] = 1.0;
   kp_calerr_d2[1] = 1.195;
   kp_calerr_d2[2] = 1.214;//1.28;
   kp_calerr_d2[3] = 1.0;
   kp_calerr_d2[4] = 1.156;

   km_calerr_d2[0] = 1.0;
   km_calerr_d2[1] = 1.221;//1.204;
   km_calerr_d2[2] = 1.238;//1.236;
   km_calerr_d2[3] = 1.0;
   km_calerr_d2[4] = 1.256;

   prtp_calerr_d2[0] = 1.0;
   prtp_calerr_d2[1] = 1.125;
   prtp_calerr_d2[2] = 1.168;
   prtp_calerr_d2[3] = 1.0;
   prtp_calerr_d2[4] = 1.090;
 
   prtm_calerr_d2[0] = 1.0;
   prtm_calerr_d2[1] = 1.107;
   prtm_calerr_d2[2] = 1.112;
   prtm_calerr_d2[3] = 1.0;
   prtm_calerr_d2[4] = 1.062;

   HepVector pip_calmean_d2(5,0);
   HepVector pim_calmean_d2(5,0);
   HepVector kp_calmean_d2(5,0);
   HepVector km_calmean_d2(5,0);
   HepVector prtp_calmean_d2(5,0);
   HepVector prtm_calmean_d2(5,0);
   //mean 
   pip_calmean_d2[0] = 0;
   pip_calmean_d2[1] = -0.253;
   pip_calmean_d2[2] = 0.092;
   pip_calmean_d2[3] = 0;
   pip_calmean_d2[4] = 0.232;//0.454;

   pim_calmean_d2[0] = 0;
   pim_calmean_d2[1] =-0.083;
   pim_calmean_d2[2] = -0.0911;
   pim_calmean_d2[3] = 0;
   pim_calmean_d2[4] = 0.2674;//0.484;
 
   kp_calmean_d2[0] = 0;
   kp_calmean_d2[1] = -0.031;
   kp_calmean_d2[2] = 0.246;
   kp_calmean_d2[3] = 0;
   kp_calmean_d2[4] = 0.177;//0.34;

   km_calmean_d2[0] = 0;
   km_calmean_d2[1] = -0.058;
   km_calmean_d2[2] = -0.093;//-0.15;
   km_calmean_d2[3] = 0;
   km_calmean_d2[4] = 0.0284;//0.14;
   prtp_calmean_d2[0] = 0;
   prtp_calmean_d2[1] = 0;
   prtp_calmean_d2[2] = 0;
   prtp_calmean_d2[3] = 0;
   prtp_calmean_d2[4] = 0;
   prtm_calmean_d2[0] = 0;
   prtm_calmean_d2[1] = 0;
   prtm_calmean_d2[2] = 0;
   prtm_calmean_d2[3] = 0;
   prtm_calmean_d2[4] = 0;
 
   if(trk->charge()>0 && n==0){
      //pip
      HepSymMatrix wpip_zerr(5,0);
      wpip_zerr = trk->getZError();
      HepSymMatrix wpip_zcal(3,0);
 
      wpip_zcal[0][0] = (pip_calerr_d2[1]*pip_calerr_d2[1]-1)*wpip_zerr[1][1];
      wpip_zcal[1][1] = (pip_calerr_d2[2]*pip_calerr_d2[2]-1)*wpip_zerr[2][2];
      wpip_zcal[2][2] = (pip_calerr_d2[4]*pip_calerr_d2[4]-1)*wpip_zerr[4][4];
 
      HepMatrix wpip_zerrc(3,3,0);
      helixPhiEtap::corset(wpip_zcal,wpip_zerrc,3);
      HepVector wpip_zgen(3,0);
      helixPhiEtap::corgen(wpip_zerrc,wpip_zgen,3);
 
      wtrk_zHel[0] = trk->getZHelix()[0];
      wtrk_zHel[1] = trk->getZHelix()[1]+pip_calmean_d2[1]*sqrt(wpip_zerr[1][1])+wpip_zgen[0];
      wtrk_zHel[2] = trk->getZHelix()[2]+pip_calmean_d2[2]*sqrt(wpip_zerr[2][2])+wpip_zgen[1];
      wtrk_zHel[3] = trk->getZHelix()[3];
      wtrk_zHel[4] = trk->getZHelix()[4]+pip_calmean_d2[4]*sqrt(wpip_zerr[4][4])+wpip_zgen[2];
 
   }
 
   if(trk->charge()<0 && n==0)
   {
      //pim
      HepSymMatrix wpim_zerr(5,0);
      wpim_zerr = trk->getZError();
 
      HepSymMatrix wpim_zcal(3,0);
 
      wpim_zcal[0][0] = (pim_calerr_d2[1]*pim_calerr_d2[1]-1)*wpim_zerr[1][1];
      wpim_zcal[1][1] = (pim_calerr_d2[2]*pim_calerr_d2[2]-1)*wpim_zerr[2][2];
      wpim_zcal[2][2] = (pim_calerr_d2[4]*pim_calerr_d2[4]-1)*wpim_zerr[4][4];
 
      HepMatrix wpim_zerrc(3,3,0);
      helixPhiEtap::corset(wpim_zcal,wpim_zerrc,3);
      HepVector wpim_zgen(3,0);
      helixPhiEtap::corgen(wpim_zerrc,wpim_zgen,3);
 
 
      wtrk_zHel[0] = trk->getZHelix()[0];
      wtrk_zHel[1] = trk->getZHelix()[1]+pim_calmean_d2[1]*sqrt(wpim_zerr[1][1])+wpim_zgen[0];
      wtrk_zHel[2] = trk->getZHelix()[2]+pim_calmean_d2[2]*sqrt(wpim_zerr[2][2])+wpim_zgen[1];
      wtrk_zHel[3] = trk->getZHelix()[3];
      wtrk_zHel[4] = trk->getZHelix()[4]+pim_calmean_d2[4]*sqrt(wpim_zerr[4][4])+wpim_zgen[2];
 
   }
 
   if(trk->charge()>0 && n==1)
   {
      //kp
      HepSymMatrix wkp_zerr(5,0);
      wkp_zerr = trk->getZErrorK();
 
      HepSymMatrix wkp_zcal(3,0);
 
      wkp_zcal[0][0] = (kp_calerr_d2[1]*kp_calerr_d2[1]-1)*wkp_zerr[1][1];
      wkp_zcal[1][1] = (kp_calerr_d2[2]*kp_calerr_d2[2]-1)*wkp_zerr[2][2];
      wkp_zcal[2][2] = (kp_calerr_d2[4]*kp_calerr_d2[4]-1)*wkp_zerr[4][4];
 
      HepMatrix wkp_zerrc(3,3,0);
      helixPhiEtap::corset(wkp_zcal,wkp_zerrc,3);
      HepVector wkp_zgen(3,0);
      helixPhiEtap::corgen(wkp_zerrc,wkp_zgen,3);
 
      wtrk_zHel[0] = trk->getZHelixK()[0];
      wtrk_zHel[1] = trk->getZHelixK()[1]+kp_calmean_d2[1]*sqrt(wkp_zerr[1][1])+wkp_zgen[0];
      wtrk_zHel[2] = trk->getZHelixK()[2]+kp_calmean_d2[2]*sqrt(wkp_zerr[2][2])+wkp_zgen[1];
      wtrk_zHel[3] = trk->getZHelixK()[3];
      wtrk_zHel[4] = trk->getZHelixK()[4]+kp_calmean_d2[4]*sqrt(wkp_zerr[4][4])+wkp_zgen[2];
 
   }
 
   if(trk->charge()<0 && n==1)
   {
      //km
      HepSymMatrix wkm_zerr(5,0);
      wkm_zerr = trk->getZErrorK();
 
      HepSymMatrix wkm_zcal(3,0);
 
      wkm_zcal[0][0] = (km_calerr_d2[1]*km_calerr_d2[1]-1)*wkm_zerr[1][1];
      wkm_zcal[1][1] = (km_calerr_d2[2]*km_calerr_d2[2]-1)*wkm_zerr[2][2];
      wkm_zcal[2][2] = (km_calerr_d2[4]*km_calerr_d2[4]-1)*wkm_zerr[4][4];
 
      HepMatrix wkm_zerrc(3,3,0);
      helixPhiEtap::corset(wkm_zcal,wkm_zerrc,3);
      HepVector wkm_zgen(3,0);
      helixPhiEtap::corgen(wkm_zerrc,wkm_zgen,3);
 
      wtrk_zHel[0] = trk->getZHelixK()[0];
      wtrk_zHel[1] = trk->getZHelixK()[1]+km_calmean_d2[1]*sqrt(wkm_zerr[1][1])+wkm_zgen[0];
      wtrk_zHel[2] = trk->getZHelixK()[2]+km_calmean_d2[2]*sqrt(wkm_zerr[2][2])+wkm_zgen[1];
      wtrk_zHel[3] = trk->getZHelixK()[3];
      wtrk_zHel[4] = trk->getZHelixK()[4]+km_calmean_d2[4]*sqrt(wkm_zerr[4][4])+wkm_zgen[2];
 
   }
 
   if(trk->charge()>0 && n==2)
   {
      //prtp
      HepSymMatrix wprtp_zerr(5,0);
      wprtp_zerr = trk->getZErrorP();
 
      HepSymMatrix wprtp_zcal(3,0);
 
      wprtp_zcal[0][0] = (prtp_calerr_d2[1]*prtp_calerr_d2[1]-1)*wprtp_zerr[1][1];
      wprtp_zcal[1][1] = (prtp_calerr_d2[2]*prtp_calerr_d2[2]-1)*wprtp_zerr[2][2];
      wprtp_zcal[2][2] = (prtp_calerr_d2[4]*prtp_calerr_d2[4]-1)*wprtp_zerr[4][4];
 
      HepMatrix wprtp_zerrc(3,3,0);
      helixPhiEtap::corset(wprtp_zcal,wprtp_zerrc,3);
      HepVector wprtp_zgen(3,0);
      helixPhiEtap::corgen(wprtp_zerrc,wprtp_zgen,3);
 
      wtrk_zHel[0] = trk->getZHelixP()[0];
      wtrk_zHel[1] = trk->getZHelixP()[1]+prtp_calmean_d2[1]*sqrt(wprtp_zerr[1][1])+wprtp_zgen[0];
      wtrk_zHel[2] = trk->getZHelixP()[2]+prtp_calmean_d2[2]*sqrt(wprtp_zerr[2][2])+wprtp_zgen[1];
      wtrk_zHel[3] = trk->getZHelixP()[3];
      wtrk_zHel[4] = trk->getZHelixP()[4]+prtp_calmean_d2[4]*sqrt(wprtp_zerr[4][4])+wprtp_zgen[2];
 
   }
 
   if(trk->charge()<0 && n==2)
   {
      //prtm
      HepSymMatrix wprtm_zerr(5,0);
      wprtm_zerr = trk->getZErrorP();
 
      HepSymMatrix wprtm_zcal(3,0);
 
      wprtm_zcal[0][0] = (prtm_calerr_d2[1]*prtm_calerr_d2[1]-1)*wprtm_zerr[1][1];
      wprtm_zcal[1][1] = (prtm_calerr_d2[2]*prtm_calerr_d2[2]-1)*wprtm_zerr[2][2];
      wprtm_zcal[2][2] = (prtm_calerr_d2[4]*prtm_calerr_d2[4]-1)*wprtm_zerr[4][4];
 
      HepMatrix wprtm_zerrc(3,3,0);
      helixPhiEtap::corset(wprtm_zcal,wprtm_zerrc,3);
      HepVector wprtm_zgen(3,0);
      helixPhiEtap::corgen(wprtm_zerrc,wprtm_zgen,3);
 
      wtrk_zHel[0] = trk->getZHelixP()[0];
      wtrk_zHel[1] = trk->getZHelixP()[1]+prtm_calmean_d2[1]*sqrt(wprtm_zerr[1][1])+wprtm_zgen[0];
      wtrk_zHel[2] = trk->getZHelixP()[2]+prtm_calmean_d2[2]*sqrt(wprtm_zerr[2][2])+wprtm_zgen[1];
      wtrk_zHel[3] = trk->getZHelixP()[3];
      wtrk_zHel[4] = trk->getZHelixP()[4]+prtm_calmean_d2[4]*sqrt(wprtm_zerr[4][4])+wprtm_zgen[2];
 
   }
 
 }
