// Include files

// from Gaudi
#include "GaudiKernel/ToolFactory.h"

// local
#include "TupleToolApplyIsolation.h"

#include <Kernel/GetIDVAlgorithm.h>
#include <Kernel/IDVAlgorithm.h>
#include <Kernel/IDistanceCalculator.h>
#include <Kernel/IVertexFit.h>

#include "GaudiAlg/Tuple.h"
#include "GaudiAlg/TupleObj.h"

#include "Event/Particle.h"
#include "Event/MCParticle.h"
#include "GaudiKernel/DeclareFactoryEntries.h"
#include <functional>
#include "TrackInterfaces/ITrackVertexer.h"
#include "Linker/LinkerTable.h"
#include "Kernel/IPVReFitter.h"

#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TPluginManager.h"
#include "TStopwatch.h"


//#include "TMVA/Reader.h"
//#include "TMVA/Config.h"
//#include "TMVA/Tools.h"


//-----------------------------------------------------------------------------
// Implementation file for class : TupleToolVtxIsoln
//
// @author Mitesh Patel, Patrick Koppenburg
// @date   2008-04-15
//-----------------------------------------------------------------------------

// Declaration of the Tool Factory
// actually acts as a using namespace TupleTool
DECLARE_TOOL_FACTORY( TupleToolApplyIsolation );

using namespace LHCb;
//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolApplyIsolation::TupleToolApplyIsolation( const std::string& type,
                                      const std::string& name,
                                      const IInterface* parent )
  : TupleToolBase ( type, name , parent )
    , m_dva(0) 
    , m_dist(0)
    , m_pVertexFit(0)
{
  declareInterface<IParticleTupleTool>(this);

  m_inputParticles.push_back("/Event/Phys/StdAllNoPIDsPions");
  m_inputParticles.push_back("/Event/Phys/StdNoPIDsUpPions");
  m_inputParticles.push_back("Phys/StdNoPIDsVeloPions");
  //m_inputParticles.push_back("/Event/Phys/StdNoPIDsVeloElectrons");

  
  //havent removed / added any of this yet
  declareProperty( "MaxDeltaChi2", m_deltaChi2 = 9.0);
  declareProperty( "MaxChi2", m_Chi2 = 9.0);
  declareProperty( "VertexFit", m_typeVertexFit = "default");
  declareProperty("InputParticles", m_inputParticles );;
  declareProperty("OutputSuffix", m_outputSuffix = "" );
  declareProperty("WeightsFile", m_weightsName = "weights.xml" );

}

//=============================================================================

StatusCode TupleToolApplyIsolation::initialize() {
  if( ! TupleToolBase::initialize() ) return StatusCode::FAILURE;
  
  m_dva = Gaudi::Utils::getIDVAlgorithm ( contextSvc() ) ;
  if (0==m_dva) return Error("Couldn't get parent DVAlgorithm", 
                             StatusCode::FAILURE);
  m_dist       = tool<IDistanceCalculator>("LoKi::DistanceCalculator",this);
  m_p2mcAssoc = tool<IParticle2MCAssociator>("DaVinciSmartAssociator", this);
  if( !m_dist ){
    Error("Unable to retrieve the IDistanceCalculator tool");
    return StatusCode::FAILURE;
  }
  m_pvReFitter = tool<IPVReFitter>("AdaptivePVReFitter", this );
  m_pVertexFit= m_dva->vertexFitter();
  //m_pVertexFit= tool<ITrackVertexer>

  if( !m_pVertexFit ){
    Error("Unable to retrieve the IVertexFit tool");
    return StatusCode::FAILURE;
  }
  
  
  m_Reader = new TMVA::Reader( "!Silent" );
  m_Reader->AddSpectator( "Track_TYPE",&type);
  m_Reader->AddVariable( "Track_MINIPCHI2",&minipchi2);
  m_Reader->AddVariable( "Track_PT",&pt);
  m_Reader->AddVariable( "Track_OPENING",&opening);
  m_Reader->AddVariable( "Track_IPCHI2",&chi2);
  m_Reader->AddVariable( "Track_FLIGHT",&newfdchi2);
  m_Reader->AddVariable( "Track_DELTAFLIGHT",&deltafd);
  //dummy variables
  //m_Reader->AddVariable( "Bplus_PT",&dummy);
  //m_Reader->AddVariable( "Dst_PT",&Dst_PT);
  //m_Reader->AddVariable( "Bplus_ENDVERTEX_CHI2",&vertexchi2);
  //m_Reader->AddVariable( "Dst_ENDVERTEX_CHI2",&dummy);
  //m_Reader->AddVariable( "Dst_FDCHI2_OWNPV",&dummy);
  //m_Reader->AddVariable( "log(1-D_DIRA_OWNPV)",&dummy);
  //m_Reader->AddVariable( "log(1-Bplus_DIRA_OWNPV)",&dummy);
  
  

  
   

  //reader->AddVariable("Bplus_PT",&Bplus_PTf);
  //reader->AddVariable("Dst_PT",&Dst_PTf);
  m_Reader->BookMVA( "BDT method", m_weightsName );
  
    if( !m_Reader ){
    Error("Unable to retrieve the IVertexFit tool");
    return StatusCode::FAILURE;
  }
  return StatusCode::SUCCESS;
}

//=============================================================================

StatusCode TupleToolApplyIsolation::fill( const Particle* mother
                                    , const Particle* P
                                    , const std::string& head
                                    , Tuples::Tuple& tuple )
{
  
  const std::string prefix=fullName(head);
  Assert( P && mother && m_dist
          , "This should not happen, you are inside TupleToolVtxIsoln.cpp :(" );

  
  bool test=true;
  float charge = 0;
  float type = 0;
  float px = 0;
  float py = 0;
  float pz = 0;
  float pe = 0;
  float pidk = 0;
  float nnk = 0;
  float nnpi = 0;
  float nng = 0;
  float ismuon = 0;

  /*
  const LHCb::Vertex* vtx;
  if (P->isBasicParticle()){
    vtx = mother->endVertex(); 
  }
  else{
    vtx = P->endVertex();
    
  }
  debug()<<"vertex for P, ID " <<P->particleID().pid()<<" = " <<vtx<<" at "<<vtx->position()<<  endmsg;
  if( !vtx ){
    Error("Can't retrieve the  vertex for " + prefix );
    return StatusCode::FAILURE;
  }
  */  
  std::vector<const LHCb::Track*> daughtertracks;
  daughtertracks.clear();
  LHCb::Particle::ConstVector source;
  LHCb::Particle::ConstVector target;
  LHCb::Particle::ConstVector finalStates;
  LHCb::Particle::ConstVector parts2Vertex;
  LHCb::Particle::ConstVector parts2VertexD;
  double maxbdt = -2;
  double bdt2 = -2;
  double bdt3 = -2;
  const LHCb::Particle* maxpart;
  const LHCb::Particle* part2;
  const LHCb::Particle* part3;
  vertexchi2 = P->endVertex()->chi2();
  parts2Vertex.clear();
  parts2VertexD.clear();
  
  //   const LHCb::Particle* prefix = P;
  if (P->isBasicParticle()){
    source.push_back(mother);
  }
  else{
    source.push_back(P);
  }
  LHCb::Vertex dv2;

  do {
    target.clear();
    for(LHCb::Particle::ConstVector::const_iterator isource = source.begin(); 
        isource != source.end(); isource++){
      
      if(!((*isource)->daughters().empty())){
        
        LHCb::Particle::ConstVector tmp = (*isource)->daughtersVector();
        
        for( LHCb::Particle::ConstVector::const_iterator itmp = tmp.begin(); 
             itmp!=tmp.end(); itmp++){
          target.push_back(*itmp);
          // Add the final states, i.e. particles with proto and ignoring gammas
          if((*itmp)->proto() && 22 != (*itmp)->particleID().pid()){
           	finalStates.push_back(*itmp);
           	daughtertracks.push_back((*itmp)->proto()->track());
           	if((*itmp)->particleID().abspid() == 413) Dst_PT = (*itmp)->pt();
           }
        }
      } // if endVertex
    } // isource
    source = target;
  } while(target.size() > 0);
  if (msgLevel(MSG::DEBUG)) debug() << "Final states size= " <<  finalStates.size()  << endreq;
  //warning() << " D VERTEX CHI2 " << dv2.chi2() << " NDOF " << dv2.nDoF() << endreq;
  //warning() << "DAUGHTER SIZE " << daughtertracks.size() << endreq;
  //default gives best tracks (why would you default to anything less than the best?)  


  

  LHCb::Vertex v;
  //double chi2ndof = 0;//oldvtx->chi2();
  //int ndof = 0;//oldvtx->nDoF();
  
  
    if (P->isBasicParticle()){
    parts2Vertex.push_back(P);
  }
  else{
    parts2Vertex = P->daughtersVector();
    StatusCode sc = m_pVertexFit->fit(v,parts2Vertex);
  }

  //hack due to lack of programming skill
  //v=new LHCb::Vertex(v2);
  //********Loop through tracks************
  
  //number below am IPCHI2 threshold isnt that useful, will probably remove it 

  
  LHCb::Particle::ConstVector theParts;
    

       
  for(std::vector<std::string>::iterator i = m_inputParticles.begin();
      i !=m_inputParticles.end(); ++i){

    if (!exist<LHCb::Particle::Range>(*i+"/Particles")){
      if (msgLevel(MSG::DEBUG)) debug() << "No particles at " << *i << " !!!!!" << endreq;
      continue;
    }

    LHCb::Particle::Range parts = get<LHCb::Particle::Range>(*i+"/Particles");
    if (msgLevel(MSG::DEBUG)) debug() << "Getting particles from " << *i
                                      << " with " << (parts).size() << " particles" << endreq;
    //warning() << "Getting particles from " << *i
    //                                  << " with " << (parts).size() << " particles" << endreq;
    for(LHCb::Particle::Range::const_iterator iparts = (parts).begin();
        iparts != (parts).end(); ++iparts)
    {
    const LHCb::Particle* part = (*iparts);
    
    
    //if(isTrackInDecay(part->proto()->track(),daughtertracks)) warning() << "FOUND DAUGHTER TRACK" << endreq;
    if(part->proto()->track()->type() < 5 && !isTrackInDecay(part->proto()->track(),daughtertracks)){
    LHCb::Vertex vtxWithExtraTrack;
    parts2Vertex.push_back(*iparts);
    StatusCode sc3 = m_pVertexFit->fit (vtxWithExtraTrack,parts2Vertex);
    parts2Vertex.pop_back();
    opening = getopening(part->proto()->track(),P);
    minipchi2 = getminipchi(part);
    newfdchi2 = getfdchi2(part->proto()->track(),vtxWithExtraTrack);
    oldfdchi2 = getfdchi2(part->proto()->track(),v);
    ghostprob = part->proto()->track()->ghostProbability();
    trackchi2 = part->proto()->track()->chi2PerDoF();
    deltafd = log10(newfdchi2-oldfdchi2) - 7;
    type = part->proto()->track()->type();
    if (newfdchi2-oldfdchi2 < 0) deltafd = deltafd * -1.;
    newfdchi2=log10(newfdchi2);
    if(part->proto()->track()->type() == 1) pt = part->proto()->track()->momentum().z();
    else pt = part->proto()->track()->pt();



    if(ghostprob > 0.5){continue;}


    //warning() << "type " << type << " opening " << opening << " pt " << pt << endreq;

    if(part->proto()->track()->type() == 3 && !(opening > 0.994 )){continue;}
    if(part->proto()->track()->type() == 4 && !(opening > 0.98)){continue;}
    if(part->proto()->track()->type() == 1 && !(opening > 0.98)){continue;}    
    //if(track->info(LHCb::Track::CloneDist, -1.) > 0){continue;}
    StatusCode sc = StatusCode::SUCCESS;
    double tmpip, tmpchi2;
    StatusCode dump = m_dist->distance((const LHCb::Particle *) part,(const LHCb::Vertex *)&v,tmpip,tmpchi2);
    chi2=tmpchi2;
	//StatusCode dump2 = m_dist->distance((const LHCb::Particle *) part,(const LHCb::Vertex *)vd,D_ip,D_chi2);
	    
	    if(chi2 < 50){
	        dummy = 4000;            float bdtval = m_Reader->EvaluateMVA( "BDT method" );
            //warning() << "bdtval " << bdtval << " old maxbdt " << maxbdt << endreq;	        if (bdtval > maxbdt) {
	            bdt3 = bdt2;
	            bdt2 = maxbdt;
	            maxbdt = bdtval;
	            part3=part2;
	            part2=maxpart;
	            maxpart = part;
	        }
	        else if (bdtval > bdt2) {
	            bdt3 = bdt2;
	            bdt2 = bdtval;
	            part3=part2;
	            part2=part;
	        }
	        else if (bdtval > bdt3) {
	            bdt3 = bdtval;
	            part3=part;
	        }
	        //warning() << "new max bdtval " << maxbdt << endreq;
	    }
	    

    }
  } // end particles loop
 }//end particle types loop

          
  if(maxbdt > -1){
	            pe = maxpart->momentum().E();
	            px = maxpart->momentum().Px();
	            py = maxpart->momentum().Py();
	            pz = maxpart->momentum().Pz();
	            pidk = maxpart->proto()->info(LHCb::ProtoParticle::CombDLLk,-1000);
	            nnk = maxpart->proto()->info(LHCb::ProtoParticle::ProbNNk,-1000);
	            nnpi = maxpart->proto()->info(LHCb::ProtoParticle::ProbNNpi,-1000);
	            nng = maxpart->proto()->info(LHCb::ProtoParticle::ProbNNghost,-1000);
	            if(maxpart->proto()->track()->type() ==1){charge=0;}
	            else{charge=maxpart->proto()->track()->charge();}        
		    type = maxpart->proto()->track()->type();
	            const MuonPID * muonPID = maxpart->proto()->muonPID();
		    ismuon = muonPID ? muonPID->IsMuon() : false;
          }

  tuple->column(prefix + "_ISOLATION_BDT" + m_outputSuffix,  maxbdt );
  tuple->column(prefix + "_ISOLATION_CHARGE" + m_outputSuffix,  charge );
  tuple->column(prefix + "_ISOLATION_Type" + m_outputSuffix,  type );
  tuple->column(prefix + "_ISOLATION_PE" + m_outputSuffix,  pe );
  tuple->column(prefix + "_ISOLATION_PX" + m_outputSuffix,  px );
  tuple->column(prefix + "_ISOLATION_PY" + m_outputSuffix,  py );
  tuple->column(prefix + "_ISOLATION_PZ" + m_outputSuffix,  pz );
  tuple->column(prefix + "_ISOLATION_PIDK" + m_outputSuffix,  pidk );
  tuple->column(prefix + "_ISOLATION_NNk" + m_outputSuffix,  nnk );
  tuple->column(prefix + "_ISOLATION_NNpi" + m_outputSuffix,  nnpi );
  tuple->column(prefix + "_ISOLATION_IsMuon" + m_outputSuffix,  ismuon );
  tuple->column(prefix + "_ISOLATION_NNghost" + m_outputSuffix,  nng );
  
  if(bdt2 > -1){
	            pe = part2->momentum().E();
	            px = part2->momentum().Px();
	            py = part2->momentum().Py();
	            pz = part2->momentum().Pz();
	            pidk = part2->proto()->info(LHCb::ProtoParticle::CombDLLk,-1000);
	            nnk = part2->proto()->info(LHCb::ProtoParticle::ProbNNk,-1000);
	            nnpi = part2->proto()->info(LHCb::ProtoParticle::ProbNNpi,-1000);
	            nng = part2->proto()->info(LHCb::ProtoParticle::ProbNNghost,-1000);
	            if(part2->proto()->track()->type() ==1){charge=0;}
	            else{charge=part2->proto()->track()->charge();}        
		    type = part2->proto()->track()->type();
	            const MuonPID * muonPID = part2->proto()->muonPID();
		    ismuon = muonPID ? muonPID->IsMuon() : false;        
          }

  tuple->column(prefix + "_ISOLATION_BDT2" + m_outputSuffix,  bdt2 );
  tuple->column(prefix + "_ISOLATION_CHARGE2" + m_outputSuffix,  charge );
  tuple->column(prefix + "_ISOLATION_Type2" + m_outputSuffix,  type );
  tuple->column(prefix + "_ISOLATION_PE2" + m_outputSuffix,  pe );
  tuple->column(prefix + "_ISOLATION_PX2" + m_outputSuffix,  px );
  tuple->column(prefix + "_ISOLATION_PY2" + m_outputSuffix,  py );
  tuple->column(prefix + "_ISOLATION_PZ2" + m_outputSuffix,  pz );
  tuple->column(prefix + "_ISOLATION_NNk2" + m_outputSuffix,  nnk );
  tuple->column(prefix + "_ISOLATION_NNpi2" + m_outputSuffix,  nnpi );
  tuple->column(prefix + "_ISOLATION_IsMuon2" + m_outputSuffix,  ismuon );
  tuple->column(prefix + "_ISOLATION_NNghost2" + m_outputSuffix,  nng );
  
  if(bdt3 > -1){
	            pe = part3->momentum().E();
	            px = part3->momentum().Px();
	            py = part3->momentum().Py();
	            pz = part3->momentum().Pz();
	            pidk = part3->proto()->info(LHCb::ProtoParticle::CombDLLk,-1000);
	            nnk = part3->proto()->info(LHCb::ProtoParticle::ProbNNk,-1000);
	            nnpi = part3->proto()->info(LHCb::ProtoParticle::ProbNNpi,-1000);
	            nng = part3->proto()->info(LHCb::ProtoParticle::ProbNNghost,-1000);
	            if(part3->proto()->track()->type() ==1){charge=0;}
	            else{charge=part3->proto()->track()->charge();}        
		    type = part3->proto()->track()->type();  
	            const MuonPID * muonPID = part3->proto()->muonPID();
		    ismuon = muonPID ? muonPID->IsMuon() : false;
          }

  tuple->column(prefix + "_ISOLATION_BDT3" + m_outputSuffix,  bdt3 );
  tuple->column(prefix + "_ISOLATION_CHARGE3" + m_outputSuffix,  charge );
  tuple->column(prefix + "_ISOLATION_Type3" + m_outputSuffix,  type );
  tuple->column(prefix + "_ISOLATION_PE3" + m_outputSuffix,  pe );
  tuple->column(prefix + "_ISOLATION_PX3" + m_outputSuffix,  px );
  tuple->column(prefix + "_ISOLATION_PY3" + m_outputSuffix,  py );
  tuple->column(prefix + "_ISOLATION_PZ3" + m_outputSuffix,  pz );
  tuple->column(prefix + "_ISOLATION_NNk3" + m_outputSuffix,  nnk );
  tuple->column(prefix + "_ISOLATION_NNpi3" + m_outputSuffix,  nnpi );
  tuple->column(prefix + "_ISOLATION_IsMuon3" + m_outputSuffix,  ismuon );
  tuple->column(prefix + "_ISOLATION_NNghost3" + m_outputSuffix,  nng );
  
  
	    

  return StatusCode(test);
}

//=========================================================================
//  
//=========================================================================
const Vertex* TupleToolApplyIsolation::originVertex( const Particle* top
                                               , const Particle* P ) const {
  if( top == P || P->isBasicParticle() ) return 0;
  
  const SmartRefVector< LHCb::Particle >& dau = top->daughters ();
  if( dau.empty() ){
    // if (msgLevel(MSG::DEBUG)) debug() << " Particle has no daughters! "  << endreq;
    return 0;
  }
  
  SmartRefVector< LHCb::Particle >::const_iterator it;
  for( it = dau.begin(); dau.end()!=it; ++it ){
    if( P == *it ){ // I found the daughter
      return top->endVertex();
    }
  }
  
  // vertex not yet found, get deeper in the decay:
  for( it = dau.begin(); dau.end()!=it; ++it ){
    if( P != *it && !(*it)->isBasicParticle() ){
      const Vertex* vv = originVertex( *it, P );
      if( vv ) return vv;
    }
  }
  return 0;
}


//=============================================================================
// Check if the track is already in the decay
//=============================================================================
bool TupleToolApplyIsolation::isTrackInDecay(const LHCb::Track* track, std::vector<const LHCb::Track*> daughters){
  bool isInDecay = false;
        //loop over daughters
        for(std::vector<const LHCb::Track*>::iterator it = daughters.begin(); it != daughters.end(); ++it){
        const LHCb::Track* itrack = (*it);
        if(itrack){
         if(itrack == track){
          if ( msgLevel(MSG::DEBUG) ) debug() << "Track is in decay, skipping it" << endmsg;
          isInDecay = true;
          }
         }        
  }//end daughter loop
      
  return isInDecay;
}

//=============================================================================
// MINIPCHI2 for a track
//=============================================================================
double TupleToolApplyIsolation::getminipchi(const LHCb::Particle* track){

double minchi2 = -1 ;
const RecVertex::Range PV = m_dva->primaryVertices();
  if ( !PV.empty() ){
    for ( RecVertex::Range::const_iterator pv = PV.begin() ; pv!=PV.end() ; ++pv){
      double ip, chi2;
      StatusCode test2 = m_dist->distance ( (const LHCb::Particle*)track, *pv, ip, chi2 );
        if ((chi2<minchi2) || (minchi2<0.)) 
        {
	LHCb::RecVertex newPV(**pv);
        StatusCode scfit = m_pvReFitter->remove(track, &newPV);
	LHCb::RecVertex* newPVPtr = (LHCb::RecVertex*)&newPV;
	test2 = m_dist->distance ((LHCb::Particle *)track, (LHCb::VertexBase*) newPVPtr, ip, chi2 );
        minchi2 = chi2 ;       
        } 
    }

}

return minchi2;
}

double TupleToolApplyIsolation::getfdchi2(const LHCb::Track* track, LHCb::Vertex Vtx){

double minchi2 = -1 ;
double fdchi2 = -1;
double fd;
const RecVertex::Range PV = m_dva->primaryVertices();
  if ( !PV.empty() ){
    for ( RecVertex::Range::const_iterator pv = PV.begin() ; pv!=PV.end() ; ++pv){
      double ip, chi2;
      StatusCode test2 = m_dist->distance ( (const LHCb::Track *)track, *pv, ip, chi2 );
        if ((chi2<minchi2) || (minchi2<0.)) 
        {
          minchi2 = chi2 ;
	  StatusCode test2 = m_dist->distance ( *pv, &Vtx, fd, fdchi2 ); 
        } 
    }

}

return fdchi2;
}



//=============================================================================
// Opening angle for a track and particle
//=============================================================================
double TupleToolApplyIsolation::getopening(const LHCb::Track* track,const  LHCb::Particle* P){
    Gaudi::XYZVector A = P->momentum().Vect();
    Gaudi::XYZVector B = track->momentum();
    double cosopening = A.Dot( B ) / std::sqrt( A.Mag2()*B.Mag2() );
    return cosopening;
}
