// $Id: TupleToolAllParticles.cpp,v 1.6 2010-09-09 12:22:42 pkoppenb Exp $
// Include files

// from Gaudi
#include "GaudiKernel/ToolFactory.h"
#include "Event/Track.h"
#include "Event/RecSummary.h"
#include "Event/VertexBase.h"
#include "Event/RecVertex.h"
#include "Kernel/DefaultDVToolTypes.h"

// local
#include "TupleToolAllParticles.h"
#include "TrackInterfaces/ITrackExtrapolator.h"
#include "TrackInterfaces/ITrackStateProvider.h"
#include "Kernel/IParticle2MCAssociator.h"

//-----------------------------------------------------------------------------
// Implementation file for class : TupleToolAllParticles
//
// 2009-02-11 : Patrick Koppenburg
//-----------------------------------------------------------------------------

//=============================================================================
// Standard constructor, initializes variables
//=============================================================================
TupleToolAllParticles::TupleToolAllParticles( const std::string& type,
                                        const std::string& name,
                                        const IInterface* parent )
  : TupleToolBase ( type, name, parent )
  , m_dist(0)
  , m_pvf(0)
  , m_extrapolator(0)
{
  declareInterface<IEventTupleTool>(this);

  declareProperty("Location", m_location = "/Event/Phys/StdAllNoPIDsPions/Particles", "Location of particles");
  declareProperty("Max", m_max = 1000, "Max size of array");
  declareProperty("AllChi2", m_allChi2 = false , "Fill all Chi2?");
  // MC associators to try, in order
  m_p2mcAssocTypes.push_back( "DaVinciSmartAssociator" );
  m_p2mcAssocTypes.push_back( "MCMatchObjP2MCRelator"  );
  declareProperty( "IP2MCPAssociatorTypes", m_p2mcAssocTypes );
}

//=============================================================================
// Destructor
//=============================================================================
TupleToolAllParticles::~TupleToolAllParticles() {}

//=============================================================================
// Destructor
//=============================================================================
StatusCode TupleToolAllParticles::initialize(){
  StatusCode sc =  TupleToolBase::initialize();
  if ( sc.isFailure() ) return sc;
  info() << "Will fill all particles from " << m_location << ". WARNING! This will make your tuple HUGE!" << endmsg ;
  m_pvf = tool<IRelatedPVFinder>( DaVinci::DefaultTools::PVRelator, this );
  if (!m_pvf) return Error("Couldn't get PVRelator");

  m_dist = tool<IDistanceCalculator>( DaVinci::DefaultTools::Distance, this );
  if ( !m_dist )
  {
    return Error("Unable to retrieve the IDistanceCalculator tool");
  }
  m_extrapolator=tool<ITrackStateProvider>("TrackStateProvider",this);
  // the MC associators
  for ( std::vector<std::string>::const_iterator iMCAss = m_p2mcAssocTypes.begin();
        iMCAss != m_p2mcAssocTypes.end(); ++iMCAss )
  {
    m_p2mcAssocs.push_back( tool<IParticle2MCAssociator>(*iMCAss,this) );
  }

  return sc ;

}

//=============================================================================
StatusCode TupleToolAllParticles::fill( Tuples::Tuple& tup )
{
  const std::string prefix = fullName();

  LHCb::Particle::Range parts ;
  if (exist<LHCb::Particle::Range>(m_location)){
    parts = get<LHCb::Particle::Range>(m_location);
  } else return Warning("Nothing found at "+m_location,StatusCode::SUCCESS,1);
  // Fill the tuple
  bool test = true;
  
  std::vector<double> PX,PY,PZ,CHI2,PIDK,PIDMU,PIDP,PIDE,IPCHI2,IP,GHOST,PIDPI,PIDGH;
  std::vector<double> CHI2Velo, CHI2T, CHI2Match, XATT1, YATT1;
  std::vector<int> NDOFVelo, NDOFT;
  std::vector<bool> ISMUON;
  std::vector<int> PID, TCAT, TPID;
  
  const LHCb::VertexBase* aPV = NULL;
  double ip, chi2;
  for ( LHCb::Particle::Range::const_iterator p = parts.begin() ; p!=parts.end() ; ++p){
    if (PID.size()==m_max) {
      Warning("Reached maximum size of array",StatusCode::SUCCESS).ignore();
      break ;
    }
    PID.push_back((*p)->particleID().pid());
    TPID.push_back(TID(*p));
    PX.push_back((*p)->momentum().x());
    PY.push_back((*p)->momentum().y());
    PZ.push_back((*p)->momentum().z());
    const LHCb::ProtoParticle* pp = (*p)->proto();
    if ( pp ) {
      if ( pp->track()) {
        CHI2.push_back( pp->track()->chi2PerDoF()) ;
        GHOST.push_back( pp->track()->ghostProbability()) ;
        if ( m_allChi2){
          CHI2Velo.push_back( pp->track()->info(LHCb::Track::FitVeloChi2,-1.)) ;
          NDOFVelo.push_back( pp->track()->info(LHCb::Track::FitVeloNDoF,-1.)) ;
          CHI2T.push_back( pp->track()->info(LHCb::Track::FitTChi2,-1.)) ;
          NDOFT.push_back( pp->track()->info(LHCb::Track::FitTNDoF,-1.)) ;
          CHI2Match.push_back( pp->track()->info(LHCb::Track::FitMatchChi2,-1.)) ;          
          bool hasIT = false ;
          bool hasOT = false ;
          bool hasTT = false ;
          const std::vector<LHCb::LHCbID> mes = pp->track()->lhcbIDs();
          for ( std::vector< LHCb::LHCbID >::const_iterator m = mes.begin() ; m!=mes.end() ; ++m){
            if (m->isIT()) hasIT=true ;
            if (m->isOT()) hasOT=true ;
            if (m->isTT()) hasTT=true ;
          }
          int cat = (hasIT?1:0); // IT is one
          if (hasOT) cat += 2;   // OT is two. overlap is 3.
          if (!hasTT) cat = -cat; // when there's no TT it's -cat
          TCAT.push_back( cat ) ;          
          LHCb::State aState=LHCb::State();
          StatusCode sc=m_extrapolator->stateFromTrajectory(aState,*(pp->track()),7800.);
          if (!sc) return sc;
          XATT1.push_back(aState.x());
          YATT1.push_back(aState.y());
        }
      }
      /*
      PIDK.push_back( pp->info(LHCb::ProtoParticle::CombDLLk,-1000));
      PIDMU.push_back( pp->info(LHCb::ProtoParticle::CombDLLmu,-1000));
      PIDP.push_back( pp->info(LHCb::ProtoParticle::CombDLLp,-1000));
      PIDE.push_back( pp->info(LHCb::ProtoParticle::CombDLLe,-1000));
      */
      PIDK.push_back( pp->info(LHCb::ProtoParticle::ProbNNk,-1000));
      PIDMU.push_back( pp->info(LHCb::ProtoParticle::ProbNNmu,-1000));
      PIDP.push_back( pp->info(LHCb::ProtoParticle::ProbNNp,-1000));
      PIDE.push_back( pp->info(LHCb::ProtoParticle::ProbNNe,-1000));
      PIDPI.push_back( pp->info(LHCb::ProtoParticle::ProbNNpi,-1000));
      PIDGH.push_back( pp->info(LHCb::ProtoParticle::ProbNNghost,-1000));
    }
    aPV = m_pvf->relatedPV ( *p, LHCb::RecVertexLocation::Primary  );
    if (aPV){
      m_dist->distance ( *p, aPV, ip, chi2 ).ignore();
      IP.push_back(ip);
      IPCHI2.push_back(chi2);      
    } else {
      IP.push_back(-999.);
      IPCHI2.push_back(-999.);
    }
  }
  unsigned int siz = PID.size();
  test &= ( CHI2.empty() || CHI2.size()==siz ) ; // check
  test &= ( PIDK.empty() || PIDK.size()==siz ) ; // check
  if (!test){
    err() << "Inconsistent array sizes " << siz << " and " << CHI2.size() 
          << " and " << PIDK.size() << endmsg ;
    return StatusCode(test);
  }

  test &= tup->farray( prefix+"PARTS_PID", PID, prefix+"nParts", m_max );
  test &= tup->farray( prefix+"PARTS_TPID", TPID, prefix+"nParts", m_max );
  test &= tup->farray( prefix+"PARTS_PX", PX, prefix+"nParts", m_max );
  test &= tup->farray( prefix+"PARTS_PY", PY, prefix+"nParts", m_max );
  test &= tup->farray( prefix+"PARTS_PZ", PZ, prefix+"nParts", m_max );
  if ( !CHI2.empty() || 0==siz ){
    test &= tup->farray( prefix+"PARTS_CHI2", CHI2, prefix+"nParts", m_max );
    test &= tup->farray( prefix+"PARTS_GHOST", GHOST, prefix+"nParts", m_max );
  }
  if ( !CHI2Velo.empty() || 0==siz ){
    test &= tup->farray( prefix+"PARTS_CHI2VELO", CHI2Velo, prefix+"nParts", m_max );
    test &= tup->farray( prefix+"PARTS_NDOFVELO", NDOFVelo, prefix+"nParts", m_max );
    test &= tup->farray( prefix+"PARTS_CHI2T", CHI2T, prefix+"nParts", m_max );
    test &= tup->farray( prefix+"PARTS_NDOFT", NDOFT, prefix+"nParts", m_max );
    test &= tup->farray( prefix+"PARTS_CHI2MATCH", CHI2Match, prefix+"nParts", m_max );
    test &= tup->farray( prefix+"PARTS_TCAT", TCAT, prefix+"nParts", m_max );
    test &= tup->farray( prefix+"PARTS_XATT1", XATT1, prefix+"nParts", m_max );
    test &= tup->farray( prefix+"PARTS_YATT1", YATT1, prefix+"nParts", m_max );
  }
  if ( !PIDK.empty() || 0==siz ){
    test &= tup->farray( prefix+"PARTS_PIDK", PIDK, prefix+"nParts", m_max );
    test &= tup->farray( prefix+"PARTS_PIDMU", PIDMU, prefix+"nParts", m_max );
    test &= tup->farray( prefix+"PARTS_PIDP", PIDP, prefix+"nParts", m_max );
    test &= tup->farray( prefix+"PARTS_PIDE", PIDE, prefix+"nParts", m_max );
    test &= tup->farray( prefix+"PARTS_PIDPI", PIDPI, prefix+"nParts", m_max );
    test &= tup->farray( prefix+"PARTS_PIDGH", PIDGH, prefix+"nParts", m_max );
  }
  test &= tup->farray( prefix+"PARTS_IP", IP, prefix+"nParts", m_max );
  test &= tup->farray( prefix+"PARTS_IPCHI2", IPCHI2, prefix+"nParts", m_max );
  
  return StatusCode(test);
}
//============================================================================
int TupleToolAllParticles::TID(const LHCb::Particle* P)
{
  Assert( !m_p2mcAssocs.empty(),
          "The DaVinci smart associator(s) have not been initialized!");
  const LHCb::MCParticle* mcp(NULL);
  if ( P )
  {
    //assignedPid = P->particleID().pid();
    if (msgLevel(MSG::VERBOSE)) verbose() << "Getting related MCP to " << P << endmsg ; 
    for ( std::vector<IParticle2MCAssociator*>::const_iterator iMCAss = m_p2mcAssocs.begin();
          iMCAss != m_p2mcAssocs.end(); ++iMCAss ){
      mcp = (*iMCAss)->relatedMCP(P); 
      if ( mcp ) break;
    }
    if (msgLevel(MSG::VERBOSE)) verbose() << "Got mcp " << mcp << endmsg ;
  }
  // pointer is ready, prepare the values:
  if( mcp ) return mcp->particleID().pid(); 
  else return 0 ;
}


// Declaration of the Tool Factory
DECLARE_TOOL_FACTORY( TupleToolAllParticles )
