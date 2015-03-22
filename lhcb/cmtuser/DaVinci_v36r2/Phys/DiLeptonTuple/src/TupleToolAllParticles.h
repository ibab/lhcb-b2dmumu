// $Id: TupleToolAllParticles.h,v 1.7 2010-09-09 12:22:42 pkoppenb Exp $
#ifndef TUPLETOOLALLTRACKS_H
#define TUPLETOOLALLTRACKS_H 1

// Include files
// from Gaudi
#include "DecayTreeTupleBase/TupleToolBase.h"
#include "Kernel/IEventTupleTool.h"            // Interface
#include "Kernel/IRelatedPVFinder.h"
#include <Kernel/IDistanceCalculator.h>

class ITrackExtrapolator;
class ITrackStateProvider;
class IParticle2MCAssociator;

/** @class TupleToolAllParticles TupleToolAllParticles.h
 *
 *  Fills minimal info for a set of particles in an EventTupleTool
 *  Warning: This will explode your tuple.
 *
 *  @author Patrick Koppenburg
 *  @date   2013-02-28
 */
class TupleToolAllParticles : public TupleToolBase, 
                              virtual public IEventTupleTool
{

public:

  /// Standard constructor
  TupleToolAllParticles( const std::string& type,
                      const std::string& name,
                      const IInterface* parent);

  virtual ~TupleToolAllParticles( ); ///< Destructor
  virtual StatusCode fill( Tuples::Tuple& );///< Fill tuple
  StatusCode initialize();///< init
  
private:
  int TID(const LHCb::Particle* p);
  std::string m_location ;
  unsigned int m_max ;
  const IDistanceCalculator* m_dist;
  IRelatedPVFinder* m_pvf;
  bool m_allChi2 ; ///< fill all chi2
  const ITrackStateProvider* m_extrapolator; ///<pointer tot the track extrapolator
  std::vector<IParticle2MCAssociator*> m_p2mcAssocs;  
  std::vector<std::string> m_p2mcAssocTypes;
  

};

#endif // TUPLETOOLALLTRACKS_H

