/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELANNULUSCLUSTERIMPL_H
#define EUTELANNULUSCLUSTERIMPL_H

#include "EUTelGenericSparseClusterImpl.h"
#include "EUTelAnnulusPixel.h"

namespace eutelescope {
  
  //! Implementation of a cluster made of sparsified pixels
  class EUTelAnnulusClusterImpl : public EUTelGenericSparseClusterImpl<EUTelAnnulusPixel> {
    
  public:
    //! Default constructor
    EUTelAnnulusClusterImpl(IMPL::TrackerDataImpl* data);
    
    //! Destructor
    virtual ~EUTelAnnulusClusterImpl();
    
    void getClusterGeomInfo(float& FxCoord, float& FyCoord, float& xPos, float& yPos, float& rPos, float& phiPos, float& rSize, float& phiSize) const;
    
    void getGeometricCenterOfGravity(float& xCoG, float& yCoG) const;
    
    void getAnnulusCenterOfGravity(float& FxCoord, float& FyCoord, float& xCoG, float& yCoG) const;
    
    void getHitCovMatrix(double& Rres, double& Thetares, double (&cov)[4] ) const;
    
  private:
    DISALLOW_COPY_AND_ASSIGN(EUTelAnnulusClusterImpl)
      };
}
#endif
