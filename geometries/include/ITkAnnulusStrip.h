#ifndef ITKANNULUSSTRIP_H
#define ITKANNULUSSTRIP_H

/** @class ITkAnnulusStrip
      * This class is the implementation of  @class EUTelGenericPixGeoDescr
      * for the ITkAnnulusStrip layout.
  */

// EUTELESCOPE
#include "EUTelGenericPixGeoDescr.h"

// ROOT
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TGeoVolume.h"

namespace eutelescope {
  namespace geo {

    class ITkAnnulusStrip : public EUTelGenericPixGeoDescr {

    public:
      ITkAnnulusStrip( int xPixel, int yPixel, double xSize, double ySize, double zSize, double pitchPhi, double stereoAngle, double rmin, double rmax, double rCentre, int order, double radLength);
      ~ITkAnnulusStrip();

      void createRootDescr(char const *);
      std::string getPixName(int, int);
      std::pair<int, int> getPixIndex(char const *);

      void getPitchesAndStereo(double&, double&);
      void getRadii(double&, double&, double&);
      void getStripLength(double&);
      void getStripOrder(int&);

    protected:
      TGeoMaterial *matSi;
      TGeoMedium *Si;
      TGeoVolume *plane,*rowstrip;
      double _pitchPhi, _stereoAngle;
      double _rMin, _rMax, _rCentre;
      int _stripOrder;
      double _stripLength;
    };

    extern "C" {
    EUTelGenericPixGeoDescr *maker();
    }

  } // namespace geo
} // namespace eutelescope

#endif // ANNULUSSTRIP_H
