#include "ITkAnnulusStrip.h"
#include <iostream>

namespace eutelescope {
  namespace geo {

    ITkAnnulusStrip::ITkAnnulusStrip( int xPixel, int yPixel, double xSize, double ySize, double zSize, double pitchPhi, double stereoAngle, 
				      double rmin, double rmax, double rCentre, int order, double radLength)
      : EUTelGenericPixGeoDescr( xSize, ySize, zSize,
				 0, xPixel-1, 0, yPixel-1,
				 radLength), _pitchPhi(pitchPhi), _stereoAngle(stereoAngle), _rMin(rmin), _rMax(rmax), _rCentre(rCentre), _stripOrder(order), _stripLength(rmax-rmin)
    {
      static const double PI = 3.141592653589793;
      static const double DEG = 180. / PI;

      // Create the material for the sensor
      matSi =  new TGeoMaterial("Si", 28.0855, 14.0, 2.33, -_radLength, 45.753206);
      Si = new TGeoMedium("AnnulusStripSilicon", 1, matSi);

      // Create a plane for the sensitive area
      plane = _tGeoManager->MakeBox("sns_annulusstrip", Si, xSize/2, ySize/2, zSize/2);

      // Create a strip for the sensitive area
      rowstrip = _tGeoManager->MakeTubs("strip_annulusstrip" , Si, rmin, rmax, zSize/2, 90+(-pitchPhi/2)*DEG, 90+(pitchPhi/2)*DEG);

      // Calculate strip position and coordinates (according to XXX)
      Double_t phi_i, b, r, c, x, y, theta;

      theta = pitchPhi*( (xPixel-1)/2.0 ) + stereoAngle + PI/2;

      //strip index counting from right to left
      if(order == -1) {
	for( int i = xPixel-1; i >=0; i-- ) {
	  TGeoCombiTrans* transform = new TGeoCombiTrans(0,0,0, new TGeoRotation("rot",0,0,0) );

	  phi_i = (i -(xPixel-1)/2.0)*pitchPhi;  
	  b = -2*(2*rCentre*sin(stereoAngle/2))*sin(stereoAngle/2+phi_i);
	  c = pow((2*rCentre*sin(stereoAngle/2)),2)-pow(rmin,2);
	  r = 0.5*(-b+sqrt(pow(b,2)-4*c));
	  y = r*cos(phi_i+stereoAngle) - rCentre*cos(stereoAngle);
	  x = -r*sin(phi_i+stereoAngle) + rCentre*sin(stereoAngle);

	  //create first transform (rotate to get correct angle of strip)
	  transform->RotateZ(theta*DEG-90);
	  transform->SetTranslation(x-rmin*cos(theta), y-rmin*sin(theta), 0);
 
	  //add the nodes
	  plane->AddNode(rowstrip, i+1, transform);

	  //get angle of next strip
	  theta -= pitchPhi;
	}	
      }
      //strip index counting from left to right
      else if(order == 1) {
	for( int i = 0; i <=xPixel-1; i++ ) {
	  TGeoCombiTrans* transform=new TGeoCombiTrans(0,0,0,new TGeoRotation("rot",0,0,0));

	  phi_i = ((xPixel -1 - i) - (xPixel-1)/2.0)*pitchPhi;  
	  b = -2*(2*rCentre*sin(stereoAngle/2))*sin(stereoAngle/2+phi_i);
	  c = pow((2*rCentre*sin(stereoAngle/2)),2)-pow(rmin,2);
	  r = 0.5*(-b+sqrt(pow(b,2)-4*c));
	  y = r*cos(phi_i+stereoAngle) - rCentre*cos(stereoAngle);
	  x = -r*sin(phi_i+stereoAngle) + rCentre*sin(stereoAngle);

	  //create first transform (rotate to get correct angle of strip)
	  transform->RotateZ(theta*DEG-90);
	  transform->SetTranslation(x-rmin*cos(theta), y-rmin*sin(theta), 0);

	  //add the nodes
	  plane->AddNode(rowstrip, i+1, transform);

	  //get angle of next strip
	  theta -= pitchPhi;
	}
      }
      else {
	std::cout << "Strip index order must be 1 or -1! Please check!" << std::endl;
      }
    }

    ITkAnnulusStrip::~ITkAnnulusStrip() {
      // delete matSi,
      // delete Si;
    }

    void ITkAnnulusStrip::createRootDescr(char const *planeVolume) {
      // Get the plane as provided by the EUTelGeometryTelescopeGeoDescription
      TGeoVolume *topplane = _tGeoManager->GetVolume(planeVolume);
      // Add the sensitive area to the plane
      topplane->AddNode(plane, 1, new TGeoTranslation(0,0,0) );
    }

    /*TODO: needs to be changed to include the information of rows and indice of pixels in x/y direction */
    std::string ITkAnnulusStrip::getPixName(int x, int y) {
      char buffer[100];
      // return path to the pixel, don't forget to shift indices by +1+
      if(x<512) {
	snprintf(buffer, 100, "/sns_annulusstrip_1/strip_annulusstrip_%d_%d", x+1, y+1);
      }
      else if(x>=512 && x<=_maxIndexX+256) {
	snprintf(buffer, 100, "/sns_annulusstrip_1/strip_annulusstrip_%d_%d", x+1-256, y+1);
      }
      else {
	snprintf(buffer, 100, "invalid index");
      }

      return std::string(buffer);
    }

    /*TODO*/ std::pair<int, int> ITkAnnulusStrip::getPixIndex(char const *) {
      return std::make_pair(0, 0);
    }

    void ITkAnnulusStrip::getPitchesAndStereo(double& pitchPhi,  double& stereoAngle){
      pitchPhi = _pitchPhi;
      stereoAngle = _stereoAngle;
    }

    void ITkAnnulusStrip::getRadii(double& rMin, double& rMax, double& rCentre){
      rMin = _rMin;
      rMax = _rMax;
      rCentre = _rCentre;
    }

    void ITkAnnulusStrip::getStripLength(double & stripLength){
      stripLength = _stripLength;
    }

    void ITkAnnulusStrip::getStripOrder(int & stripOrder){
      stripOrder = _stripOrder;
    }

  } // namespace geo
} // namespace eutelescope
