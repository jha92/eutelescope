#include "EUTelAnnulusClusterImpl.h"
#include "EUTelAnnulusPixel.h"

using namespace eutelescope;

EUTelAnnulusClusterImpl::EUTelAnnulusClusterImpl(IMPL::TrackerDataImpl *data)
  : EUTelGenericSparseClusterImpl<EUTelAnnulusPixel>(data) {
}

EUTelAnnulusClusterImpl::~EUTelAnnulusClusterImpl() {
}

void EUTelAnnulusClusterImpl::getClusterGeomInfo(float& FxCoord, float& FyCoord, 
						 float& xPos, float& yPos, 
						 float& rPos, float& phiPos, 
						 float& rSize, float& phiSize) const {
  float rMin = std::numeric_limits<float>::max();
  float phiMin = rMin;
  float rMax = -rMin;
  float phiMax = rMax;
  float rMinBoundary = 0;
  float rMaxBoundary = 0;
  float phiMinBoundary = 0;
  float phiMaxBoundary = 0;

  // Loop over all the pixels in the cluster
  auto &pixelVec = getPixels();
  for(auto &pixel : pixelVec) {
  
    float rlength = pixel.getR();
    float ang = pixel.getAngle();
    
    if(ang < phiMin) {
      phiMin = ang;
      phiMinBoundary = pixel.getDphi();
    }
    if(ang > phiMax) {
      phiMax = ang;
      phiMaxBoundary = pixel.getDphi();
    }
    if(rlength < rMin) {
      rMin = rlength;
      rMinBoundary = pixel.getRmax() - pixel.getRmin();
    }
    if(rlength > rMax) {
      rMax = rlength;
      rMaxBoundary = pixel.getRmax() - pixel.getRmin();
    }
  }
  
  rSize = rMax - rMin + rMinBoundary*0.5 + rMaxBoundary*0.5;  //cluster size in radius r
  phiSize = phiMax - phiMin + phiMaxBoundary*0.5 + phiMinBoundary*0.5; //cluster width in angle phi
  
  rPos = rMax + rMaxBoundary*0.5 - rSize*0.5;   //cluster center in r
  phiPos = phiMax + phiMaxBoundary*0.5 - phiSize*0.5; //cluster center in phi

  xPos = rPos*cos(phiPos) + FxCoord;
  yPos = rPos*sin(phiPos) + FyCoord;
}

void EUTelAnnulusClusterImpl::getGeometricCenterOfGravity(float &xCoG, float &yCoG) const
{
  xCoG = 0;
  yCoG = 0;
  
  double totalCharge = 0;

  // Loop over all the pixels in the cluster
  auto &pixelVec = getPixels();
  for( auto &pixel : pixelVec ) {
    double curSignal = pixel.getSignal();
    xCoG += (pixel.getPosX()) * curSignal;
    yCoG += (pixel.getPosY()) * curSignal;
    totalCharge += curSignal;
  }

  xCoG /= totalCharge;
  yCoG /= totalCharge;
}

void EUTelAnnulusClusterImpl::getAnnulusCenterOfGravity(float& FxCoord, float& FyCoord, float &xCoG, float &yCoG) const
{
  xCoG = 0;
  yCoG = 0;

  double rCoG = 0;
  double phiCoG = 0;
  double totalCharge = 0;

  // Loop over all the pixels in the cluster
  auto &pixelVec = getPixels();
  for( auto &pixel : pixelVec ) {
    double curSignal = pixel.getSignal();
    rCoG += (pixel.getR()) * curSignal;
    phiCoG += (pixel.getAngle()) * curSignal;
    totalCharge += curSignal;
  }

  rCoG /= totalCharge;
  phiCoG /= totalCharge;
  xCoG = rCoG*cos(phiCoG) + FxCoord;
  yCoG = rCoG*sin(phiCoG) + FyCoord;
}

void EUTelAnnulusClusterImpl::getHitCovMatrix(double& res_phi, double& res_r, double (&cov)[4] ) const
{
  double xCoG = 0;
  double yCoG = 0;
  double totalCharge = 0;
  
  auto &pixelVec = getPixels();
  for( auto &pixel : pixelVec ) {
    double curSignal = pixel.getSignal();
    double theta_sin = sin(pixel.getAngle());
    double theta_cos = cos(pixel.getAngle());
    double radius = pixel.getR();
    
    xCoG += radius * theta_cos * curSignal;
    yCoG += radius * theta_sin * curSignal;
    totalCharge += curSignal;
  }
  xCoG /= totalCharge;
  yCoG /= totalCharge;

  for( auto &pixel : pixelVec ) {
    double curSignal = pixel.getSignal();
    double theta_sin = sin(pixel.getAngle());
    double theta_cos = cos(pixel.getAngle());
    double radius = pixel.getR();

    cov[0] += res_r*res_r*theta_cos*theta_cos*curSignal*curSignal + res_phi*res_phi*theta_sin*theta_sin*radius*radius*curSignal*curSignal;
    cov[1] += res_r*res_r*theta_sin*theta_sin*curSignal*curSignal + res_phi*res_phi*theta_cos*theta_cos*radius*radius*curSignal*curSignal;
    cov[2] += res_r*res_r*(xCoG*theta_cos+yCoG*theta_sin)*(xCoG*theta_cos+yCoG*theta_sin)*curSignal*curSignal 
      + res_phi*res_phi*(xCoG*theta_sin-yCoG*theta_cos)*(xCoG*theta_sin-yCoG*theta_cos)*radius*radius*curSignal*curSignal;
    cov[3] += res_r*res_r*(xCoG*theta_sin-yCoG*theta_cos)*(xCoG*theta_sin-yCoG*theta_cos)*curSignal*curSignal 
      + res_phi*res_phi*(xCoG*theta_cos+yCoG*theta_sin)*(xCoG*theta_cos+yCoG*theta_sin)*radius*radius*curSignal*curSignal ;
  }

  cov[0] = sqrt(cov[0]) / totalCharge; //resolution of x0
  cov[1] = sqrt(cov[1]) / totalCharge; //resolution of y0
  cov[2] = sqrt(cov[2])/sqrt(xCoG*xCoG + yCoG*yCoG)/totalCharge;  //resolution of r
  cov[3] = sqrt(cov[3])/(xCoG*xCoG + yCoG*yCoG)/totalCharge;  //resolution of phi
}
  
  
