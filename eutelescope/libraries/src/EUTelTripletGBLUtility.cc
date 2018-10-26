// -*- mode: c++; mode: auto-fill; mode: flyspell-prog; -*-

/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#if defined USE_GEAR

// EUTelescope includes:
#include "EUTelTripletGBLUtility.h"
#include "EUTelGeometryTelescopeGeoDescription.h"
//#include "EUTelTripletGBLDUTscatInstance.h"

#include "EUTELESCOPE.h"
#include <cmath>
#include <memory>
#include <type_traits>

// AIDA includes <.h>
#if defined(USE_AIDA) || defined(MARLIN_USE_AIDA)
#include <AIDA/IBaseHistogram.h>
#include <marlin/AIDAProcessor.h>
#include <AIDA/IHistogramFactory.h>
#include <AIDA/IAxis.h>
#include <AIDA/IHistogram1D.h>
#include <AIDA/IHistogram2D.h>
#include <AIDA/IProfile1D.h>
#include <AIDA/IProfile2D.h>
#include <AIDA/ITree.h>
#endif

// marlin includes ".h"
#include "marlin/Processor.h"
#include "marlin/Exceptions.h"
#include "marlin/ProcessorMgr.h"
#include "marlin/Global.h"

using namespace std;
using namespace eutelescope;
using namespace marlin;


EUTelTripletGBLUtility::EUTelTripletGBLUtility(){}

Eigen::Matrix<double, 5,5> EUTelTripletGBLUtility::JacobianPointToPoint( double ds ) {
  /* for GBL:
     Jacobian for straight line track
     track = q/p, x', y', x, y
     0,   1,  2,  3, 4
     */
  Eigen::Matrix<double, 5,5> jac = Eigen::Matrix<double, 5,5>::Identity();
  jac(3,1) = ds; // x = x0 + xp * ds
  jac(4,2) = ds; // y = y0 + yp * ds
  return jac;
}

void EUTelTripletGBLUtility::bookHistos(){
    
  marlin::AIDAProcessor::tree(parent)->mkdir("GBLUtility");

  sixkxHisto = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "GBLUtility/sixkx", 100, -10, 10 );
  sixkxHisto->setTitle( "kink x;kink x [mrad];triplet-driplet pairs" );

  sixkyHisto = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "GBLUtility/sixky", 100, -10, 10 );
  sixkyHisto->setTitle( "kink y;kink y [mrad];triplet-driplet pairs" );

  sixdxHisto = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "GBLUtility/sixdx", 100, -1, 1 );
  sixdxHisto->setTitle( "six match x;match x [mm];triplet-driplet pairs" );

  sixdyHisto = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "GBLUtility/sixdy", 100, -1, 1 );
  sixdyHisto->setTitle( "six match y;match y [mm];triplet-driplet pairs" );

  sixdxcHisto = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "GBLUtility/sixdxc", 100, -250, 250 );
  sixdxcHisto->setTitle( "six match x;track #Deltax[#mum];triplet-driplet pairs" );

  sixdycHisto = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "GBLUtility/sixdyc", 100, -250, 250 );
  sixdycHisto->setTitle( "six match y;track #Deltay[#mum];triplet-driplet pairs" );

  sixkxcHisto = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "GBLUtility/sixkxc", 100, -10, 10 );
  sixkxcHisto->setTitle( "kink x, x-y matched;kink x [mrad];triplet-driplet pairs" );

  sixkycHisto = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "GBLUtility/sixkyc", 100, -10, 10 );
  sixkycHisto->setTitle( "kink y, x-y matched;kink y [mrad];triplet-driplet pairs" );

  sixxHisto = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "GBLUtility/sixx", 240, -12, 12 );
  sixxHisto->setTitle( "six x at DUT;six x_{out} at DUT [mm];six-plane tracks" );

  sixyHisto = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "GBLUtility/sixy", 120, -6, 6 );
  sixyHisto->setTitle( "six y at DUT;six y_{up} at DUT [mm];six-plane tracks" );

  sixxyHisto = AIDAProcessor::histogramFactory(parent)->createHistogram2D( "GBLUtility/sixxy", 240, -12, 12, 120, -6, 6 );
  sixxyHisto->setTitle( "six at DUT;six x_{out} at DUT [mm];six y_{up} at DUT [mm];six-plane tracks" );

  sixxycHisto = AIDAProcessor::histogramFactory(parent)->createHistogram2D( "GBLUtility/sixxyc", 240, -12, 12, 120, -6, 6 );
  sixxycHisto->setTitle( "six large kink;six x_{out} at DUT [mm];six y_{up} at DUT [mm];large kink tracks" );

  kinkxvsxy = AIDAProcessor::histogramFactory(parent)->createProfile2D( "GBLUtility/kinkxvsxy", 120, -12, 12, 60, -6, 6, 0, 100 );
  kinkxvsxy->setTitle( "kink x;six x_{out} at DUT [mm];six y_{up} at DUT [mm];sqrt(<kink^{2}>) [mrad]" );

  kinkyvsxy = AIDAProcessor::histogramFactory(parent)->createProfile2D( "GBLUtility/kinkyvsxy", 120, -12, 12, 60, -6, 6, 0, 100 );
  kinkyvsxy->setTitle( "kink y;six x_{out} at DUT [mm];six y_{up} at DUT [mm];sqrt(<kink^{2}>) [mrad]" );

  kinkxyvsxy = AIDAProcessor::histogramFactory(parent)->createProfile2D( "GBLUtility/kinkxyvsxy", 120, -12, 12, 60, -6, 6, 0, 100 );
  kinkxyvsxy->setTitle( "kink <|x|+|y|>;six x_{out} at DUT [mm];six y_{up} at DUT [mm]; (|x| + |y|)/2 [mrad]" );

  kinkx = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "GBLUtility/kinkx", 500, -5, 5 );
  kinkx->setTitle( "triplet kink x angle at DUT;x angle at DUT [mrad];tracks" );

  kinky = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "GBLUtility/kinky", 500, -5, 5 );
  kinky->setTitle( "triplet kink y angle at DUT;y angle at DUT [mrad];tracks" );

  kinkxy = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "GBLUtility/kinkxy", 500, -5, 5 );
  kinkxy->setTitle( "triplet kink xy angle at DUT;xy angle at DUT [mrad];tracks" );
  
  //cut plots
  marlin::AIDAProcessor::tree(parent)->mkdir("Cuts");
  
  upstreamTripletSlopeX = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/upstreamTripletSlopeCutX", 1000, -3, 3 ); //binning to be reviewed
  upstreamTripletSlopeX->setTitle( "Upstream Triplet Slope X;Upstream Triplet Slope X [mrad];Counts" );
  
  upstreamTripletSlopeY = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/upstreamTripletSlopeCutY", 1000, -3, 3 ); //binning to be reviewed
  upstreamTripletSlopeY->setTitle( "Upstream Triplet Slope Y;Upstream Triplet Slope Y [mrad];Counts" );
  
  downstreamTripletSlopeX = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/downstreamTripletSlopeCutX", 1000, -3, 3 ); //binning to be reviewed
  downstreamTripletSlopeX->setTitle( "Downstream Triplet Slope X;Downstream Triplet Slope X [mrad];Counts" );
  
  downstreamTripletSlopeY = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/downstreamTripletSlopeCutY", 1000, -3, 3 ); //binning to be reviewed
  downstreamTripletSlopeY->setTitle( "Downstream Triplet Slope Y;Downstream Triplet Slope Y [mrad];Counts" );
  
  upstreamTripletResidualX = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/upstreamTripletResidualCutX", 1000, -3, 3 ); //binning to be reviewed
  upstreamTripletResidualX->setTitle( "Upstream Triplet Residual X;Upstream Triplet Residual X [mm];Counts" );
  
  upstreamTripletResidualY = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/upstreamTripletResidualCutY", 1000, -3, 3 ); //binning to be reviewed
  upstreamTripletResidualY->setTitle( "Upstream Triplet Residual Y;Upstream Triplet Residual Y [mm];Counts" );
  
  downstreamTripletResidualX = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/downstreamTripletResidualCutX", 1000, -3, 3 ); //binning to be reviewed
  downstreamTripletResidualX->setTitle( "Downstream Triplet Residual X;Downstream Triplet Residual X [mm];Counts" );
  
  downstreamTripletResidualY = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/downstreamTripletResidualCutY", 1000, -3, 3 ); //binning to be reviewed
  downstreamTripletResidualY->setTitle( "Downstream Triplet Residual Y;Downstream Triplet Residual Y [mm];Counts" );
  
  tripletMatchingResidualX = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/tripletMatchingResidualCutX", 1000, -3, 3 ); //binning to be reviewed
  tripletMatchingResidualX->setTitle( "Triplet Matching Residual X;Triplet Matching Residual X [mm];Counts" );
  
  tripletMatchingResidualY = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/tripletMatchingResidualCutY", 1000, -3, 3 ); //binning to be reviewed
  tripletMatchingResidualY->setTitle( "Triplet Matching Residual Y;Triplet Matching Residual Y [mm];Counts" );
  
  DUTMatchingResidualX = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/DUTMatchingResidualCutX", 1000, -3, 3 ); //binning to be reviewed
  DUTMatchingResidualX->setTitle( "DUT Matching Residual local X;DUT Matching Residual local X [mm];Counts" );
  
  DUTMatchingResidualY = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/DUTMatchingResidualCutY", 1000, -3, 3 ); //binning to be reviewed
  DUTMatchingResidualY->setTitle( "DUT Matching Residual local Y;DUT Matching Residual local Y [mm];Counts" );
  
  DUTHitNumber = AIDAProcessor::histogramFactory(parent)->createHistogram1D( "Cuts/DUTHitNumber", 21, -0.5, 20.5 ); //binning to be reviewed
  DUTHitNumber->setTitle( "Number of Hits matched to a track per DUT ID;DUT ID;Number of Hits matched to a track" );

}

void EUTelTripletGBLUtility::MatchTriplets(std::vector<triplet> const & up, std::vector<EUTelTripletGBLUtility::triplet> const & down, double z_match, double trip_matching_cut, std::vector<EUTelTripletGBLUtility::track> &tracks) {

  // Cut on the matching of two triplets [mm]

  for( auto& trip: up ){

    // Track impact position at Matching Point from Upstream:
    double xA = trip.getx_at(z_match); // triplet impact point at matching position
    double yA = trip.gety_at(z_match);

    // check if trip is isolated. use at least double the trip_machting_cut for isolation in order to avoid double matching
    bool IsolatedTrip = IsTripletIsolated(trip, up, z_match, trip_matching_cut*2.0001);
    streamlog_out(DEBUG4) << "  Is triplet isolated? " << IsolatedTrip << std::endl;

    for( auto& drip: down ){

      // Track impact position at Matching Point from Downstream:
      double xB = drip.getx_at(z_match); // triplet impact point at matching position
      double yB = drip.gety_at(z_match);

      // check if drip is isolated
      bool IsolatedDrip = IsTripletIsolated(drip, down, z_match, trip_matching_cut*2.0001);
      streamlog_out(DEBUG4) << "  Is driplet isolated? " << IsolatedDrip << std::endl;


      // Build a track candidate from one upstream and one downstream triplet:
      EUTelTripletGBLUtility::track newtrack(trip,drip);

      // Track kinks as difference in triplet slopes:
      double kx = newtrack.kink_x();
      double ky = newtrack.kink_y();

      // driplet - triplet
      double dx = xB - xA; 
      double dy = yB - yA;

      sixkxHisto->fill( kx*1E3 );
      sixkyHisto->fill( ky*1E3 );
      sixdxHisto->fill( dx );
      sixdyHisto->fill( dy );

      if( abs(dy) < 0.5 ) sixdxcHisto->fill( dx*1E3 );
      if( abs(dx) < 0.5 ) sixdycHisto->fill( dy*1E3 );
      
      //cut plots
      tripletMatchingResidualX->fill(dx);
      tripletMatchingResidualY->fill(dy);
      // match driplet and triplet:
      streamlog_out(DEBUG4) << "  Distance for matching x: " << fabs(dx)<< std::endl;
      streamlog_out(DEBUG4) << "  Distance for matching y: " << fabs(dy)<< std::endl;
      if( fabs(dx) > trip_matching_cut) continue;
      if( fabs(dy) > trip_matching_cut) continue;
      streamlog_out(DEBUG4) << "  Survived matching " << std::endl;

      // check isolation
      if( !IsolatedTrip || !IsolatedDrip ) {
	continue;
      }
      streamlog_out(DEBUG4) << "  Trip and Drip isolated " << std::endl;      

      sixkxcHisto->fill( kx*1E3 );
      sixkycHisto->fill( ky*1E3 );
      sixxHisto->fill( -xA ); // -xA = x_DP = out
      sixyHisto->fill( -yA ); // -yA = y_DP = up
      sixxyHisto->fill( -xA, -yA ); // DP: x_out, y_up
      // Fill kink map histogram:
      if( abs( kx ) > 0.002 || abs( ky ) > 0.002 ) sixxycHisto->fill( -xA, -yA );      

      kinkx->fill( kx*1E3 ); //sqrt(<kink^2>) [mrad]
      kinky->fill( ky*1E3 ); //sqrt(<kink^2>) [mrad]
      kinkxy->fill( (fabs(kx)+fabs(ky))/2*1E3 ); // [mrad]
      kinkxvsxy->fill( -xA, -yA, fabs(kx)*1E3 ); //sqrt(<kink^2>) [mrad]
      kinkyvsxy->fill( -xA, -yA, fabs(ky)*1E3 ); //sqrt(<kink^2>) [mrad]
      kinkxyvsxy->fill( -xA, -yA, (fabs(kx) + fabs(ky))/2*1E3 ); // [mrad]

      // Add the track to the vector if trip/drip are isolated, the triplets are matched, and all other cuts are passed
      tracks.push_back(newtrack);

    } // Downstream
  } // Upstream

  streamlog_out(DEBUG2) << "Found " << tracks.size() << " tracks from matched triplets." << std::endl;
  //return tracks;
}

bool EUTelTripletGBLUtility::IsTripletIsolated(EUTelTripletGBLUtility::triplet const & it, std::vector<EUTelTripletGBLUtility::triplet> const & trip, double z_match, double isolation_cut) { // isolation_cut is defaulted to 0.3 mm
  bool IsolatedTrip = true;

  // check first if trip is isolated
  double xA = it.getx_at(z_match); // triplet impact point at matching position
  double yA = it.gety_at(z_match);

  double ddAMin = -1.0;
  for( auto& tripIsoCheck: trip ) {
    if( &it != &tripIsoCheck ){
      double xAIsoCheck = tripIsoCheck.getx_at(z_match);
      double yAIsoCheck = tripIsoCheck.gety_at(z_match);
      double ddA = sqrt( fabs(xAIsoCheck - xA)*fabs(xAIsoCheck - xA) 
	  + fabs(yAIsoCheck - yA)*fabs(yAIsoCheck - yA) );
      if(ddAMin < 0 || ddA < ddAMin) ddAMin = ddA;
	}
  }

  if(ddAMin < isolation_cut && ddAMin > -0.5) IsolatedTrip = false; // if there is only one triplet, ddAmin is still -1.

  return IsolatedTrip;
}

bool EUTelTripletGBLUtility::AttachDUT(EUTelTripletGBLUtility::triplet & triplet, std::vector<EUTelTripletGBLUtility::hit> const & hits, unsigned int dutID,  std::vector<float> dist_cuts){

	auto zPos = geo::gGeometry().getPlaneZPosition(dutID);
	int minHitIx = -1;
	double minDist = std::numeric_limits<float>::max();

	auto trX = triplet.getx_at(zPos);
	auto trY = triplet.gety_at(zPos);
	double trglpos[3] = {trX,trY,zPos};
	double trlocpos[3];
    geo::gGeometry().master2Local(dutID,trglpos,trlocpos);
    
	size_t ix = 0;	
	for(auto& hit: hits) {
		if(hit.plane == dutID) {
			auto hitX = hit.x;
			auto hitY = hit.y;
			//track and DUT hit positions are transformed in the local coordinate frame of the DUT, where the cuts are applied
			//it is horrible to convert the hits from global to local, where the latter is actually basically the input of the hitmaker
			//this is a first step toward a more general refactoring to work mainly in local coordinate systems 
			double hitglpos[3] = {hitX,hitY,zPos};
	        double hitlocpos[3];
            geo::gGeometry().master2Local(dutID,hitglpos,hitlocpos);
			auto distX = fabs(trlocpos[0]-hitlocpos[0]);
			auto distY = fabs(trlocpos[1]-hitlocpos[1]);
                        double dist = distX*distX + distY*distY;
			DUTMatchingResidualX->fill(trlocpos[0]-hitlocpos[0]);
			DUTMatchingResidualY->fill(trlocpos[1]-hitlocpos[1]);
			if(distX <= dist_cuts.at(0) && distY <= dist_cuts.at(1) && dist < minDist ){
				minHitIx = static_cast<int>(ix);
				DUTHitNumber->fill(dutID);
				minDist = dist;
			}
		}
		++ix;
	}

	if(minHitIx != -1) {
		triplet.push_back_DUT(hits[minHitIx].plane, hits[minHitIx]);
		return true;
	}
return false;
}

double EUTelTripletGBLUtility::PlaneEfficiency(std::vector<EUTelTripletGBLUtility::triplet> &eff_triplets_UP, std::vector<EUTelTripletGBLUtility::triplet> &eff_triplets_DOWN, std::vector<EUTelTripletGBLUtility::hit> &hits, unsigned int PUT, double track_match_z, double DUTz, double track_match_cut, double eff_radius, std::vector<AIDA::IProfile1D*> &profile) {

  std::vector<EUTelTripletGBLUtility::triplet> eff_triplets;

  for( auto& trip: eff_triplets_UP ) {

    // Track impact position at Matching Point from Upstream:
    double xA = trip.getx_at(track_match_z); // triplet impact point at matching position
    double yA = trip.gety_at(track_match_z);

    // check if trip is isolated
    bool IsolatedTrip = IsTripletIsolated(trip, eff_triplets_UP, track_match_z, track_match_cut*2.0001);

    for( auto& drip: eff_triplets_DOWN ){

      // Track impact position at Matching Point from Downstream:
      double xB = drip.getx_at(track_match_z); // triplet impact point at matching position
      double yB = drip.gety_at(track_match_z);

      // check if drip is isolated
      bool IsolatedDrip = IsTripletIsolated(drip, eff_triplets_DOWN, track_match_z, track_match_cut*2.0001);

      // driplet - triplet
      double dx = xB - xA; 
      double dy = yB - yA;

      // match driplet and triplet:
      if( fabs(dx) > track_match_cut) continue;
      if( fabs(dy) > track_match_cut) continue;

      //std::cout << " intersec ";

      // check isolation
      if( !IsolatedTrip || !IsolatedDrip ) continue;
      //std::cout << " , isolated ";

      eff_triplets.emplace_back(trip);


    } // Downstream
  } // Upstream

  //std::cout << " n eff triplets = " << eff_triplets.size() << std::endl;

  int n_sum = 0;
  int n_matched = 0;
  for( auto& trip: eff_triplets ) {

    double ddAMin = -1.0;
    // extrapolate triplet to plane  under test
    double xTrip = trip.getx_at(DUTz);
    double yTrip = trip.gety_at(DUTz);

    for( auto& lhit: hits ){

      // Fill residuals of triplet and hit in the selected plane:
      if( lhit.plane == PUT ) {
	double xHit = lhit.x;
	double yHit = lhit.y;

	double ddA = sqrt( fabs(xHit - xTrip)*fabs(xHit - xTrip) 
	    + fabs(yHit - yTrip)*fabs(yHit - yTrip) );
	if(ddAMin < 0 || ddA < ddAMin) ddAMin = ddA;
      }
    } // end loop over hits

    // if distance is smaller then limit, accept this as matched Hit
    if(fabs(ddAMin) < eff_radius) {
      //n_matched_trips++;
      profile.at(0)->fill(-xTrip, 1.);
      profile.at(1)->fill(-yTrip, 1.);
      n_matched++;
      n_sum++;
    } else {
      profile.at(0)->fill(-xTrip, 0.);
      profile.at(1)->fill(-yTrip, 0.);
      n_sum++;
    }

  }
  
  eff_triplets_UP.clear();
  eff_triplets_DOWN.clear();

  double eff = static_cast<double>(n_matched)/static_cast<double>(n_sum);
  return eff;

}

EUTelTripletGBLUtility::track::track(triplet up, triplet down) : upstream(up), downstream(down) {}

double EUTelTripletGBLUtility::track::kink_x() {
  return (downstream.slope().x - upstream.slope().x);
}

double EUTelTripletGBLUtility::track::kink_y() {
  return (downstream.slope().y - upstream.slope().y);
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::track::intersect() {
  hit inter;
  // Re-check what this actually is...
  // and simplify using triplet class members...
  inter.x = ( upstream.base().x - upstream.slope().x * upstream.base().z - downstream.base().x + downstream.slope().x * downstream.base().z ) / kink_x();
  inter.y = ( upstream.base().y - upstream.slope().y * upstream.base().z - downstream.base().y + downstream.slope().y * downstream.base().z ) / kink_y();
  return inter;
}

EUTelTripletGBLUtility::triplet& EUTelTripletGBLUtility::track::get_upstream() {
  return upstream;
}

EUTelTripletGBLUtility::triplet& EUTelTripletGBLUtility::track::get_downstream() {
  return downstream;
}

EUTelTripletGBLUtility::hit const & EUTelTripletGBLUtility::track::gethit(int plane) {
  if(plane < 3) return upstream.gethit(plane);
  else return downstream.gethit(plane);
}



EUTelTripletGBLUtility::triplet::triplet() : linked_dut(false), hits() {
  // Empty default constructor
}

EUTelTripletGBLUtility::triplet::triplet(hit hit0, hit hit1, hit hit2) : linked_dut(false), hits() {
  triplet();
  filltriplet(hit0, hit1, hit2);
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::triplet::getpoint_at(double z) const{
  hit impact;
  impact.z = z - base().z;
  impact.x = base().x + slope().x * impact.z;
  impact.y = base().y + slope().y * impact.z;
  return impact;
}

double EUTelTripletGBLUtility::triplet::getx_at(double z) const {
  return base().x + slope().x * (z - base().z);
}

double EUTelTripletGBLUtility::triplet::getdx() const {
  return hits.rbegin()->second.x - hits.begin()->second.x;
}

double EUTelTripletGBLUtility::triplet::getdx(int ipl) const {
  return hits.at(ipl).x - base().x - slope().x * (hits.at(ipl).z - base().z);
}

double EUTelTripletGBLUtility::triplet::getdx(hit point) const {
  return point.x - base().x - slope().x * (point.z - base().z);
}

double EUTelTripletGBLUtility::triplet::gety_at(double z) const {
  return base().y + slope().y * (z - base().z);
}

double EUTelTripletGBLUtility::triplet::getdy() const {
  return hits.rbegin()->second.y - hits.begin()->second.y;
}

double EUTelTripletGBLUtility::triplet::getdy(int ipl) const {
  return hits.at(ipl).y - base().y - slope().y * (hits.at(ipl).z - base().z);
}

double EUTelTripletGBLUtility::triplet::getdy(hit point) const {
  return point.y - base().y - slope().y * (point.z - base().z);
}

double EUTelTripletGBLUtility::triplet::getdz() const {
  return hits.rbegin()->second.z - hits.begin()->second.z;
}

EUTelTripletGBLUtility::hit const & EUTelTripletGBLUtility::triplet::gethit(int plane) const {
  return hits.at(plane);
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::triplet::base() const {
  hit center;
  center.x = 0.5*( hits.begin()->second.x + hits.rbegin()->second.x );
  center.y = 0.5*( hits.begin()->second.y + hits.rbegin()->second.y );
  center.z = 0.5*( hits.begin()->second.z + hits.rbegin()->second.z );
  return center;
}

EUTelTripletGBLUtility::hit EUTelTripletGBLUtility::triplet::slope() const {
  hit sl;
  double dz = (hits.rbegin()->second.z - hits.begin()->second.z);
  sl.x = (hits.rbegin()->second.x - hits.begin()->second.x) / dz;
  sl.y = (hits.rbegin()->second.y - hits.begin()->second.y) / dz;
  return sl;
}

std::pair<double,double> EUTelTripletGBLUtility::doIterativeGaussianFit(AIDA::IHistogram1D* in_hist, int need_rebin) const {
    //--- First convert from IHistogram to TH1 so that we could use ROOT's fitters
    int nb_bins = in_hist->axis().bins(); //, nb_entries = in_hist->allEntries();
    streamlog_out (MESSAGE4) << in_hist->title() << " has " << in_hist->allEntries() << " entries" << std::endl;
    
    TH1D current(in_hist->title().data(), in_hist->title().data(), nb_bins, in_hist->axis().lowerEdge(), in_hist->axis().upperEdge());
    for(int id = 0 ; id < nb_bins ; id++) current.SetBinContent(id+1, in_hist->binEntries(id));
    current.Rebin(need_rebin);
    
    double startForMaxFraction = 0.2;
    double max = current.GetMaximum();
    int min_bin = current.FindFirstBinAbove(startForMaxFraction*max), max_bin = current.FindLastBinAbove(startForMaxFraction*max);
    if(min_bin == 1) {
        startForMaxFraction = 0.35;
        min_bin = current.FindFirstBinAbove(startForMaxFraction*max);
        max_bin = current.FindLastBinAbove(startForMaxFraction*max);
    }
    if(min_bin == 1) {
        startForMaxFraction = 0.5;
        min_bin = current.FindFirstBinAbove(startForMaxFraction*max);
        max_bin = current.FindLastBinAbove(startForMaxFraction*max);
    }
    for(int id = min_bin ; id < max_bin-1 ; id++) {
        if( (current.GetBinContent(id) > 0.7*current.GetBinContent(id+1))&&(current.GetBinContent(id+2) > 0.7*current.GetBinContent(id+1)) ) {
            current.Smooth(2);
            min_bin = current.FindFirstBinAbove(startForMaxFraction*max);
            max_bin = current.FindLastBinAbove(startForMaxFraction*max);
        }
    }
    
    double bound_low = current.GetBinCenter( min_bin );
    double bound_hi = current.GetBinCenter( max_bin );
    double mean = (bound_hi-bound_low)/2.;
    double sigma = bound_hi-bound_low;
    //int nb_iter = 0;
    
    TF1 gausfit("gausfit", "([3]+[0]*exp(-0.5*( ((x-[1])/[2])*((x-[1])/[2]) ) ))", bound_low, bound_hi);
    streamlog_out (MESSAGE4) << "First fit between " << bound_low << " " << bound_hi << std::endl;
    
    gausfit.SetParameter(0, max); //-current.GetBinContent(1));
    gausfit.SetParameter(1, mean);
    gausfit.SetParameter(2, sigma);
    gausfit.SetParameter(3, current.GetBinContent(1));
    
    gausfit.SetParLimits(0, 0.1*max, 2*max);
    gausfit.SetParLimits(1, bound_low, bound_hi);
    gausfit.SetParLimits(2, 0.1*(bound_hi-bound_low), 1.5*(bound_hi-bound_low));
    gausfit.SetParLimits(3, 0, 0.2*max);
             
    TFitResultPtr fitresult = current.Fit(&gausfit,"S","",bound_low,bound_hi);
        
    mean = fitresult->GetParams()[1];
    sigma = fitresult->GetParams()[2];
    
    TCanvas can; can.cd();
    current.Draw();
    can.SaveAs(TString(in_hist->title())+".pdf");
    
    return std::pair<double,double> (mean, sigma);
}

void EUTelTripletGBLUtility::determineBestCuts() const {
    
    std::pair<double, double> res_upstreamTripletResidualX   = doIterativeGaussianFit(upstreamTripletResidualX, 5);
    std::pair<double, double> res_upstreamTripletSlopeX      = doIterativeGaussianFit(upstreamTripletSlopeX, 10);
    std::pair<double, double> res_upstreamTripletSlopeY      = doIterativeGaussianFit(upstreamTripletSlopeY, 10);
    std::pair<double, double> res_downstreamTripletSlopeX    = doIterativeGaussianFit(downstreamTripletSlopeX, 10);
    std::pair<double, double> res_downstreamTripletSlopeY    = doIterativeGaussianFit(downstreamTripletSlopeY, 10);
    std::pair<double, double> res_upstreamTripletResidualY   = doIterativeGaussianFit(upstreamTripletResidualY, 5);
    std::pair<double, double> res_downstreamTripletResidualX = doIterativeGaussianFit(downstreamTripletResidualX, 5);
    std::pair<double, double> res_downstreamTripletResidualY = doIterativeGaussianFit(downstreamTripletResidualY, 5);
    std::pair<double, double> res_tripletMatchingResidualX   = doIterativeGaussianFit(tripletMatchingResidualX, 5);
    std::pair<double, double> res_tripletMatchingResidualY   = doIterativeGaussianFit(tripletMatchingResidualY, 5);
    std::pair<double, double> res_DUTMatchingResidualX       = doIterativeGaussianFit(DUTMatchingResidualX, 5);
    std::pair<double, double> res_DUTMatchingResidualY       = doIterativeGaussianFit(DUTMatchingResidualY, 5);
    
    for(int id = 0 ; id < 10 ; id++) streamlog_out (MESSAGE4) << std::endl;
    streamlog_out (MESSAGE4) << "__________ *** An iterative gaussian fit determined the following parameters for the distributions we will cut on : __________" << std::endl;
    streamlog_out (MESSAGE4) << "upstreamTripletResidualX --- Mean = " << res_upstreamTripletResidualX.first << " ; Sigma = " << res_upstreamTripletResidualX.second << std::endl;
    streamlog_out (MESSAGE4) << "upstreamTripletSlopeX --- Mean = " << res_upstreamTripletSlopeX.first << " ; Sigma = " << res_upstreamTripletSlopeX.second << std::endl;
    streamlog_out (MESSAGE4) << "upstreamTripletSlopeY --- Mean = " << res_upstreamTripletSlopeY.first << " ; Sigma = " << res_upstreamTripletSlopeY.second << std::endl;
    streamlog_out (MESSAGE4) << "downstreamTripletSlopeX --- Mean = " << res_downstreamTripletSlopeX.first << " ; Sigma = " << res_downstreamTripletSlopeX.second << std::endl;
    streamlog_out (MESSAGE4) << "downstreamTripletSlopeY --- Mean = " << res_downstreamTripletSlopeY.first << " ; Sigma = " << res_downstreamTripletSlopeY.second << std::endl;
    streamlog_out (MESSAGE4) << "upstreamTripletResidualY --- Mean = " << res_upstreamTripletResidualY.first << " ; Sigma = " << res_upstreamTripletResidualY.second << std::endl;
    streamlog_out (MESSAGE4) << "downstreamTripletResidualX --- Mean = " << res_downstreamTripletResidualX.first << " ; Sigma = " << res_downstreamTripletResidualX.second << std::endl;
    streamlog_out (MESSAGE4) << "downstreamTripletResidualY --- Mean = " << res_downstreamTripletResidualY.first << " ; Sigma = " << res_downstreamTripletResidualY.second << std::endl;
    streamlog_out (MESSAGE4) << "tripletMatchingResidualX --- Mean = " << res_tripletMatchingResidualX.first << " ; Sigma = " << res_tripletMatchingResidualX.second << std::endl;
    streamlog_out (MESSAGE4) << "tripletMatchingResidualY --- Mean = " << res_tripletMatchingResidualY.first << " ; Sigma = " << res_tripletMatchingResidualY.second << std::endl;
    streamlog_out (MESSAGE4) << "DUTMatchingResidualX --- Mean = " << res_DUTMatchingResidualX.first << " ; Sigma = " << res_DUTMatchingResidualX.second << std::endl;
    streamlog_out (MESSAGE4) << "DUTMatchingResidualY --- Mean = " << res_DUTMatchingResidualY.first << " ; Sigma = " << res_DUTMatchingResidualY.second << std::endl;
    for(int id = 0 ; id < 10 ; id++) streamlog_out (MESSAGE4) << std::endl;
    streamlog_out (MESSAGE4) << "__________ *** We then recommend the following cuts, although they should be checked ! : __________" << std::endl;
    streamlog_out (MESSAGE4) << "UpstreamTripletCut 	= " << std::max( fabs(res_upstreamTripletResidualX.first)+(3.)*res_upstreamTripletResidualX.second, fabs(res_upstreamTripletResidualY.first)+(3.)*res_upstreamTripletResidualY.second) << std::endl;
    streamlog_out (MESSAGE4) << "DownstreamTripletCut 	= " << std::max( fabs(res_downstreamTripletResidualX.first)+(3.)*res_downstreamTripletResidualX.second, fabs(res_downstreamTripletResidualY.first)+(3.)*res_downstreamTripletResidualY.second) << std::endl;
    streamlog_out (MESSAGE4) << "UpstreamSlopeCut   	= " << std::max( fabs(res_upstreamTripletSlopeX.first)+(3.)*res_upstreamTripletSlopeX.second, fabs(res_upstreamTripletSlopeY.first)+(3.)*res_upstreamTripletSlopeY.second) << std::endl;
    streamlog_out (MESSAGE4) << "DownstreamSlopeCut 	= " << std::max( fabs(res_downstreamTripletSlopeX.first)+(3.)*res_downstreamTripletSlopeX.second, fabs(res_downstreamTripletSlopeY.first)+(3.)*res_downstreamTripletSlopeY.second) << std::endl;
    streamlog_out (MESSAGE4) << "TripletMatchingCut 	= " << std::max( fabs(res_tripletMatchingResidualX.first)+(3.)*res_tripletMatchingResidualX.second, fabs(res_tripletMatchingResidualY.first)+(3.)*res_tripletMatchingResidualY.second) << std::endl;
    streamlog_out (MESSAGE4) << "DUTCuts 		= " << fabs(res_DUTMatchingResidualX.first)+(3.)*res_DUTMatchingResidualX.second << " " <<  fabs(res_DUTMatchingResidualY.first)+(3.)*res_DUTMatchingResidualY.second << std::endl;
    for(int id = 0 ; id < 10 ; id++) streamlog_out (MESSAGE4) << std::endl;
}

#endif
