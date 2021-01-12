/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

// personal includes ".h"
#include "EUTelAnnulusPixel.h"
#include "EUTELESCOPE.h"

// system includes <>
#include <iomanip>
#include <iostream>

using namespace eutelescope;

// Default constructor, returns pixel with all fields set to zero
EUTelAnnulusPixel::EUTelAnnulusPixel()
  : EUTelGenericSparsePixel(), _posX(0), _posY(0), _rmin(0), _rmax(0), _dphi(0), _angle(0), _rlength(0) {
  _noOfElementsDerived = 11;
  _typeDerived = kEUTelAnnulusPixel;
}

// Constructor taking all possible arguments
EUTelAnnulusPixel::EUTelAnnulusPixel(short xCoord, short yCoord, float signal, short time, 
				     float posX, float posY, float rmin, float rmax, 
				     float dphi, float rlength, float angle)
    : EUTelGenericSparsePixel(xCoord, yCoord, signal, time), 
      _posX(posX), _posY(posY), _rmin(rmin), _rmax(rmax), _dphi(dphi), _angle(angle), _rlength(rlength) {
  _noOfElementsDerived = 11;
  _typeDerived = kEUTelAnnulusPixel;
}

// Constructor taking a EUTelGenericSparsePixel, all geometry entries are set to zero
EUTelAnnulusPixel::EUTelAnnulusPixel(EUTelGenericSparsePixel const &genericPixel)
  : EUTelGenericSparsePixel(genericPixel), _posX(0), _posY(0), _rmin(0), _rmax(0), _dphi(0), _angle(0), _rlength(0) {
  _noOfElementsDerived = 11;
  _typeDerived = kEUTelAnnulusPixel;
}

// Constructor taking a EUTelGenericSparsePixel and all geometry related parameters
EUTelAnnulusPixel::EUTelAnnulusPixel(EUTelGenericSparsePixel const &genericPixel, 
				     float posX, float posY, float rmin, float rmax,
                                     float dphi, float rlength, float angle)
  : EUTelGenericSparsePixel(genericPixel), _posX(posX), _posY(posY), _rmin(rmin), _rmax(rmax), 
    _dphi(dphi), _angle(angle), _rlength(rlength) {
  _noOfElementsDerived = 11;
  _typeDerived = kEUTelAnnulusPixel;
}

unsigned int EUTelAnnulusPixel::getNoOfElements() const {
  return _noOfElementsDerived;
}

SparsePixelType EUTelAnnulusPixel::getSparsePixelType() const {
  return _typeDerived;
}

void EUTelAnnulusPixel::print(std::ostream &os) const {
  int bigWidth = 50;
  for(int i = 0; i < bigWidth; ++i) {
    os << "-";
  }
  os << std::endl;
  int width = 20;
  os << std::setw(width) << std::setiosflags(std::ios::left) << "Type: " << _type << std::endl
     << std::setw(width) << "Elements: " << _noOfElements << std::endl
     << std::setw(width) << "x coord: "  << _xCoord << std::endl
     << std::setw(width) << "y coord: "  << _yCoord << std::endl
     << std::setw(width) << "signal: "  << _signal << std::endl
     << std::setw(width) << "time: "  << _time << std::endl
     << std::setw(width) << "posX: "  << _posX << std::endl
     << std::setw(width) << "posY: "  << _posY << std::endl
     << std::setw(width) << "rmin: "  << _rmin << std::endl
     << std::setw(width) << "rmax: "  << _rmax << std::endl
     << std::setw(width) << "dphi: "  << _dphi << std::endl
     << std::setw(width) << "rlength: "  << _rlength << std::endl
     << std::setw(width) << "angle: "  << _angle << std::endl;
  for(int i = 0; i < bigWidth; ++i ) {
    os << "-";
  }
}
