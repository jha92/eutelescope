/*
 *   This source code is part of the Eutelescope package of Marlin.
 *   You are free to use this source files for your own development as
 *   long as it stays in a public research context. You are not
 *   allowed to use it for commercial purpose. You must put this
 *   header with author names in all development based on this file.
 *
 */

#ifndef EUTELANNULUSPIXEL_H
#define EUTELANNULUSPIXEL_H

// personal includes ".h"
#include "EUTELESCOPE.h"
#include "EUTelGenericSparsePixel.h"

// marlin includes ".h"

// lcio includes <.h>

// system includes <>

namespace eutelescope {

  //! Helper class decribing a annulus pixel

  class EUTelAnnulusPixel : public EUTelGenericSparsePixel {

  public:
    //! Default constructor with all arguments (individually)
    EUTelAnnulusPixel(short xCoord, short yCoord, float signal, short time,
		      float posX, float posY, float rmin, float rmax,
		      float dphi, float rlength, float angle);

    //! Default constructor with all arguments (EUTelGenericSparsePixel + geometry)
    EUTelAnnulusPixel(EUTelGenericSparsePixel const &genericPixel, 
		      float posX, float posY, float rmin, float rmax,
                      float dphi, float rlength, float angle);

    //! Default constructor with only EUTelescopeGenericSparsePixel (others are zero)
    EUTelAnnulusPixel(EUTelGenericSparsePixel const &genericPixel);

    //! Default constructor with no args (others are zero)
    EUTelAnnulusPixel();

    //! Destructor
    virtual ~EUTelAnnulusPixel() {}

    //! Get the number of elements in the data structure
    /*! This method returns the number of elements the sparse pixel
     *  contains.
     *
     *  @return The number of elements in the data structure
     */
    virtual unsigned int getNoOfElements() const;

    //! Get the sparse pixel type using the enumerator
    /*! This methods returns the sparse pixel type using the
     *  enumerator defined in EUTELESCOPE.h
     *
     *  Overloaded for derived class, since downcast should
     *  yield the sparse pixel type of the base class.
     *
     *  @return The sparse pixel type using the enumerator
     */
    virtual SparsePixelType getSparsePixelType() const;

    //! Print method
    /*! This method is used to print out the contents of the sparse
     *  pixel
     *
     *  @param os The input output stream
     */
    virtual void print(std::ostream &os) const;

    //! Setter functions for the geometry parameters
    void setRmin(float rmin) { _rmin = rmin; }
    void setRmax(float rmax) { _rmax = rmax; }
    void setDphi(float dphi) { _dphi = dphi; }
    void setAngle(float angle) { _angle = angle; }
    void setR(float rlength) { _rlength = rlength; }
    void setPosX(float posX) { _posX = posX; }
    void setPosY(float posY) { _posY = posY; }

    //! Getter functions for the geometry parameters
    inline float getRmin() const { return _rmin; }
    inline float getRmax() const { return _rmax; }
    inline float getDphi() const { return _dphi; }
    inline float getAngle() const { return _angle; }
    inline float getR() const { return _rlength; }
    inline float getPosX() const { return _posX; }
    inline float getPosY() const { return _posY; }

  protected:
    //! X position on the sensor plane
    float _posX;

    //! Y position on the sensor plane
    float _posY;

    //! Minimal radius of the sensor plane
    float _rmin;

    //! Maximal radius of the sensor plane
    float _rmax;

    //! Phi difference of the sensor plane
    float _dphi;

    //! Angle of the sensor plane
    float _angle;

    //! Radial length of the strips
    float _rlength;

    //! Number of elements in the data structure
    unsigned int _noOfElementsDerived;

    //! The sparse pixel type enumerator for the derived type
    /*! Required since a downcast to EUTelGenericSparsePixel
     *  needs to return the downcast type and member variable
     *  overloading is not possible in C++.
     */
    SparsePixelType _typeDerived;
  };
} // namespace eutelescope

#endif
