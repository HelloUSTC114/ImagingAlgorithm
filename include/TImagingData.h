#ifndef TIMAGINGDATA_H
#define TIMAGINGDATA_H

#include "TObject.h"
#include "TVector2.h"
#include "TVector3.h"
#include <iostream>
#include "TMath.h"

namespace TImagingType
{

typedef TVector2 Point2D;
typedef TVector3 Point3D;
typedef TVector2 Vector2D;
typedef TVector3 Vector3D;
}; // namespace TImagingType
using namespace std;

class TPreImagingDataPart : public TObject
{
public:
    TPreImagingDataPart(double x1, double hx1, double y1, double hy1, double x2, double hx2, double y2, double hy2);
    TPreImagingDataPart(const TImagingType::Point2D &x1, const TImagingType::Point2D &x2, const TImagingType::Point2D &y1, const TImagingType::Point2D &y2);

    bool SetPoint(double x1, double hx1, double y1, double hy1, double x2, double hx2, double y2, double hy2);
    bool SetPoint(const TImagingType::Point2D &x1, const TImagingType::Point2D &x2, const TImagingType::Point2D &y1, const TImagingType::Point2D &y2);

public:
    TImagingType::Point2D fX1;
    TImagingType::Point2D fY1;
    TImagingType::Point2D fX2;
    TImagingType::Point2D fY2;

    ClassDef(TPreImagingDataPart, 1);
};

class TImagingData : public TObject
{
public:
    bool SetPoint(const TPreImagingDataPart &prei, const TPreImagingDataPart &preo);
    bool SetPoint(const TImagingType::Point2D &pxi1, const TImagingType::Point2D &pxi2, const TImagingType::Point2D &pxo1, const TImagingType::Point2D &pxo2, const TImagingType::Point2D &pyi1, const TImagingType::Point2D &pyi2, const TImagingType::Point2D &pyo1, const TImagingType::Point2D &pyo2);

    const TImagingType::Point3D &PI1() const { return fP3I1; }
    const TImagingType::Point3D &PI2() const { return fP3I2; }
    const TImagingType::Point3D &PO1() const { return fP3O1; }
    const TImagingType::Point3D &PO2() const { return fP3O2; }

    double GetScatteringAngle() const;
    ostream &Show(ostream &os) const;

    // Calculate corrected points
    static bool Combine3DPoints(const TImagingType::Point2D &px1, const TImagingType::Point2D &px2, const TImagingType::Point2D &py1, const TImagingType::Point2D &py2, TImagingType::Point3D &p1, TImagingType::Point3D &p2);
    static bool Combine3DPoints(const TPreImagingDataPart &pre, TImagingType::Point3D &p1, TImagingType::Point3D &p2);

    const TVector3 &I1() const { return fP3I1; }
    const TVector3 &I2() const { return fP3I2; }
    const TVector3 &O1() const { return fP3O1; }
    const TVector3 &O2() const { return fP3O2; }
    const TVector3 *p0() const { return &fP3I1; }
    const TVector3 *p1() const { return &fP3I2; }
    const TVector3 *q0() const { return &fP3O1; }
    const TVector3 *q1() const { return &fP3O2; }

private:
    TImagingType::Point2D fPXI1; // In 1 X
    TImagingType::Point2D fPXI2; // In 2 X
    TImagingType::Point2D fPXO1; // Out 1 X
    TImagingType::Point2D fPXO2; // Out 2 X
    TImagingType::Point2D fPYI1; // In 1 Y
    TImagingType::Point2D fPYI2; // In 2 Y
    TImagingType::Point2D fPYO1; // Out 1 Y
    TImagingType::Point2D fPYO2; // Out 2 Y

    TImagingType::Point3D fP3I1; //! Calculated 3D Points
    TImagingType::Point3D fP3I2; //!
    TImagingType::Point3D fP3O1; //!
    TImagingType::Point3D fP3O2; //!
    void Calculate3DPoint();
    bool CheckPointData();

    ClassDef(TImagingData, 1);
};

#ifdef __ROOTCLING__
#pragma link C++ class TImagingData;
#pragma link C++ class TPreImagingDataPart;
#endif

#endif
