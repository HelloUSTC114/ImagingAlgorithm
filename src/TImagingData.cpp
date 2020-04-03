#include "TImagingData.h"
#include "PoCA.h"

using namespace TImagingType;

ClassImp(TImagingData);
ClassImp(TPreImagingDataPart);

bool TImagingData::SetPoint(const TPreImagingDataPart &prei, const TPreImagingDataPart &preo)
{
     SetPoint(prei.fX1, prei.fX2, preo.fX1, preo.fX2, prei.fY1, prei.fY2, preo.fY1, preo.fY2);
}

bool TImagingData::SetPoint(const Point2D &pxi1, const Point2D &pxi2, const Point2D &pxo1, const Point2D &pxo2, const Point2D &pyi1, const Point2D &pyi2, const Point2D &pyo1, const Point2D &pyo2)
{
     double z1 = pxi1.Y();
     double z2 = pxi2.Y();
     double z3 = pxo1.Y();
     double z4 = pxo2.Y();

     // Rearrange all detectors
     if (z1 > z4)
     {
          if (z1 > z2)
          {
               fPXI1 = pxi1;
               fPXI2 = pxi2;
               fPYI1 = pyi1;
               fPYI2 = pyi2;
          }
          else
          {
               fPXI1 = pxi2;
               fPXI2 = pxi1;
               fPYI1 = pyi2;
               fPYI2 = pyi1;
          }
          if (z3 > z4)
          {
               fPXO1 = pxo1;
               fPXO2 = pxo2;
               fPYO1 = pyo1;
               fPYO2 = pyo2;
          }
          else
          {
               fPXO1 = pxo2;
               fPXO2 = pxo1;
               fPYO1 = pyo2;
               fPYO2 = pyo1;
          }
     }

     if (z1 < z4)
     {
          if (z1 < z2)
          {
               fPXI1 = pxi1;
               fPXI2 = pxi2;
               fPYI1 = pyi1;
               fPYI2 = pyi2;
          }
          else
          {
               fPXI1 = pxi2;
               fPXI2 = pxi1;
               fPYI1 = pyi2;
               fPYI2 = pyi1;
          }
          if (z3 < z4)
          {
               fPXO1 = pxo1;
               fPXO2 = pxo2;
               fPYO1 = pyo1;
               fPYO2 = pyo2;
          }
          else
          {
               fPXO1 = pxo2;
               fPXO2 = pxo1;
               fPYO1 = pyo2;
               fPYO2 = pyo1;
          }
     }

     auto dataFlag = CheckPointData();
     if (!dataFlag)
     {
          return false;
     }
     Calculate3DPoint();
     return true;
}

bool TImagingData::CheckPointData()
{
     if (TMath::Abs(fPXI1.X()) < 1000 && TMath::Abs(fPXI1.Y()) < 1000)
     {
          if (TMath::Abs(fPXI2.X()) < 1000 && TMath::Abs(fPXI2.Y()) < 1000)
          {
               if (TMath::Abs(fPXO1.X()) < 1000 && TMath::Abs(fPXO1.Y()) < 1000)
               {
                    if (TMath::Abs(fPXO2.X()) < 1000 && TMath::Abs(fPXO2.Y()) < 1000)
                    {
                         if (TMath::Abs(fPYI1.X()) < 1000 && TMath::Abs(fPYI1.Y()) < 1000)
                         {
                              if (TMath::Abs(fPYI2.X()) < 1000 && TMath::Abs(fPYI2.Y()) < 1000)
                              {
                                   if (TMath::Abs(fPYO1.X()) < 1000 && TMath::Abs(fPYO1.Y()) < 1000)
                                   {
                                        if (TMath::Abs(fPYO2.X()) < 1000 && TMath::Abs(fPYO2.Y()) < 1000)
                                        {
                                             return true;
                                        }
                                   }
                              }
                         }
                    }
               }
          }
     }

     return false;
}

void TImagingData::Calculate3DPoint()
{
     Combine3DPoints(fPXI1, fPXI2, fPYI1, fPYI2, fP3I1, fP3I2);
     Combine3DPoints(fPXO1, fPXO2, fPYO1, fPYO2, fP3O1, fP3O2);
}

bool TImagingData::Combine3DPoints(const TImagingType::Point2D &px1, const TImagingType::Point2D &px2, const TImagingType::Point2D &py1, const TImagingType::Point2D &py2, TImagingType::Point3D &p1, TImagingType::Point3D &p2)
{
     // Calculate line vec
     Vector2D vecIX, vecOX, vecIY, vecOY;
     vecIX = px2 - px1;
     vecIY = py2 - py1;

     // Calculate mean value of height, set as 3D point Z coordinate value
     double zi1, zi2;
     zi1 = px1.Y() / 2.0 + py1.Y() / 2.0;
     zi2 = px2.Y() / 2.0 + py2.Y() / 2.0;
     p1.SetZ(zi1);
     p2.SetZ(zi2);

     // Calulate other 2 coordinates
     Point2D pxi1 = px1 + vecIX.Unit() * (zi1 - px1.Y()) / TMath::Abs(vecIX.Unit() * TVector2(0, 1));
     Point2D pyi1 = py1 + vecIY.Unit() * (zi1 - py1.Y()) / TMath::Abs(vecIY.Unit() * TVector2(0, 1));
     p1.SetX(pxi1.X());
     p1.SetY(pyi1.X());
     Point2D pxi2 = px2 + vecIX.Unit() * (zi2 - px2.Y()) / TMath::Abs(vecIX.Unit() * TVector2(0, 1));
     Point2D pyi2 = py2 + vecIY.Unit() * (zi2 - py2.Y()) / TMath::Abs(vecIY.Unit() * TVector2(0, 1));
     p2.SetX(pxi2.X());
     p2.SetY(pyi2.X());

     return true;
}

bool TImagingData::Combine3DPoints(const TPreImagingDataPart &pre, TImagingType::Point3D &p1, TImagingType::Point3D &p2)
{
     Combine3DPoints(pre.fX1, pre.fX2, pre.fY1, pre.fY2, p1, p2);
}

ostream &TImagingData::Show(ostream &os) const
{
     cout << "################################" << endl;
     cout << "Incident 1: " << endl;
     cout << "(" << fPXI1.X() << ", ?, " << fPXI1.Y() << ") + ( ?, " << fPYI1.X() << ", " << fPYI1.Y() << ") -> "
          << "( " << fP3I1.X() << ", " << fP3I1.Y() << ", " << fP3I1.Z() << ")" << endl;
     cout << "Incident 2: " << endl;
     cout << "(" << fPXI2.X() << ", ?, " << fPXI2.Y() << ") + ( ?, " << fPYI2.X() << ", " << fPYI2.Y() << ") -> "
          << "( " << fP3I2.X() << ", " << fP3I2.Y() << ", " << fP3I2.Z() << ")" << endl;
     cout << "Out 1: " << endl;
     cout << "(" << fPXO1.X() << ", ?, " << fPXO1.Y() << ") + ( ?, " << fPYO1.X() << ", " << fPYO1.Y() << ") -> "
          << "( " << fP3O1.X() << ", " << fP3O1.Y() << ", " << fP3O1.Z() << ")" << endl;
     cout << "Out 2: " << endl;
     cout << "(" << fPXO2.X() << ", ?, " << fPXO2.Y() << ") + ( ?, " << fPYO2.X() << ", " << fPYO2.Y() << ") -> "
          << "( " << fP3O2.X() << ", " << fP3O2.Y() << ", " << fP3O2.Z() << ")" << endl;
     cout << "################################" << endl;

     return os;
}

double TImagingData::GetScatteringAngle() const
{
     return (fP3I2 - fP3I1).Angle(fP3O2 - fP3O1);
}

TPreImagingDataPart::TPreImagingDataPart(double x1, double hx1, double y1, double hy1, double x2, double hx2, double y2, double hy2)
{
     SetPoint(x1, hx1, y1, hy1, x2, hx2, y2, hy2);
}

TPreImagingDataPart::TPreImagingDataPart(const Point2D &x1, const Point2D &x2, const Point2D &y1, const Point2D &y2)
{
     SetPoint(x1, y1, x2, y2);
}

bool TPreImagingDataPart::SetPoint(double x1, double hx1, double y1, double hy1, double x2, double hx2, double y2, double hy2)
{
     fX1.SetX(x1);
     fX1.SetY(hx1);
     fY1.SetX(y1);
     fY1.SetY(hy1);
     fX2.SetX(x2);
     fX2.SetY(hx2);
     fY2.SetX(y2);
     fY2.SetY(hy2);
}

bool TPreImagingDataPart::SetPoint(const Point2D &x1, const Point2D &x2, const Point2D &y1, const Point2D &y2)
{
     fX1 = x1;
     fY1 = y1;
     fX2 = x2;
     fY2 = y2;
}

double TImagingDataProcessor::CalculatePathLength(const TImagingData &imageData, double z1, double z2, double zTargetZone)
{
     double dz1, dz2, dz3, dz4;
     dz1 = imageData.I1().Z();
     dz2 = imageData.I2().Z();
     dz3 = imageData.O1().Z();
     dz4 = imageData.O2().Z();

     if (zTargetZone < TMath::Min(dz2, dz3) || zTargetZone > TMath::Max(dz2, dz3))
     {
          cout << "Error while input target zone. Set target zone automatically." << endl;
          zTargetZone = (dz2 + dz3) / 2.0;
     }

     TVector3 pPoCA, pc, qc;
     auto flag = PoCAManager::sCalcPoCA(pPoCA, pc, qc, imageData);
     double length;

     TVector3 zUnit = {0, 0, 1};

     TVector3 line1;
     TVector3 line2;
     TVector3 pointOrigin1;
     TVector3 pointOrigin2;
     TVector3 pcut1;
     TVector3 pcut2;

     double zChange;

     // if PoCA is calculated to be valid, then PoCA is change point; else point a change point (usually is set to be object region)
     if (!flag)
     {
          zChange = zTargetZone;
     }
     else if (flag)
     {
          zChange = pPoCA.z();
     }

     if ((z1 - zChange) * (dz1 - zChange) >= 0)
     {
          line1 = imageData.I2() - imageData.I1();
          line1 = line1.Unit();
          pointOrigin1 = imageData.I2();
          pcut1 = pc;
     }
     else
     {
          line1 = imageData.O2() - imageData.O1();
          line1 = line1.Unit();
          pointOrigin1 = imageData.O2();
          pcut1 = qc;
     }

     if ((z2 - zChange) * (dz1 - zChange) >= 0)
     {
          line2 = imageData.I2() - imageData.I1();
          line2 = line2.Unit();
          pointOrigin2 = imageData.I2();
          pcut2 = pc;
     }
     else
     {
          line2 = imageData.O2() - imageData.O1();
          line2 = line2.Unit();
          pointOrigin2 = imageData.O2();
          pcut2 = qc;
     }

     // If PoCA is invalid, reset cut points
     if (!flag)
     {
          pcut1 = pointOrigin1 + (zChange - pointOrigin1.Z()) / TMath::Abs(line1.Dot(zUnit)) * line1;
          pcut2 = pointOrigin2 + (zChange - pointOrigin2.Z()) / TMath::Abs(line2.Dot(zUnit)) * line2;
     }

     TVector3 pz1 = pointOrigin1 + (z1 - pointOrigin1.Z()) / TMath::Abs(line1.Dot(zUnit)) * line1;
     TVector3 pz2 = pointOrigin1 + (z2 - pointOrigin1.Z()) / TMath::Abs(line1.Dot(zUnit)) * line1;

     // z1 and z2 are in the same side of PoCA
     if ((z1 - pcut1.Z()) * (z2 - pcut2.Z()) >= 0)
     {
          length = (pz1 - pz2).Mag();
     }
     else
     {
          length = (pz1 - pcut1).Mag() + (pz2 - pcut2).Mag() + (pcut1 - pcut2).Mag();
          // length = (pz1 - pcut1).Mag() + (pz2 - pcut2).Mag();
          // length = (pz1 - pz2).Mag();
     }


     return length;
}

double TImagingDataProcessor::CalculateVelocity(const TImagingData &imageData, double z1, double z2, double zTargetZone, double t1, double t2)
{
     double length = CalculatePathLength(imageData, z1, z2, zTargetZone);
     double dt = TMath::Abs(t1 - t2);
     double v = TMath::Abs(length / dt);
     return v;
}

double TImagingDataProcessor::CalculateGamma(const TImagingData &imageData, double z[4], double t[4], double zTargetZone, double &betaIn)
{
     double v1 = CalculateVelocity(imageData, z[0], z[2], zTargetZone, t[0], t[2]);
     double v2 = CalculateVelocity(imageData, z[1], z[3], zTargetZone, t[1], t[3]);
     double v3 = CalculateVelocity(imageData, z[0], z[3], zTargetZone, t[0], t[3]);

     double vAve = (v1 + v2 + v3) / 3.0;
     double beta = vAve * 1e7 / TMath::C();
     betaIn = beta;
     double gamma = 1.0 / TMath::Sqrt(1 - beta * beta);
     return gamma;
}

double TImagingDataProcessor::CalculateGamma(const TImagingData &imageData, double t[4], double zTargetZone, double &betaIn)
{
     double z[4];
     z[0] = imageData.I1().Z();
     z[1] = imageData.I2().Z();
     z[2] = imageData.O1().Z();
     z[3] = imageData.O2().Z();
     return CalculateGamma(imageData, z, t, zTargetZone, betaIn);
}