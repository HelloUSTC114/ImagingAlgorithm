#ifndef POCA_H
#define POCA_H

#include "TMath.h"
#include "TVector2.h"
#include "TVector3.h"
#include "TH3F.h"
#include "TH3I.h"
#include "TFile.h"

#include <map>
#include <iostream>
#include "TImagingData.h"

#define MAX(x, y) (x > y ? x : y)
#define MIN(x, y) (x < y ? x : y)

namespace POCAType
{
typedef TVector3 Index3;
typedef TVector3 Point3; // Point Coordinate
typedef TVector3 Vector3;

inline double max(double x, double y) { return (x > y) ? x : y; }
inline double min(double x, double y) { return (x < y) ? x : y; }
}; // namespace POCAType

using namespace std;

const static double PI = TMath::Pi();
// const static double PI = 3.14159265358;

const Float_t ANGLE_LIMITED = 0.01;             // rad, Angle Limit
const Int_t BIN_X = 30, BIN_Y = 10, BIN_Z = 10; //
const Float_t X_START = -20., X_END = 20.0;     // cm
const Float_t Y_START = -20., Y_END = 20.;      // cm
const Float_t Z_START = -20., Z_END = 20.;      // cm
const Float_t Z_POS[4] = {-90., -20., 20., 90.};

class PoCAManager
{
public:
    PoCAManager();
    ~PoCAManager();
    bool Init(const POCAType::Point3 &pi, const POCAType::Point3 &pf, int nbinsX, int nbinsY, int nbinsZ); // Unit: cm
    void Clear();

    long GetVoxelID(POCAType::Index3 &idx, const POCAType::Point3 &pNow);
    long GetVoxelID(const POCAType::Point3 &pNow);
    POCAType::Index3 ConvertIndex(long arrayIndex);

    // pi is incident point, pf is shooting point
    bool CalcPath(const POCAType::Point3 &pi, const POCAType::Point3 &pf, float *ev_path);
    double CalcPoCA(TVector3 &pPoCA, const TVector3 *vPos);
    double CalcPoCA(TVector3 &pPoCA, const TImagingData &posData);

    static bool sCalcPoCA(TVector3 &pPoCA, TVector3 &pcout, TVector3 &qcout, const TVector3 *vPos);
    static bool sCalcPoCA(TVector3 &pPoCA, TVector3 &pcout, TVector3 &qcout, const TImagingData &posData);

    bool CheckPointInside(const POCAType::Point3 &pNow);

    // Size for each voxel
    double dX() { return TMath::Abs((fXf - fXi) / (double)fNbinsX); }
    double dY() { return TMath::Abs((fYf - fYi) / (double)fNbinsY); }
    double dZ() { return TMath::Abs((fZf - fZi) / (double)fNbinsZ); }
    double dL() { return TMath::Sqrt(dX() * dX() + dY() * dY() + dZ() * dZ()); }
    // Minimum/Maximum for each dimension
    double XMax() { return POCAType::max(fXi, fXf); }
    double YMax() { return POCAType::max(fYi, fYf); }
    double ZMax() { return POCAType::max(fZi, fZf); }
    double XMin() { return POCAType::min(fXi, fXf); }
    double YMin() { return POCAType::min(fYi, fYf); }
    double ZMin() { return POCAType::min(fZi, fZf); }

    double Xi() { return fXi; }
    double Yi() { return fYi; }
    double Zi() { return fZi; }
    double Xf() { return fXf; }
    double Yf() { return fYf; }
    double Zf() { return fZf; }

    int GetNbinsX() { return fNbinsX; }
    int GetNbinsY() { return fNbinsY; }
    int GetNbinsZ() { return fNbinsZ; }
    long GetVoxelCount() { return fVoxelCount; }

    const float *GetPathArray() { return fPathArray; }
    const float *GetDensityArray() { return fDensityArray; }
    const long *GetMuCountArray() { return fMuCountArray; }

    ostream &Show(ostream &os);

    static PoCAManager *&CurrentPoCAManager();

    TFile *WriteResult(TFile *file = NULL);
    TFile *WritePath(TFile *file);

private:
    // POCAType::Point3
    double fXi; // Imaging start
    double fYi;
    double fZi;
    double fXf;
    double fYf;
    double fZf;

    int fNbinsX; // Total X bin number
    int fNbinsY;
    int fNbinsZ;
    long fVoxelCount;

    float *fPathArray;    // Penetration length through all voxels for one muon event
    float *fDensityArray; // Save summary of density for all voxels
    long *fMuCountArray;  // Save count of muons for all voxels
};

#define gPoCA (PoCAManager::CurrentPoCAManager())

#endif