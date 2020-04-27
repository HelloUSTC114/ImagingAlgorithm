#include "PoCA.h"

using namespace POCAType;

PoCAManager::PoCAManager()
{
    fXi = fYi = fZi = fXf = fYf = fZf = 0;
    fNbinsX = fNbinsY = fNbinsZ = 0;
    fVoxelCount = 0;

    fPathArray = NULL;
    fDensityArray = NULL;
    fMuCountArray = NULL;
}

PoCAManager::~PoCAManager()
{
    Clear();
}

void PoCAManager::Clear()
{
    fXi = fYi = fZi = fXf = fYf = fZf = 0;
    fNbinsX = fNbinsY = fNbinsZ = 0;
    fVoxelCount = 0;

    delete[] fPathArray;
    delete[] fDensityArray;
    delete[] fMuCountArray;
    fPathArray = NULL;
    fDensityArray = NULL;
    fMuCountArray = NULL;
}

bool PoCAManager::Init(const POCAType::Point3 &pi, const POCAType::Point3 &pf, int nbinsX, int nbinsY, int nbinsZ)
{
    Clear();
    fXi = TMath::Min(pi.X(), pf.X());
    fYi = TMath::Min(pi.Y(), pf.Y());
    fZi = TMath::Min(pi.Z(), pf.Z());
    fXf = TMath::Max(pi.X(), pf.X());
    fYf = TMath::Max(pi.Y(), pf.Y());
    fZf = TMath::Max(pi.Z(), pf.Z());

    fNbinsX = nbinsX;
    fNbinsY = nbinsY;
    fNbinsZ = nbinsZ;
    if (fNbinsX < 0 || fNbinsY < 0 || fNbinsZ < 0)
    {
        Clear();
        return false;
    }

    fVoxelCount = fNbinsX * fNbinsY * fNbinsZ;
    fMuCountArray = new long[fVoxelCount];
    fDensityArray = new float[fVoxelCount];
    fPathArray = new float[fVoxelCount];

    for (int i = 0; i < fVoxelCount; i++)
    {
        fPathArray[i] = 0;
        fMuCountArray[i] = 0;
        fDensityArray[i] = 0;
    }
    return true;
}

long PoCAManager::GetVoxelID(Index3 &idx, const Point3 &pNow)
{
    bool flagCheck = CheckPointInside(pNow);
    if (!flagCheck)
    {
        return -1;
    }

    long vID;
    int idxX, idxY, idxZ;
    idxX = (pNow.X() - fXi) * fNbinsX / (fXf - fXi);
    idxY = (pNow.Y() - fYi) * fNbinsY / (fYf - fYi);
    idxZ = (pNow.Z() - fZi) * fNbinsZ / (fZf - fZi);

    // Avoid boundary
    // if the point is at the boundary, count it as idxX--
    if (pNow.X() == fXf)
    {
        idxX--;
    }
    if (pNow.Y() == fYf)
    {
        idxY--;
    }
    if (pNow.Z() == fZf)
    {
        idxZ--;
    }
    idx.SetXYZ(idxX, idxY, idxZ);
    vID = idxX + idxY * fNbinsX + idxZ * fNbinsX * fNbinsY;
    return vID;
}

long PoCAManager::GetVoxelID(const Point3 &pNow)
{
    Index3 idxTemp;
    return GetVoxelID(idxTemp, pNow);
}

bool PoCAManager::CheckPointInside(const Point3 &pNow)
{
    if (pNow.X() < XMin() || pNow.X() > XMax())
    {
        return false;
    }
    if (pNow.Y() < YMin() || pNow.Y() > YMax())
    {
        return false;
    }
    if (pNow.Z() < ZMin() || pNow.Z() > ZMax())
    {
        return false;
    }
    return true;
}

bool PoCAManager::CalcPath(const Point3 &pi, const Point3 &pf, float *ev_path)
{

    Index3 idx_pi, idx_pf;
    long voxelID = 0;
    voxelID = GetVoxelID(idx_pi, pi);
    if (voxelID < 0)
        return false;
    voxelID = GetVoxelID(idx_pf, pf);
    if (voxelID < 0)
        return false;

    // line vec: Unit vector
    TVector3 lineVec = pf - pi;
    lineVec = lineVec.Unit();

    map<float, const Point3 *> mapNode;
    mapNode[pi.Z()] = &pi;
    mapNode[pf.Z()] = &pf;

    // Locate all points through which line crosses voxel boundary
    const Point3 *pNow = NULL;
    for (int idxNow = min(idx_pi.X(), idx_pf.X()) + 1; idxNow <= max(idx_pi.X(), idx_pf.X()); idxNow++)
    {
        float xNow = fXi + idxNow * dX();
        // 3D position for this voxel
        auto pTemp = new Point3;
        *pTemp = (pi + (xNow - pi.X()) / lineVec.X() * lineVec);
        pNow = pTemp;
        mapNode[pNow->Z()] = pNow;
    }
    for (int idxNow = min(idx_pi.Y(), idx_pf.Y()) + 1; idxNow <= max(idx_pi.Y(), idx_pf.Y()); idxNow++)
    {
        float yNow = fYi + idxNow * dY();
        // 3D position for this voxel
        auto pTemp = new Point3;
        *pTemp = (pi + (yNow - pi.Y()) / lineVec.Y() * lineVec);
        pNow = pTemp;
        mapNode[pNow->Z()] = pNow;
    }
    for (int idxNow = min(idx_pi.Z(), idx_pf.Z()) + 1; idxNow <= max(idx_pi.Z(), idx_pf.Z()); idxNow++)
    {
        float zNow = fZi + idxNow * dZ();
        // 3D position for this voxel
        auto pTemp = new Point3;
        *pTemp = (pi + (zNow - pi.Z()) / lineVec.Z() * lineVec);
        pNow = pTemp;
        mapNode[pNow->Z()] = pNow;
    }

    // Calc Path Segement
    const Point3 *pNext = NULL;
    for (auto itr = mapNode.begin(); itr != mapNode.end() && pNext != &pf; itr++)
    {
        if (itr == mapNode.begin())
        {
            pNow = &pi;
            continue;
        }
        pNext = itr->second;
        voxelID = GetVoxelID(0.5 * (*pNow + *pNext));
        if (voxelID < 0)
            return false;
        ev_path[voxelID] += TMath::Abs((pNext->Z() - pNow->Z()) / lineVec.Z());
        pNow = pNext;
    }
    for (auto itr = mapNode.begin(); itr != mapNode.end(); itr++)
    {
        if (itr->second == &pi || itr->second == &pf)
        {
            continue;
        }
        delete itr->second;
        itr->second = NULL;
    }

    return true;
}

double PoCAManager::CalcPoCA(TVector3 &pPoCA, const TImagingData &posData)
{
    const TVector3 vPos[4]{posData.I1(), posData.I2(), posData.O1(), posData.O2()};
    return CalcPoCA(pPoCA, vPos);
}

bool PoCAManager::sCalcPoCA(TVector3 &pPoCA, TVector3 &pcout, TVector3 &qcout, const TVector3 *vPos)
{
    const TVector3 &p0 = vPos[0];
    const TVector3 &p1 = vPos[1];
    const TVector3 &q0 = vPos[2];
    const TVector3 &q1 = vPos[3];

    TVector3 u = p1 - p0;
    TVector3 v = q1 - q0;
    TVector3 w = p0 - q0;

    Float_t angle = u.Angle(v);
    angle = std::atan(v.X() / v.Z()) - std::atan(u.X() / u.Z());
    // if (fabs(angle) > ANGLE_LIMITED)
    Float_t a = u * u;
    Float_t b = u * v;
    Float_t c = v * v;
    Float_t d = u * w;
    Float_t e = v * w;

    TVector3 pc = p0 + ((b * e - c * d) / (a * c - b * b)) * u;
    TVector3 qc = q0 + ((a * e - b * d) / (a * c - b * b)) * v;

    pcout = pc;
    qcout = qc;

    pPoCA.SetX(pc.x() / 2.0 + qc.x() / 2.0);
    pPoCA.SetY(pc.y() / 2.0 + qc.y() / 2.0);
    pPoCA.SetZ(pc.z() / 2.0 + qc.z() / 2.0);

    bool flag = 1;

    // Judge whether z of pPoCA is inside valid region
    double z1, z2;
    if (vPos[0].z() > vPos[3].z())
    {
        z1 = TMath::Min(vPos[0].z(), vPos[1].z());
        z2 = TMath::Max(vPos[2].z(), vPos[3].z());
    }
    else if (vPos[0].z() < vPos[3].z())
    {
        z1 = TMath::Max(vPos[0].z(), vPos[1].z());
        z2 = TMath::Min(vPos[2].z(), vPos[3].z());
    }
    else
    {
        flag = 0;
        return flag;
    }

    if (pPoCA.z() < TMath::Max(z1, z2) && pPoCA.z() > TMath::Min(z1, z2))
        flag = 1;
    else
        flag = 0;

    return flag;
}

bool PoCAManager::sCalcPoCA(TVector3 &pPoCA, TVector3 &pcout, TVector3 &qcout, const TImagingData &posData)
{
    const TVector3 vPos[4]{posData.I1(), posData.I2(), posData.O1(), posData.O2()};
    return sCalcPoCA(pPoCA, pcout, qcout, vPos);
}

double PoCAManager::CalcPoCA(TVector3 &pPoCA, const TVector3 *vPos)
{
    const TVector3 &p0 = vPos[0];
    const TVector3 &p1 = vPos[1];
    const TVector3 &q0 = vPos[2];
    const TVector3 &q1 = vPos[3];

    TVector3 u = p1 - p0;
    TVector3 v = q1 - q0;
    TVector3 w = p0 - q0;

    for (Int_t i = 0; i < fVoxelCount; i++)
        fPathArray[i] = 0.0;

    Float_t angle = u.Angle(v);
    angle = std::atan(v.X() / v.Z()) - std::atan(u.X() / u.Z());
    if (fabs(angle) > ANGLE_LIMITED)
    {
        Float_t a = u * u;
        Float_t b = u * v;
        Float_t c = v * v;
        Float_t d = u * w;
        Float_t e = v * w;

        TVector3 pc = p0 + ((b * e - c * d) / (a * c - b * b)) * u;
        TVector3 qc = q0 + ((a * e - b * d) / (a * c - b * b)) * v;

        pPoCA.SetX(pc.x() / 2.0 + qc.x() / 2.0);
        pPoCA.SetY(pc.y() / 2.0 + qc.y() / 2.0);
        pPoCA.SetZ(pc.z() / 2.0 + qc.z() / 2.0);

        Long_t voxelID = GetVoxelID(pPoCA);
        // pc.Print();
        // qc.Print();
        // p0.Print();
        // (((b * e - c * d) / (a * c - b * b)) * u).Print();

        // If pc or qc or poca is outside valid zone,
        if (voxelID < 0 || GetVoxelID(pc) < 0 || GetVoxelID(qc) < 0)
            return -angle;

        CalcPath(p1, pc, fPathArray);
        // cout << "End p1 pc" << endl
        //  << endl;
        // cout << "Start pc qc" << endl;
        CalcPath(pc, qc, fPathArray);
        // cout << "End pc qc" << endl
        //  << endl;
        // cout << "Start qc q0" << endl;
        CalcPath(qc, q0, fPathArray);
        // cout << endl;

        // Float_t dL = (Z_END - Z_START) / BIN_Z; //fPathArray[voxelID];
        // Float_t dL = (fZf - fZi) / fNbinsZ; //fPathArray[voxelID];
        // if (dL < 0.1)
        //     dL = 0.1;
        if (fPathArray[voxelID] == 0)
        {
            return -angle;
        }

        // TVector3 lineTemp = p0 - q1;
        // double path = (pc-qc).Mag();
        // double path = TMath::Max((double)fPathArray[voxelID], 0.1 * dZ());
        // fDensityArray[voxelID] += angle * angle / path;
        // fDensityArray[voxelID] += angle * angle * TMath::Abs(lineTemp.Unit().Z());
        // fDensityArray[voxelID] += pow(10.0, 6.0) * angle * angle / fPathArray[voxelID];
        fDensityArray[voxelID] += pow(10.0, 6.0) * angle * angle / dZ();
    }
    else
        CalcPath(p1, q0, fPathArray);

    for (Int_t i = 0; i < fVoxelCount; i++)
    {
        if (fPathArray[i] > 0.001 * dL())
        {
            fMuCountArray[i]++;
            // cout << "Path Length: " << fPathArray[i] << endl;
        }
        // else if (fPathArray[i] != 0)
        // {
        //     cout << "dL: " << dL() << endl;
        //     cout << "Test Path length: " << fPathArray[i] << endl;
        // }
    }

    return angle;
}

ostream &PoCAManager::Show(ostream &os)
{
    os << "################POCA imaging algorithm manager################" << endl;

    os << "Imaging region: " << endl;
    os << "(" << fXi << ", " << fYi << ", " << fZi << ") ==> (" << fXf << ", " << fYf << ", " << fZf << ") cm" << endl;
    os << "Bin Number: " << endl;
    os << "(" << fNbinsX << "," << fNbinsY << "," << fNbinsZ << ")" << endl;
    os << "Voxel Number: " << endl;
    os << fVoxelCount << endl;
    os << "Bin Step: " << endl;
    os << "(" << dX() << "," << dY() << "," << dZ() << ") cm" << endl;

    os << "################POCA imaging algorithm manager################" << endl;
    return os;
}

PoCAManager *&PoCAManager::CurrentPoCAManager()
{
    static auto currentPoCAManager = new PoCAManager();
    return currentPoCAManager;
}

TFile *PoCAManager::WriteResult(TFile *file)
{
    if (!fDensityArray)
    {
        return NULL;
    }
    if (file == NULL)
    {
        file = new TFile("ImagingPoCA.root", "recreate");
    }
    auto hDensity = new TH3F("hDensity", "Density Distribution", fNbinsX, XMin(), XMax(), fNbinsY, YMin(), YMax(), fNbinsZ, ZMin(), ZMax());
    auto hMuon = new TH3I("hMuon", "Muon Count Distribution", fNbinsX, XMin(), XMax(), fNbinsY, YMin(), YMax(), fNbinsZ, ZMin(), ZMax());

    for (int i = 0; i < fVoxelCount; i++)
    {
        auto index = ConvertIndex(i);
        cout << "x: " << index.X() << '\t' << "y: " << index.Y() << '\t' << "z: " << index.Z() << '\t' << fDensityArray[i] << '\t' << fMuCountArray[i] << endl;
        double densityTemp = 0;
        if (fMuCountArray[i] > 0)
            densityTemp = fDensityArray[i] / fMuCountArray[i];
        else
            densityTemp = 0.0;

        hDensity->SetBinContent(index.X() + 1, index.Y() + 1, index.Z() + 1, fDensityArray[i]);
        hMuon->SetBinContent(index.X() + 1, index.Y() + 1, index.Z() + 1, fMuCountArray[i]);
    }

    hDensity->Write();
    hMuon->Write();

    return file;
}

TFile *PoCAManager::WritePath(TFile *file)
{
    if (!fDensityArray)
    {
        return NULL;
    }
    file->cd();
    auto hPath = new TH3F("hPath", "Path", fNbinsX, XMin(), XMax(), fNbinsY, YMin(), YMax(), fNbinsZ, ZMin(), ZMax());
    for (int i = 0; i < fVoxelCount; i++)
    {
        auto index = ConvertIndex(i);
        hPath->SetBinContent(index.X() + 1, index.Y() + 1, index.Z() + 1, fPathArray[i]);
    }
    hPath->SetDirectory(file);
    hPath->Write();
    return file;
}

Index3 PoCAManager::ConvertIndex(long arrayIndex)
{
    Index3 idxTemp;
    int idX = arrayIndex % fNbinsX;
    int idY = (arrayIndex / fNbinsX) % fNbinsY;
    int idZ = arrayIndex / fNbinsX / fNbinsY;
    idxTemp.SetXYZ(idX, idY, idZ);
    return idxTemp;
}