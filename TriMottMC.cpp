#include <fstream>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <random>
#include <cmath>
#include <cstring>
#define TEMPERATURESEQUENCEMAX 1000
using namespace std;

std::uniform_real_distribution<double> realRand(0,1);
std::uniform_int_distribution<int> intRand(0,65535);
std::default_random_engine engineTime(time(NULL));

class Lattice
{
    static const int LATTICESIZE_MAX = 100;
    static const int NEIGHBOR_MAX = 10;

    public:
        Lattice(int row, int column, string latType, double interactionJ1, double interactionJ2, double externalXElecField, double externalYElecField,
        double externalMagField, string interactionType, double DeducedIsingDirection[4]);

        int m_size1;
        int m_size2;

        void SetNeighborIndex(string latType);
        void InitSpin();
        void ShowConfig();
        void SweepFlip();

        double GetEnergy()
        {
            return configE;
        }

        double GetMagnetization()
        {
            configM = 0;
            // if (m_interactionType == "Heisenberg")
            // {
            //     for (int i=0; i<m_size1*m_size2; i++)
            //     {
            //         configM += GetSiteZSpinComponent(i);
            //     }
            // }
            // else if (m_interactionType == "Ising")
            // {
            //     for (int i=0; i<m_size1*m_size2; i++)
            //     {
            //         configM += spinConfig[i][0];
            //     }
            // }
            for (int i=0; i<m_size1*m_size2; i++)
            {
                configM += GetSiteZSpinComponent(i);
            }
            return configM;
        }

        double GetPolarizationX()
        {
            configPx = 0;
            for (int i=0; i<m_size1*m_size2; i++)
            {
                configPx += GetSiteXPseudospinComponent(i);
            }
            return configPx;
        }

        double GetPolarizationY()
        {
            configPy = 0;
            for (int i=0; i<m_size1*m_size2; i++)
            {
                configPy += GetSiteYPseudospinComponent(i);
            }
            return configPy;
        }

        void SetTemperature(double assignedTemperature)
        {
            m_temperature = assignedTemperature;
        }

        void DebugInput(int oritemp)//test hamiltonian
        {
            cout << "Debug Information:" << endl;
            for (int i=0; i<floor(log2(oritemp)/2)+1; i++)
            {
                if (oritemp>>(2*i)&1)
                {
                    spinConfig[i][0] = m_deducedIsingDirection[0];
                    spinConfig[i][1] = m_deducedIsingDirection[1];
                }
                else
                {
                    spinConfig[i][0] = M_PI - m_deducedIsingDirection[0];
                    spinConfig[i][1] = M_PI + m_deducedIsingDirection[1];
                }
                if (oritemp>>(2*i+1)&1)
                {
                    spinConfig[i][2] = m_deducedIsingDirection[2];
                    spinConfig[i][3] = m_deducedIsingDirection[3];
                }
                else
                {
                    spinConfig[i][2] = M_PI - m_deducedIsingDirection[2];
                    spinConfig[i][3] = M_PI + m_deducedIsingDirection[3];
                }
            }
            cout << spinConfig[0][0] << spinConfig[0][1] << spinConfig[0][2] << spinConfig[0][3] << endl;
            cout << GetSiteHamiltonian(0) << endl;
        }

        void DebugOutput()
        {
            ShowConfig();
            // cout << GetSiteInnerProduct(0,1) << endl;
        }

        void GetBinaryConfigIsingEMP(int binaryCode, double result[2]);

        // ~Lattice()
        // {
        //     delete temp_stackedIndex;
        //     delete temp_randomSphericalDirection;
        //     delete temp_selectedSiteSpin;
        // }

    private:
        int m_neighborNum;
        double m_interactionJ1;
        double m_interactionJ2;
        double m_temperature = 1;
        double m_externalMagField;
        double m_externalXElecField;
        double m_externalYElecField;
        string m_interactionType;
        double m_deducedIsingDirection[4];
        double m_s = 0.5;
        double m_tau = 0.5;

        double configE;
        double deltaE;
        double configM;
        double configPx;
        double configPy;

        double spinConfig[LATTICESIZE_MAX * LATTICESIZE_MAX][4];
        int neighbor[LATTICESIZE_MAX * LATTICESIZE_MAX][NEIGHBOR_MAX];
        double temp_anisotropic[3][2] = {M_PI/2, 0, M_PI/2, M_PI*4/3, M_PI/2, M_PI*8/3};

        int temp_stackedIndex[2];
        double temp_siteE;
        int temp_orientation;
        double temp_randomSphericalDirection[2];
        int temp_flipPos;
        double temp_selectedSiteSpin[4];
        bool isComponentSFlipped;


        int GetFlattenedCoordinate(int stackedIndex1, int stackedIndex2)
        {
            if (stackedIndex1 < 0)  stackedIndex1 += m_size1;
            if (stackedIndex1 >= m_size1)  stackedIndex1 -= m_size1;
            if (stackedIndex2 < 0)  stackedIndex2 += m_size2;
            if (stackedIndex2 >= m_size2)  stackedIndex2 -= m_size2;

            return stackedIndex1 * m_size1 + stackedIndex2;
        }

        void GetStackedCoordinate(int flattenedIndex)
        {
            temp_stackedIndex[0] = int(flattenedIndex / m_size1);
            temp_stackedIndex[1] = flattenedIndex % m_size1;
        }

        void SetRandomSphericalDirection()
        {
            temp_randomSphericalDirection[0] = acos(1-2*realRand(engineTime));
            temp_randomSphericalDirection[1] = 2 * M_PI * realRand(engineTime);
        }

        double GetSiteHamiltonian(int site)
        {
            temp_siteE = 0;
            for (int i=0; i<3; i++)
            {
                temp_siteE +=
                    m_s * m_s * GetSphericalInnerProduct(spinConfig[site][0],spinConfig[site][1],
                                spinConfig[neighbor[site][2*i]][0],spinConfig[neighbor[site][2*i]][1])
                    * (1 - 4 * m_tau * GetSphericalInnerProduct(spinConfig[site][2],spinConfig[site][3],
                                        temp_anisotropic[i][0],temp_anisotropic[i][1])
                        + 2 * m_tau * GetSphericalInnerProduct(spinConfig[neighbor[site][2*i]][2],spinConfig[neighbor[site][2*i]][3],
                                        temp_anisotropic[i][0],temp_anisotropic[i][1])
                        - 8 * m_tau * m_tau * GetSphericalInnerProduct(spinConfig[neighbor[site][2*i]][2],spinConfig[neighbor[site][2*i]][3],
                                        temp_anisotropic[i][0],temp_anisotropic[i][1])
                            * GetSphericalInnerProduct(spinConfig[site][2],spinConfig[site][3],temp_anisotropic[i][0],temp_anisotropic[i][1]));
            }
            temp_siteE *= 2.0 / 9.0 * m_interactionJ1;

            // temp_siteE += m_interactionJ2 / 3 + interlayer
            temp_siteE -= m_externalMagField * m_s * GetSiteZSpinComponent(site); // g\mu_B = 1
            return temp_siteE;
        }

        double GetSiteZSpinComponent(int site)
        {
            return cos(spinConfig[site][0]);
        }

        double GetSiteXPseudospinComponent(int site)
        {
            return sin(spinConfig[site][2]) * cos(spinConfig[site][3]);
        }

        double GetSiteYPseudospinComponent(int site)
        {
            return sin(spinConfig[site][2]) * sin(spinConfig[site][3]);
        }

        double GetSiteInnerProduct(int site1, int site2)
        {
            return sin(spinConfig[site1][0]) * sin(spinConfig[site2][0]) * cos(spinConfig[site1][1]-spinConfig[site2][1])
                 + cos(spinConfig[site1][0]) * cos(spinConfig[site2][0]);
        }

        double GetSphericalInnerProduct(double theta1, double phi1, double theta2, double phi2)
        {
            // cout << sin(theta1) * sin(theta2) * cos(phi1-phi2) + cos(theta1) * cos(theta2) << endl;
            return sin(theta1) * sin(theta2) * cos(phi1-phi2) + cos(theta1) * cos(theta2);
        }
};

Lattice::Lattice(int row, int column, string latType, double interactionJ1, double interactionJ2, double externalXElecField, double externalYElecField,
     double externalMagField=0, string interactionType="Ising", double deducedIsingDirection[4]=nullptr)
    : m_size1(column), m_size2(row), m_interactionJ1(interactionJ1), m_interactionJ2(interactionJ2), m_externalMagField(externalMagField),
        m_externalXElecField(externalXElecField), m_externalYElecField(externalYElecField), m_interactionType(interactionType)
{
    if (m_deducedIsingDirection == nullptr)
    {
        double temp_Zaxis[4] = {0.0};
        *m_deducedIsingDirection = *temp_Zaxis;
    }
    else
    {
        for (int i=0; i<4; i++)
        {
            m_deducedIsingDirection[i] = deducedIsingDirection[i];
        }
    }
    SetNeighborIndex(latType);
    InitSpin();
}

void Lattice::SetNeighborIndex(string latType)  // i/j not necessarily correspond to x/y for spin/lattice?
{
    if (latType == "Square")
    {
        m_neighborNum = 4;
        for (int i=0; i<m_size1; i++)
        {
            for (int j=0; j<m_size2; j++)
            {
                neighbor[i*m_size1+j][0] = Lattice::GetFlattenedCoordinate(i,j+1);
                neighbor[i*m_size1+j][1] = Lattice::GetFlattenedCoordinate(i-1,j);
                neighbor[i*m_size1+j][2] = Lattice::GetFlattenedCoordinate(i,j-1);
                neighbor[i*m_size1+j][3] = Lattice::GetFlattenedCoordinate(i+1,j);
            }
        }
    }
    else if (latType == "Hexagon")
    {
        m_neighborNum = 3;
        for (int i=0; i<m_size1; i++)
        {
            for (int j=0; j<m_size2; j++)
            {
                neighbor[i*m_size1+j][0] = Lattice::GetFlattenedCoordinate(i,j+1);
                neighbor[i*m_size1+j][1] = Lattice::GetFlattenedCoordinate(i,j-1);
                if (i % 2 == 0)
                {
                    neighbor[i*m_size1+j][2] = Lattice::GetFlattenedCoordinate(i+1,j+1);
                }
                else
                {
                    neighbor[i*m_size1+j][2] = Lattice::GetFlattenedCoordinate(i-1,j-1);
                }
            }
        }
    }
    else if (latType == "Triangle")
    {
        m_neighborNum = 6;
        for (int i=0; i<m_size1; i++)
        {
            for (int j=0; j<m_size2; j++)
            {
                neighbor[i*m_size1+j][0] = Lattice::GetFlattenedCoordinate(i-1,j);
                neighbor[i*m_size1+j][1] = Lattice::GetFlattenedCoordinate(i,j+1);
                neighbor[i*m_size1+j][2] = Lattice::GetFlattenedCoordinate(i+1,j+1);
                neighbor[i*m_size1+j][3] = Lattice::GetFlattenedCoordinate(i+1,j);
                neighbor[i*m_size1+j][4] = Lattice::GetFlattenedCoordinate(i,j-1);
                neighbor[i*m_size1+j][5] = Lattice::GetFlattenedCoordinate(i-1,j-1);
            }
        }
    }
    // Throw exception
}

void Lattice::InitSpin()    // Heisenberg ConfigE
{
    configE = 0;
    if (m_interactionType == "Heisenberg")
    {
        for (int i=0; i<m_size1; i++)
        {
            for (int j=0; j<m_size2; j++)
            {
                Lattice::SetRandomSphericalDirection();
                spinConfig[i*m_size1+j][0] = temp_randomSphericalDirection[0];
                spinConfig[i*m_size1+j][1] = temp_randomSphericalDirection[1];
            }
        }

        for (int i=0; i<m_size1; i++)
        {
            for (int j=0; j<m_size2; j++)
            {
                for (int k=0; k<m_neighborNum; k++)
                {
                    configE += m_interactionJ1 * Lattice::GetSiteInnerProduct(i*m_size1+j,neighbor[i*m_size1+j][k]) / 2;
                }
                configE -= m_externalMagField * Lattice::GetSiteZSpinComponent(i*m_size1+j);
            }
        }
    }

    else if (m_interactionType == "Ising")
    {
        for (int i=0; i<m_size1; i++)
        {
            for (int j=0; j<m_size2; j++)
            {
                temp_orientation = intRand(engineTime) % 4;
                if (temp_orientation>>0&1)
                {
                    spinConfig[i*m_size1+j][0] = m_deducedIsingDirection[0];
                    spinConfig[i*m_size1+j][1] = m_deducedIsingDirection[1];
                }
                else
                {
                    spinConfig[i*m_size1+j][0] = M_PI - m_deducedIsingDirection[0];
                    spinConfig[i*m_size1+j][1] = M_PI + m_deducedIsingDirection[1];
                }
                if (temp_orientation>>1&1)
                {
                    spinConfig[i*m_size1+j][2] = m_deducedIsingDirection[2];
                    spinConfig[i*m_size1+j][3] = m_deducedIsingDirection[3];
                }
                else
                {
                    spinConfig[i*m_size1+j][2] = M_PI - m_deducedIsingDirection[2];
                    spinConfig[i*m_size1+j][3] = M_PI + m_deducedIsingDirection[3];
                }
            }
        }

        for (int i=0; i<m_size1; i++)
        {
            for (int j=0; j<m_size2; j++)
            {
                // for (int k=0; k<m_neighborNum; k++)
                // {
                //     configE += m_interactionJ1 * spinConfig[i*m_size1+j][0] * spinConfig[neighbor[i*m_size1+j][k]][0] / 2;
                // }
                // configE -= m_externalMagField * spinConfig[i*m_size1+j][0];
                configE += GetSiteHamiltonian(i*m_size1+j);
            }
        }
    }
    // Throw exception
}

void Lattice::ShowConfig()
{
    cout << m_interactionType << ", Config:" << endl;
    for (int i=0; i<m_size1; i++)
    {
        for (int j=0; j<m_size2; j++)
        {
            cout << "(" << spinConfig[i*m_size1+j][0] << "," << spinConfig[i*m_size1+j][1] << ")|";
        }
        cout << endl;
    }
    cout << GetEnergy() << endl;
}

void Lattice::SweepFlip()   // deltaE remember: 6 sites
{
    for (int i=0; i<m_size1; i++)
    {
        for (int j=0; j<m_size2; j++)
        {
            deltaE = 0;

            if (m_interactionType == "Heisenberg")
            {
                temp_flipPos = i * m_size1 + j;
                temp_selectedSiteSpin[0] = spinConfig[temp_flipPos][0];
                temp_selectedSiteSpin[1] = spinConfig[temp_flipPos][1];

                Lattice::SetRandomSphericalDirection();

                for (int k=0; k<m_neighborNum; k++)
                {
                    deltaE -= m_interactionJ1 * Lattice::GetSiteInnerProduct(temp_flipPos,neighbor[temp_flipPos][k]);
                }
                deltaE += m_externalMagField * Lattice::GetSiteZSpinComponent(temp_flipPos);

                spinConfig[temp_flipPos][0] = temp_randomSphericalDirection[0];
                spinConfig[temp_flipPos][1] = temp_randomSphericalDirection[1];

                for (int k=0; k<m_neighborNum; k++)
                {
                    deltaE += m_interactionJ1 * Lattice::GetSiteInnerProduct(temp_flipPos,neighbor[temp_flipPos][k]);
                }
                deltaE -= m_externalMagField * Lattice::GetSiteZSpinComponent(temp_flipPos);

                if (realRand(engineTime) < exp(-deltaE / m_temperature))
                {
                    // spinConfig[flipPos] = -spinConfig[flipPos];
                    configE += deltaE;
                }
                else
                {
                    spinConfig[temp_flipPos][0] = temp_selectedSiteSpin[0];
                    spinConfig[temp_flipPos][1] = temp_selectedSiteSpin[1];
                }
            }

            else if (m_interactionType == "Ising")
            {
                temp_flipPos = intRand(engineTime) % ((m_size1*m_size2)*2);
                isComponentSFlipped = temp_flipPos % 2;
                temp_flipPos /= 2;

                for (int k=0; k<4; k++)
                {
                    temp_selectedSiteSpin[k] = spinConfig[temp_flipPos][k];
                }

                deltaE -= GetSiteHamiltonian(temp_flipPos);
                for (int k=0; k<3; k++)  // remember to modify when add layers
                {
                    deltaE -= GetSiteHamiltonian(neighbor[temp_flipPos][2*k+1]);
                }

                if (isComponentSFlipped)
                {
                    spinConfig[temp_flipPos][0] = M_PI - spinConfig[temp_flipPos][0];
                    spinConfig[temp_flipPos][1] = M_PI + spinConfig[temp_flipPos][1];
                }
                else
                {
                    spinConfig[temp_flipPos][2] = M_PI - spinConfig[temp_flipPos][2];
                    spinConfig[temp_flipPos][3] = M_PI + spinConfig[temp_flipPos][3]; 
                }

                deltaE += GetSiteHamiltonian(temp_flipPos);
                for (int k=0; k<3; k++)  // remember to modify when add layers
                {
                    deltaE += GetSiteHamiltonian(neighbor[temp_flipPos][2*k+1]);
                }

                if (realRand(engineTime) < exp(-deltaE / m_temperature))
                {
                    if (spinConfig[temp_flipPos][1] > 2 * M_PI)                 // control M_PI + value
                    {
                        spinConfig[temp_flipPos][1] -= 2 * M_PI;
                    }
                    if (spinConfig[temp_flipPos][3] > 2 * M_PI)
                    {
                        spinConfig[temp_flipPos][3] -= 2 * M_PI;
                    }
                    configE += deltaE;
                }
                else
                {
                    for (int k=0; k<4; k++)
                    {
                        spinConfig[temp_flipPos][k] = temp_selectedSiteSpin[k];
                    }
                }



                // for (int k=0; k<m_neighborNum; k++)
                // {
                //     deltaE += -2 * m_interactionJ1 * spinConfig[temp_flipPos][0] * spinConfig[neighbor[temp_flipPos][k]][0];
                // }
                // deltaE -= -2 * m_externalMagField * spinConfig[temp_flipPos][0];

                // if (realRand(engineTime) < exp(-deltaE / m_temperature))
                // {
                //     spinConfig[temp_flipPos][0] = -spinConfig[temp_flipPos][0];
                //     configE += deltaE;
                // }
            }
        }
    }
}

void Lattice::GetBinaryConfigIsingEMP(int binaryCode, double result[4])
// GetSuchP
// Binarycode implementation
// \tau_lastmember,s_lastmember,\tau,s,...\tau_1,s_1
{
    double binary_configE=0;

    // for (int i=0; i<m_size1*m_size2; i++)
    // {
    //     spinConfig[i][0] = 2 * (binaryCode >> i & 1) - 1;
    // }

    for (int i=0; i<m_size1*m_size2; i++)
    {
        if ((binaryCode>>(2*i))&1)
        {
            spinConfig[i][0] = m_deducedIsingDirection[0];
            spinConfig[i][1] = m_deducedIsingDirection[1];
        }
        else
        {
            spinConfig[i][0] = M_PI - m_deducedIsingDirection[0];
            spinConfig[i][1] = M_PI + m_deducedIsingDirection[1];
        }
        if ((binaryCode>>(2*i+1))&1)
        {
            spinConfig[i][2] = m_deducedIsingDirection[2];
            spinConfig[i][3] = m_deducedIsingDirection[3];
        }
        else
        {
            spinConfig[i][2] = M_PI - m_deducedIsingDirection[2];
            spinConfig[i][3] = M_PI + m_deducedIsingDirection[3];
        }
    }

    for (int i=0; i<m_size1; i++)
    {
        for (int j=0; j<m_size2; j++)
        {
            // for (int k=0; k<m_neighborNum; k++)
            // {
            //     binary_configE += m_interactionJ1 * spinConfig[i*m_size1+j][0] * spinConfig[neighbor[i*m_size1+j][k]][0] / 2;
            // }
            // binary_configE -= m_externalMagField * spinConfig[i*m_size1+j][0];
            binary_configE += GetSiteHamiltonian(i*m_size1+j);
        }
    }

    result[0] = binary_configE;
    result[1] = GetMagnetization();
    result[2] = GetPolarizationX();
    result[3] = GetPolarizationY();
}

void CalculateTheoreticalParams(Lattice lattice, int tempNum, double temperatureSequence[TEMPERATURESEQUENCEMAX], double resultParam[8][TEMPERATURESEQUENCEMAX])
{
    double Z[tempNum] = {0.0};
    double HZ[tempNum] = {0.0};
    double HsqZ[tempNum] = {0.0};
    double MZ[tempNum] = {0.0};
    double MsqZ[tempNum] = {0.0};
    double PxZ[tempNum] = {0.0};
    double PxsqZ[tempNum] = {0.0};
    double PyZ[tempNum] = {0.0};
    double PysqZ[tempNum] = {0.0};

    double currentEMP[4];

    for (int binConfig=0; binConfig<(1<<(lattice.m_size1*lattice.m_size2*2)); binConfig++)
    {
        lattice.GetBinaryConfigIsingEMP(binConfig,currentEMP);
        for (int tempindex=0; tempindex<tempNum; tempindex++)
        {
            Z[tempindex] += exp(-currentEMP[0]/temperatureSequence[tempindex]);
            HZ[tempindex] += currentEMP[0] * exp(-currentEMP[0]/temperatureSequence[tempindex]);
            HsqZ[tempindex] += pow(currentEMP[0],2) * exp(-currentEMP[0]/temperatureSequence[tempindex]);
            MZ[tempindex] += currentEMP[1] * exp(-currentEMP[0]/temperatureSequence[tempindex]);
            MsqZ[tempindex] += pow(currentEMP[1],2) * exp(-currentEMP[0]/temperatureSequence[tempindex]);
            PxZ[tempindex] += currentEMP[2] * exp(-currentEMP[0]/temperatureSequence[tempindex]);
            PxsqZ[tempindex] += pow(currentEMP[2],2) * exp(-currentEMP[0]/temperatureSequence[tempindex]);
            PyZ[tempindex] += currentEMP[3] * exp(-currentEMP[0]/temperatureSequence[tempindex]);
            PysqZ[tempindex] += pow(currentEMP[3],2) * exp(-currentEMP[0]/temperatureSequence[tempindex]);
        }
    }
    for (int tempindex=0; tempindex<tempNum; tempindex++)   // check \epsilon equation
    {
        resultParam[0][tempindex] = HZ[tempindex] / Z[tempindex];
        resultParam[1][tempindex] = 1 / pow(temperatureSequence[tempindex],2) * (HsqZ[tempindex] * Z[tempindex] - pow(HZ[tempindex],2)) / pow(Z[tempindex],2);
        resultParam[2][tempindex] = MZ[tempindex] / Z[tempindex];
        resultParam[3][tempindex] = 1 / temperatureSequence[tempindex] * (MsqZ[tempindex] * Z[tempindex] - pow(MZ[tempindex],2)) / pow(Z[tempindex],2);
        resultParam[4][tempindex] = PxZ[tempindex] / Z[tempindex];
        resultParam[5][tempindex] = 1 / temperatureSequence[tempindex] * (PxsqZ[tempindex] * Z[tempindex] - pow(PxZ[tempindex],2)) / pow(Z[tempindex],2);
        resultParam[6][tempindex] = PyZ[tempindex] / Z[tempindex];
        resultParam[7][tempindex] = 1 / temperatureSequence[tempindex] * (PysqZ[tempindex] * Z[tempindex] - pow(PyZ[tempindex],2)) / pow(Z[tempindex],2);
    }
}

int main()
{
    int latSizeXMax = 2;
    int latSizeYMax = 2;
    string latType = "Triangle";
    double interactionJprime = -1.0;
    double interactionJprimeprime = 0.0;
    double externalXElectricField = 0;
    double externalYElectricField = 0;
    double externalZMageticField = 0;
    string interactionType = "Ising";
    double assignedDirection[4] = {M_PI/3,M_PI/3,M_PI/7,M_PI/7};    // Spherical Coordination

    double Tmin = 0.1;
    double Tmax = 5;
    double Tstep = 0.1;
    // double Tstep = 1;
    double temperature[TEMPERATURESEQUENCEMAX];

    int tempLength = (int)((Tmax-Tmin)/Tstep);
    for (int i=0; i<tempLength; i++)
    {
        temperature[i] = Tmax - Tstep * i;
    }

    // int simulationCount = 10;
    int simulationCount = 2;
    int totalSweep = 1E5;
    // int sweepPerSample = totalSweep / 50;
    int sweepPerSample = 50;
    double relaxationTime = 0.1;


    ofstream fout;
    clock_t start_time = clock();

    double resultE[simulationCount][tempLength];
    double resultEsq[simulationCount][tempLength];
    double resultM[simulationCount][tempLength];
    double resultMsq[simulationCount][tempLength];
    double resultPx[simulationCount][tempLength];
    double resultPxsq[simulationCount][tempLength];
    double resultPy[simulationCount][tempLength];
    double resultPysq[simulationCount][tempLength];

    double configEnergy;
    double sampledEnergy;
    double sampledEsq;
    double configMagnetization;
    double sampledM;
    double sampledMsq;
    double configPolarization[2];
    double sampledPx;
    double sampledPxsq;
    double sampledPy;
    double sampledPysq;

    int totalSampleNum;
    double theoreticalResult[8][TEMPERATURESEQUENCEMAX];

    // double energySequence[simulationCount][(int)(totalSweep*relaxationTime)];

    Lattice lattice(latSizeXMax,latSizeYMax,latType,interactionJprime,interactionJprimeprime,externalXElectricField,externalYElectricField,
        externalZMageticField,interactionType,assignedDirection);
    // double ans[4];
    // lattice.GetBinaryConfigIsingEMP(11,ans);
    // cout << ans[0] << " " << ans[1] << " " << ans[2] << " " << ans[3] << " " << endl;


    if (interactionType == "Ising")
    {
        CalculateTheoreticalParams(lattice, tempLength, temperature, theoreticalResult);
    }

    for (int simuindex=0; simuindex<simulationCount; simuindex++)
    {
        for (int tempindex=0; tempindex<tempLength; tempindex++)
        {
            configEnergy = 0;
            sampledEnergy = 0;
            sampledEsq = 0;
            configMagnetization = 0;
            sampledM = 0;
            sampledMsq = 0;
            configPolarization[0] = 0;
            configPolarization[1] = 0;
            sampledPx = 0;
            sampledPxsq = 0;
            sampledPy = 0;
            sampledPysq = 0;
            totalSampleNum = 0;

            // Lattice lattice(latSizeXMax,latSizeYMax,latType,interactionJ,temperature[tempindex],externalZMageticField,interactionType);
            lattice.SetTemperature(temperature[tempindex]);

            for (int i=0; i<totalSweep*relaxationTime; i++)
            {
                // configEnergy = lattice.GetEnergy();
                // energySequence[simuindex][i] = configEnergy;
                lattice.SweepFlip();
            }
            // lattice.DebugOutput();
            // lattice.SweepFlip();
            // lattice.DebugOutput();

            for (int i=0; i<totalSweep*(1-relaxationTime); i++)
            {
                lattice.SweepFlip();
                if (i % sweepPerSample == 0)
                {
                    totalSampleNum++;
                    configEnergy = lattice.GetEnergy();
                    sampledEnergy += configEnergy;
                    sampledEsq += pow(configEnergy,2);
                    configMagnetization = lattice.GetMagnetization();
                    sampledM += configMagnetization;
                    sampledMsq += pow(configMagnetization,2);
                    configPolarization[0] = lattice.GetPolarizationX();
                    sampledPx += configPolarization[0];
                    sampledPxsq += pow(configPolarization[0],2);
                    configPolarization[1] = lattice.GetPolarizationY();
                    sampledPy += configPolarization[1];
                    sampledPysq += pow(configPolarization[1],2);
                }
            }
            // cout << sampledEnergy << " " << sampledM << " " << sampledPx << " " << sampledPy << " " << endl;

            resultE[simuindex][tempindex] = sampledEnergy / totalSampleNum;
            resultEsq[simuindex][tempindex] = sampledEsq / totalSampleNum;
            resultM[simuindex][tempindex] = sampledM / totalSampleNum;
            resultMsq[simuindex][tempindex] = sampledMsq / totalSampleNum;
            resultPx[simuindex][tempindex] = sampledPx / totalSampleNum;
            resultPxsq[simuindex][tempindex] = sampledPxsq / totalSampleNum;
            resultPy[simuindex][tempindex] = sampledPy / totalSampleNum;
            resultPysq[simuindex][tempindex] = sampledPysq / totalSampleNum;

            cout << simuindex << " CoreNum " << tempindex << " temperature finished" << endl; 
        }
    }


    clock_t calc_time = clock();
    cout << "Calculation finished, elapsed time: " << (double)(calc_time - start_time) / CLOCKS_PER_SEC << "s" << endl;

    // fout.open("validateEquilibrium.csv",ios::out);
    // for (int i=0; i<simulationCount; i++)
    // {
    //     for (int j=0; j<totalSweep*relaxationTime; j++)
    //     {
    //         fout << energySequence[i][j] << ",";
    //     }
    //     fout << endl;
    // }
    // fout.close();

    fout.open("temperature.csv",ios::out);
    for (int i=0; i<tempLength; i++)
    {
        fout << temperature[i] << ",";
    }
    fout.close();

    fout.open("resultE.csv",ios::out);
    for (int i=0; i<simulationCount; i++)
    {
        for (int j=0; j<tempLength; j++)
        {
            fout << resultE[i][j] << ",";
        }
        fout << endl;
    }
    fout.close();

    fout.open("resultEsq.csv",ios::out);
    for (int i=0; i<simulationCount; i++)
    {
        for (int j=0; j<tempLength; j++)
        {
            fout << resultEsq[i][j] << ",";
        }
        fout << endl;
    }
    fout.close();

    fout.open("resultM.csv",ios::out);
    for (int i=0; i<simulationCount; i++)
    {
        for (int j=0; j<tempLength; j++)
        {
            fout << resultM[i][j] << ",";
        }
        fout << endl;
    }
    fout.close();

    fout.open("resultMsq.csv",ios::out);
    for (int i=0; i<simulationCount; i++)
    {
        for (int j=0; j<tempLength; j++)
        {
            fout << resultMsq[i][j] << ",";
        }
        fout << endl;
    }
    fout.close();

    fout.open("resultPx.csv",ios::out);
    for (int i=0; i<simulationCount; i++)
    {
        for (int j=0; j<tempLength; j++)
        {
            fout << resultPx[i][j] << ",";
        }
        fout << endl;
    }
    fout.close();

    fout.open("resultPxsq.csv",ios::out);
    for (int i=0; i<simulationCount; i++)
    {
        for (int j=0; j<tempLength; j++)
        {
            fout << resultPxsq[i][j] << ",";
        }
        fout << endl;
    }
    fout.close();

    fout.open("resultPy.csv",ios::out);
    for (int i=0; i<simulationCount; i++)
    {
        for (int j=0; j<tempLength; j++)
        {
            fout << resultPy[i][j] << ",";
        }
        fout << endl;
    }
    fout.close();

    fout.open("resultPysq.csv",ios::out);
    for (int i=0; i<simulationCount; i++)
    {
        for (int j=0; j<tempLength; j++)
        {
            fout << resultPysq[i][j] << ",";
        }
        fout << endl;
    }
    fout.close();

    fout.open("theoretical.csv",ios::out);
    for (int i=0; i<8; i++)
    {
        for (int j=0; j<tempLength; j++)
        {
            fout << theoreticalResult[i][j] << ",";
        }
        fout << endl;
    }
    fout.close();
}