#include <cstdlib>
#include <cstring>
#include <cmath>
#include <random>
#include "fftw3.h"


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
        double externalMagField, string interactionType, double* DeducedIsingDirection);

        int m_totalSite;

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
            //     for (int i=0; i<m_totalSite; i++)
            //     {
            //         configM += GetSiteZSpinComponent(i);
            //     }
            // }
            // else if (m_interactionType == "Ising")
            // {
            //     for (int i=0; i<m_totalSite; i++)
            //     {
            //         configM += spinConfig[i][0];
            //     }
            // }
            for (int i=0; i<m_totalSite; i++)
            {
                configM += GetSiteZSpinComponent(i);
            }
            return configM;
        }

        double GetPolarizationX()
        {
            configPx = 0;
            for (int i=0; i<m_totalSite; i++)
            {
                configPx += GetSiteXPseudospinComponent(i);
            }
            return configPx;
        }

        double GetPolarizationY()
        {
            configPy = 0;
            for (int i=0; i<m_totalSite; i++)
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
            for (int i=0; i<m_totalSite; i++)
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
            
            double norm[m_totalSite];

            GetFFTnorm("pseudospin","x",norm);
            
            for (int i=0; i<m_totalSite; i++)
            {
                cout << in[i][0] << " " << in[i][1] << endl;
            }
            for (int i=0; i<m_totalSite; i++)
            {
                cout << norm[i] << " ";
            }
            cout << endl;
        }

        void DebugOutput()
        {
            ShowConfig();
            cout << GetEnergy() << endl;
            // cout << GetSiteInnerProduct(0,1) << endl;
        }

        void GetBinaryConfigIsingEMP(int binaryCode, double result[2]);

        void GetFFTnorm(string component, string direction, double* norm);

        ~Lattice()
        {
            fftw_destroy_plan(p);
            fftw_free(in); fftw_free(out);
        }

    private:
        int m_size1;
        int m_size2;
        int m_neighborNum;
        double m_interactionJ1;
        double m_interactionJ2;
        double m_temperature;
        double m_externalMagField;
        double m_externalXElecField;
        double m_externalYElecField;
        string m_interactionType;
        double* m_deducedIsingDirection;
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

        fftw_complex *in, *out;
        fftw_plan p;

        int GetFlattenedCoordinate(int stackedIndex1, int stackedIndex2)
        {
            if (stackedIndex1 < 0)  stackedIndex1 += m_size1;
            if (stackedIndex1 >= m_size1)  stackedIndex1 -= m_size1;
            if (stackedIndex2 < 0)  stackedIndex2 += m_size2;
            if (stackedIndex2 >= m_size2)  stackedIndex2 -= m_size2;

            return stackedIndex1 * m_size2 + stackedIndex2;
        }

        void GetStackedCoordinate(int flattenedIndex)
        {
            temp_stackedIndex[0] = int(flattenedIndex / m_size2);
            temp_stackedIndex[1] = flattenedIndex % m_size2;
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
            temp_siteE -= m_externalMagField * GetSiteZSpinComponent(site); // g\mu_B = 1
            return temp_siteE;
        }

        double GetSiteZSpinComponent(int site)
        {
            return cos(spinConfig[site][0]) * m_s;
        }

        double GetSiteXPseudospinComponent(int site)
        {
            return sin(spinConfig[site][2]) * cos(spinConfig[site][3]) * m_tau;
        }

        double GetSiteYPseudospinComponent(int site)
        {
            return sin(spinConfig[site][2]) * sin(spinConfig[site][3]) * m_tau;
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

Lattice::Lattice(int maxX, int maxY, string latType, double interactionJ1, double interactionJ2, double externalXElecField=0, double externalYElecField=0,
     double externalMagField=0, string interactionType="Ising", double* deducedIsingDirection=nullptr)
    : m_size1(maxX), m_size2(maxY), m_interactionJ1(interactionJ1), m_interactionJ2(interactionJ2), m_externalMagField(externalMagField),
        m_externalXElecField(externalXElecField), m_externalYElecField(externalYElecField), m_interactionType(interactionType), m_deducedIsingDirection(deducedIsingDirection)
{
    if (m_deducedIsingDirection == nullptr)
    {
        double temp_Zaxis[4] = {0.0};
        m_deducedIsingDirection = temp_Zaxis;
    }
    m_totalSite = m_size1 * m_size2;
    SetNeighborIndex(latType);
    InitSpin();

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_totalSite);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * m_totalSite);
    p = fftw_plan_dft_2d(m_size1, m_size2, in, out, FFTW_FORWARD, FFTW_MEASURE);
}

void Lattice::SetNeighborIndex(string latType)  // i, j aligned to x, y. For further use and adding z, add to the 3rd index and add to 3rd loop layer
{
    if (latType == "Square")
    {
        m_neighborNum = 4;
        for (int i=0; i<m_size1; i++)
        {
            for (int j=0; j<m_size2; j++)
            {
                neighbor[i*m_size2+j][0] = Lattice::GetFlattenedCoordinate(i,j+1);
                neighbor[i*m_size2+j][1] = Lattice::GetFlattenedCoordinate(i-1,j);
                neighbor[i*m_size2+j][2] = Lattice::GetFlattenedCoordinate(i,j-1);
                neighbor[i*m_size2+j][3] = Lattice::GetFlattenedCoordinate(i+1,j);
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
                neighbor[i*m_size2+j][0] = Lattice::GetFlattenedCoordinate(i,j+1);
                neighbor[i*m_size2+j][1] = Lattice::GetFlattenedCoordinate(i,j-1);
                if (i % 2 == 0)
                {
                    neighbor[i*m_size2+j][2] = Lattice::GetFlattenedCoordinate(i+1,j+1);
                }
                else
                {
                    neighbor[i*m_size2+j][2] = Lattice::GetFlattenedCoordinate(i-1,j-1);
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
                neighbor[i*m_size2+j][0] = Lattice::GetFlattenedCoordinate(i+1,j);
                neighbor[i*m_size2+j][1] = Lattice::GetFlattenedCoordinate(i,j-1);
                neighbor[i*m_size2+j][2] = Lattice::GetFlattenedCoordinate(i-1,j-1);
                neighbor[i*m_size2+j][3] = Lattice::GetFlattenedCoordinate(i-1,j);
                neighbor[i*m_size2+j][4] = Lattice::GetFlattenedCoordinate(i,j+1);
                neighbor[i*m_size2+j][5] = Lattice::GetFlattenedCoordinate(i+1,j+1);
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
                spinConfig[i*m_size2+j][0] = temp_randomSphericalDirection[0];
                spinConfig[i*m_size2+j][1] = temp_randomSphericalDirection[1];
                Lattice::SetRandomSphericalDirection();
                spinConfig[i*m_size2+j][2] = temp_randomSphericalDirection[0];
                spinConfig[i*m_size2+j][3] = temp_randomSphericalDirection[1];
            }
        }

        for (int i=0; i<m_size1; i++)
        {
            for (int j=0; j<m_size2; j++)
            {
                // for (int k=0; k<m_neighborNum; k++)
                // {
                //     configE += m_interactionJ1 * Lattice::GetSiteInnerProduct(i*m_size2+j,neighbor[i*m_size2+j][k]) / 2;
                // }
                // configE -= m_externalMagField * Lattice::GetSiteZSpinComponent(i*m_size2+j);
                configE += GetSiteHamiltonian(i*m_size2+j);
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
                    spinConfig[i*m_size2+j][0] = m_deducedIsingDirection[0];
                    spinConfig[i*m_size2+j][1] = m_deducedIsingDirection[1];
                }
                else
                {
                    spinConfig[i*m_size2+j][0] = M_PI - m_deducedIsingDirection[0];
                    spinConfig[i*m_size2+j][1] = M_PI + m_deducedIsingDirection[1];
                }
                if (temp_orientation>>1&1)
                {
                    spinConfig[i*m_size2+j][2] = m_deducedIsingDirection[2];
                    spinConfig[i*m_size2+j][3] = m_deducedIsingDirection[3];
                }
                else
                {
                    spinConfig[i*m_size2+j][2] = M_PI - m_deducedIsingDirection[2];
                    spinConfig[i*m_size2+j][3] = M_PI + m_deducedIsingDirection[3];
                }
            }
        }

        for (int i=0; i<m_size1; i++)
        {
            for (int j=0; j<m_size2; j++)
            {
                // for (int k=0; k<m_neighborNum; k++)
                // {
                //     configE += m_interactionJ1 * spinConfig[i*m_size2+j][0] * spinConfig[neighbor[i*m_size2+j][k]][0] / 2;
                // }
                // configE -= m_externalMagField * spinConfig[i*m_size2+j][0];
                configE += GetSiteHamiltonian(i*m_size2+j);
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
            cout << "(" << spinConfig[i*m_size2+j][0] << "," << spinConfig[i*m_size2+j][1] << ") " << "(" << spinConfig[i*m_size2+j][2] << "," << spinConfig[i*m_size2+j][3] << ")|";
        }
        cout << endl;
    }
    for (int i=0; i<m_neighborNum; i++)
    {
        cout << neighbor[0][i] << endl;
    }
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
                // temp_flipPos = i * m_size2 + j;
                // temp_selectedSiteSpin[0] = spinConfig[temp_flipPos][0];
                // temp_selectedSiteSpin[1] = spinConfig[temp_flipPos][1];
                temp_flipPos = intRand(engineTime) % ((m_totalSite)*2);
                isComponentSFlipped = temp_flipPos % 2;
                temp_flipPos /= 2;

                for (int k=0; k<4; k++)
                {
                    temp_selectedSiteSpin[k] = spinConfig[temp_flipPos][k];
                }

                Lattice::SetRandomSphericalDirection();

                deltaE -= GetSiteHamiltonian(temp_flipPos);
                for (int k=0; k<3; k++)  // remember to modify when add layers
                {
                    deltaE -= GetSiteHamiltonian(neighbor[temp_flipPos][2*k+1]);
                }

                if (isComponentSFlipped)
                {
                    spinConfig[temp_flipPos][0] = temp_randomSphericalDirection[0];
                    spinConfig[temp_flipPos][1] = temp_randomSphericalDirection[1];
                }
                else
                {
                    spinConfig[temp_flipPos][2] = temp_randomSphericalDirection[0];
                    spinConfig[temp_flipPos][3] = temp_randomSphericalDirection[1]; 
                }

                deltaE += GetSiteHamiltonian(temp_flipPos);
                for (int k=0; k<3; k++)  // remember to modify when add layers
                {
                    deltaE += GetSiteHamiltonian(neighbor[temp_flipPos][2*k+1]);
                }

                if (realRand(engineTime) < exp(-deltaE / m_temperature))
                {
                    // spinConfig[flipPos] = -spinConfig[flipPos];
                    configE += deltaE;
                }
                else
                {
                    for (int k=0; k<4; k++)
                    {
                        spinConfig[temp_flipPos][k] = temp_selectedSiteSpin[k];
                    }
                }
            }

            else if (m_interactionType == "Ising")
            {
                temp_flipPos = intRand(engineTime) % ((m_totalSite)*2);
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
// \tau_lastmember,s_lastmember,\tau,s,...\tau_1,s_1
{
    double binary_configE=0;

    // for (int i=0; i<m_totalSite; i++)
    // {
    //     spinConfig[i][0] = 2 * (binaryCode >> i & 1) - 1;
    // }

    for (int i=0; i<m_totalSite; i++)
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
            //     binary_configE += m_interactionJ1 * spinConfig[i*m_size2+j][0] * spinConfig[neighbor[i*m_size2+j][k]][0] / 2;
            // }
            // binary_configE -= m_externalMagField * spinConfig[i*m_size2+j][0];
            binary_configE += GetSiteHamiltonian(i*m_size2+j);
        }
    }

    result[0] = binary_configE;
    result[1] = GetMagnetization();
    result[2] = GetPolarizationX();
    result[3] = GetPolarizationY();
}

void Lattice::GetFFTnorm(string component, string direction, double* norm) // don't forget the scaling factor
{
    if (component == "spin")
    {
        if (direction == "z")
        {
            for (int i=0; i<m_totalSite; i++)
            {
                in[i][0] = GetSiteZSpinComponent(i);
                in[i][1] = 0;
            }
        }
    }

    if (component == "pseudospin")
    {
        if (direction == "x")
        {
            for (int i=0; i<m_totalSite; i++)
            {
                in[i][0] = GetSiteXPseudospinComponent(i);
                in[i][1] = 0;
            }
        }
        else if (direction == "y")
        {
            for (int i=0; i<m_totalSite; i++)
            {
                in[i][0] = GetSiteYPseudospinComponent(i);
                in[i][1] = 0;
            }
        }
    }

    fftw_execute(p);

    for (int i=0; i<m_totalSite; i++)
    {
        norm[i] = (pow(out[i][0],2) + pow(out[i][1],2)) / pow(m_totalSite,2);
    }
}
