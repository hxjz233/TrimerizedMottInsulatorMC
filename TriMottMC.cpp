#include <fstream>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <random>
#include <cmath>
#define TEMPERATURESEQUENCEMAX 1000
using namespace std;

std::uniform_real_distribution<double> realRand(0,1);
std::uniform_int_distribution<int> intRand(0,65535);
std::default_random_engine engineTime(time(NULL));

class Lattice
{
    static const int LATTICESIZE_MAX = 100;
    static const int NEIGHBOR_MAX = 6;

    public:
        Lattice(int row, int column, string latType, double interactionJ, double temperature, double externalMagField, string interactionType);

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
            if (m_interactionType == "Heisenberg")
            {
                for (int i=0; i<m_size1*m_size2; i++)
                {
                    configM += GetSiteZComponent(i);
                }
            }
            else if (m_interactionType == "Ising")
            {
                for (int i=0; i<m_size1*m_size2; i++)
                {
                    configM += spinConfig[i][0];
                }
            }
            return configM;
        }

        void DebugOutput()
        {
            ShowConfig();
            // cout << GetSiteInnerProduct(0,1) << endl;
        }

        void GetBinaryConfigIsingEM(int binaryCode, double result[2]);

        // ~Lattice()
        // {
        //     delete temp_stackedIndex;
        //     delete temp_randomSphericalDirection;
        //     delete temp_selectedSiteSpin;
        // }

    private:
        int m_neighborNum;
        double m_interactionJ;
        double m_temperature;
        double m_externalMagField;
        string m_interactionType;

        double configE;
        double deltaE;
        double configM;

        double spinConfig[LATTICESIZE_MAX * LATTICESIZE_MAX][2];
        int neighbor[LATTICESIZE_MAX * LATTICESIZE_MAX][NEIGHBOR_MAX];
        // int *temp_stackedIndex = new int[2];
        // double *temp_randomSphericalDirection = new double[2];
        // int temp_flipPos;
        // double *temp_selectedSiteSpin = new double[2];
        int temp_stackedIndex[2];
        double temp_randomSphericalDirection[2];
        int temp_flipPos;
        double temp_selectedSiteSpin[2];


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

        double GetSiteZComponent(int site)
        {
            return cos(spinConfig[site][0]);
        }

        double GetSiteInnerProduct(int site1, int site2)
        {
            return sin(spinConfig[site1][0]) * sin(spinConfig[site2][0]) * cos(spinConfig[site1][1]-spinConfig[site2][1])
                 + cos(spinConfig[site1][0]) * cos(spinConfig[site2][0]);
        }
};

Lattice::Lattice(int row, int column, string latType, double interactionJ, double temperature, double externalMagField=0, string interactionType="Ising")
    : m_size1(column), m_size2(row), m_interactionJ(interactionJ), m_temperature(temperature), m_externalMagField(externalMagField), m_interactionType(interactionType)
{
    SetNeighborIndex(latType);
    InitSpin();
}

void Lattice::SetNeighborIndex(string latType)
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
                neighbor[i*m_size1+j][0] = Lattice::GetFlattenedCoordinate(i,j+1);
                neighbor[i*m_size1+j][1] = Lattice::GetFlattenedCoordinate(i-1,j);
                neighbor[i*m_size1+j][2] = Lattice::GetFlattenedCoordinate(i-1,j-1);
                neighbor[i*m_size1+j][3] = Lattice::GetFlattenedCoordinate(i,j-1);
                neighbor[i*m_size1+j][4] = Lattice::GetFlattenedCoordinate(i+1,j);
                neighbor[i*m_size1+j][5] = Lattice::GetFlattenedCoordinate(i+1,j+1);
            }
        }
    }
    // Throw exception
}

void Lattice::InitSpin()
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
                    configE += m_interactionJ * Lattice::GetSiteInnerProduct(i*m_size1+j,neighbor[i*m_size1+j][k]) / 2;
                }
                configE -= m_externalMagField * Lattice::GetSiteZComponent(i*m_size1+j);
            }
        }
    }

    else if (m_interactionType == "Ising")
    {
        for (int i=0; i<m_size1; i++)
        {
            for (int j=0; j<m_size2; j++)
            {
                spinConfig[i*m_size1+j][0] = (intRand(engineTime) % 2) * 2 - 1;
            }
        }

        for (int i=0; i<m_size1; i++)
        {
            for (int j=0; j<m_size2; j++)
            {
                for (int k=0; k<m_neighborNum; k++)
                {
                    configE += m_interactionJ * spinConfig[i*m_size1+j][0] * spinConfig[neighbor[i*m_size1+j][k]][0] / 2;
                }
                configE -= m_externalMagField * spinConfig[i*m_size1+j][0];
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

void Lattice::SweepFlip()
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
                    deltaE -= m_interactionJ * Lattice::GetSiteInnerProduct(temp_flipPos,neighbor[temp_flipPos][k]);
                }
                deltaE += m_externalMagField * Lattice::GetSiteZComponent(temp_flipPos);

                spinConfig[temp_flipPos][0] = temp_randomSphericalDirection[0];
                spinConfig[temp_flipPos][1] = temp_randomSphericalDirection[1];

                for (int k=0; k<m_neighborNum; k++)
                {
                    deltaE += m_interactionJ * Lattice::GetSiteInnerProduct(temp_flipPos,neighbor[temp_flipPos][k]);
                }
                deltaE -= m_externalMagField * Lattice::GetSiteZComponent(temp_flipPos);

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
                temp_flipPos = intRand(engineTime) % (m_size1*m_size2);
                for (int k=0; k<m_neighborNum; k++)
                {
                    deltaE += -2 * m_interactionJ * spinConfig[temp_flipPos][0] * spinConfig[neighbor[temp_flipPos][k]][0];
                }
                deltaE -= -2 * m_externalMagField * spinConfig[temp_flipPos][0];

                if (realRand(engineTime) < exp(-deltaE / m_temperature))
                {
                    spinConfig[temp_flipPos][0] = -spinConfig[temp_flipPos][0];
                    configE += deltaE;
                }
            }
        }
    }
}

void Lattice::GetBinaryConfigIsingEM(int binaryCode, double result[2])
{
    double binary_configE=0;

    for (int i=0; i<m_size1*m_size2; i++)
    {
        if (binaryCode>>i&1)
        {
            spinConfig[i][0] = 1;
        }
        else
        {
            spinConfig[i][0] = -1;
        }
    }

    for (int i=0; i<m_size1; i++)
    {
        for (int j=0; j<m_size2; j++)
        {
            for (int k=0; k<m_neighborNum; k++)
            {
                binary_configE += m_interactionJ * spinConfig[i*m_size1+j][0] * spinConfig[neighbor[i*m_size1+j][k]][0] / 2;
            }
            binary_configE -= m_externalMagField * spinConfig[i*m_size1+j][0];
        }
    }

    result[0] = binary_configE;
    result[1] = GetMagnetization();
}

void CalculateTheoreticalParams(Lattice lattice, int tempNum, double temperatureSequence[TEMPERATURESEQUENCEMAX], double resultParam[4][TEMPERATURESEQUENCEMAX])
{
    double Z[tempNum] = {0.0};
    double HZ[tempNum] = {0.0};
    double HsqZ[tempNum] = {0.0};
    double MZ[tempNum] = {0.0};
    double MsqZ[tempNum] = {0.0};

    double Energy_Magnetization[2];

    for (int binConfig=0; binConfig<(1<<lattice.m_size1*lattice.m_size2); binConfig++)
    {
        lattice.GetBinaryConfigIsingEM(binConfig,Energy_Magnetization);
        for (int tempindex=0; tempindex<tempNum; tempindex++)
        {
            Z[tempindex] += exp(-Energy_Magnetization[0]/temperatureSequence[tempindex]);
            HZ[tempindex] += Energy_Magnetization[0] * exp(-Energy_Magnetization[0]/temperatureSequence[tempindex]);
            HsqZ[tempindex] += pow(Energy_Magnetization[0],2) * exp(-Energy_Magnetization[0]/temperatureSequence[tempindex]);
            MZ[tempindex] += Energy_Magnetization[1] * exp(-Energy_Magnetization[0]/temperatureSequence[tempindex]);
            MsqZ[tempindex] += pow(Energy_Magnetization[1],2) * exp(-Energy_Magnetization[0]/temperatureSequence[tempindex]);
        }
    }
    for (int tempindex=0; tempindex<tempNum; tempindex++)
    {
        resultParam[0][tempindex] = HZ[tempindex] / Z[tempindex];
        resultParam[1][tempindex] = 1 / pow(temperatureSequence[tempindex],2) * (HsqZ[tempindex] * Z[tempindex] - pow(HZ[tempindex],2)) / pow(Z[tempindex],2);
        resultParam[2][tempindex] = MZ[tempindex] / Z[tempindex];
        resultParam[3][tempindex] = 1 / temperatureSequence[tempindex] * (MsqZ[tempindex] * Z[tempindex] - pow(MZ[tempindex],2)) / pow(Z[tempindex],2);
    }
}

int main()
{
    int latSizeXMax = 3;
    int latSizeYMax = 3;
    string latType = "Triangle";
    double interactionJ = -1.0;
    double externalZMageticField = 1;
    string interactionType = "Ising";

    double Tmin = 0.1;
    double Tmax = 5;
    double Tstep = 0.1;
    double temperature[TEMPERATURESEQUENCEMAX];

    int tempLength = (int)((Tmax-Tmin)/Tstep);
    for (int i=0; i<tempLength; i++)
    {
        temperature[i] = Tmin + Tstep * i;
    }

    int simulationCount = 10;
    // int simulationCount = 1;
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

    double configEnergy;
    double sampledEnergy;
    double sampledEsq;
    double configMagnetization;
    double sampledM;
    double sampledMsq;
    int totalSampleNum;
    double theoreticalResult[4][TEMPERATURESEQUENCEMAX];
    // double energySequence[simulationCount][(int)(totalSweep*relaxationTime)];

    // Lattice lattice(latSizeXMax,latSizeYMax,latType,interactionJ,1,externalZMageticField,interactionType);
    // double ans[2];
    // lattice.GetBinaryConfigIsingEM(15,ans);
    // cout << ans[0] << " " << ans[1] << endl;
    if (interactionType == "Ising")
    {
        Lattice lattice(latSizeXMax,latSizeYMax,latType,interactionJ,1,externalZMageticField,interactionType);
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
            totalSampleNum = 0;

            Lattice lattice(latSizeXMax,latSizeYMax,latType,interactionJ,temperature[tempindex],externalZMageticField,interactionType);

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
                    configEnergy = lattice.GetEnergy();
                    sampledEnergy += configEnergy;
                    sampledEsq += pow(configEnergy,2);
                    totalSampleNum++;
                    configMagnetization = lattice.GetMagnetization();
                    sampledM += configMagnetization;
                    sampledMsq += pow(configMagnetization,2);
                }
            }

            resultE[simuindex][tempindex] = sampledEnergy / totalSampleNum;
            resultEsq[simuindex][tempindex] = sampledEsq / totalSampleNum;
            resultM[simuindex][tempindex] = sampledM / totalSampleNum;
            resultMsq[simuindex][tempindex] = sampledMsq / totalSampleNum;
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

    fout.open("theoretical.csv",ios::out);
    for (int i=0; i<4; i++)
    {
        for (int j=0; j<tempLength; j++)
        {
            fout << theoreticalResult[i][j] << ",";
        }
        fout << endl;
    }
    fout.close();
}