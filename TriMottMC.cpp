#include <fstream>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <random>
#include <cmath>
using namespace std;

std::uniform_real_distribution<double> realRand(0,1);
std::uniform_int_distribution<int> intRand(0,65535);
std::default_random_engine engineTime(time(NULL));

class Lattice
{
    static const int LATTICESIZE_MAX = 100;
    static const int NEIGHBOR_MAX = 6;

    public:
        Lattice(int row, int column, int neighborNum, double interactionJ, double temperature, double externalMagField);

        void SetNeighborIndex();
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
            for (int i=0; i<m_size1*m_size2; i++)
            {
                configM += GetSiteZComponent(i);
            }
            return configM;
        }

        void DebugOutput()
        {
            ShowConfig();
            cout << GetSiteInnerProduct(0,1) << endl;
        }

        ~Lattice()
        {
            delete temp_stackedIndex;
            delete temp_randomSphericalDirection;
            delete temp_selectedSiteSpin;
        }

    private:
        int m_size1;
        int m_size2;
        int m_neighborNum;
        double m_interactionJ;
        double m_temperature;
        double m_externalMagField;

        double configE;
        double deltaE;
        double configM;

        double spinConfig[LATTICESIZE_MAX * LATTICESIZE_MAX][2];
        int neighbor[LATTICESIZE_MAX * LATTICESIZE_MAX][NEIGHBOR_MAX];
        int *temp_stackedIndex = new int[2];
        double *temp_randomSphericalDirection = new double[2];
        int temp_flipPos;
        double *temp_selectedSiteSpin = new double[2];

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

Lattice::Lattice(int row, int column, int neighborNum, double interactionJ, double temperature, double externalMagField=0)
    : m_size1(column), m_size2(row), m_neighborNum(neighborNum), m_interactionJ(interactionJ), m_temperature(temperature), m_externalMagField(externalMagField)
{
    SetNeighborIndex();
    InitSpin();
    configE = 0;
    for (int i=0; i<m_size1; i++)
    {
        for (int j=0; j<m_size2; j++)
        {
            for (int k=0; k<m_neighborNum; k++)
            {
                // configE += m_interactionJ * spinConfig[i*m_size+j] * spinConfig[neighbor[i*m_size+j][k]];
                configE += m_interactionJ * Lattice::GetSiteInnerProduct(i*m_size1+j,neighbor[i*m_size1+j][k]) / 2;
            }
            configE -= m_externalMagField * Lattice::GetSiteZComponent(i*m_size1+j);
        }
    }
}

void Lattice::SetNeighborIndex()
{
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

void Lattice::InitSpin()
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
}

void Lattice::ShowConfig()
{
    cout << "Config:" << endl;
    for (int i=0; i<m_size1; i++)
    {
        for (int j=0; j<m_size2; j++)
        {
            cout << "(" << spinConfig[i*m_size1+j][0] << "," << spinConfig[i*m_size1+j][1] << ")|";
        }
        cout << endl;
    }
}

void Lattice::SweepFlip()
{
    for (int i=0; i<m_size1; i++)
    {
        for (int j=0; j<m_size2; j++)
        {
            deltaE = 0;
            // int flipPos = intRand(engineTime) % (int) pow(m_size,2);
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
    }
}


int main()
{
    int latSizeRow = 1;
    int latSizeColumn = 2;
    int latNeighbor = 4;
    double interactionJ = -1.0;
    double externalZMageticField = 0;

    double Tmin = 0.1;
    double Tmax = 5;
    double Tstep = 0.1;
    double temperature[(int)((Tmax-Tmin)/Tstep)];

    int tempLength = sizeof(temperature) / sizeof(temperature[0]);
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
    // double energySequence[simulationCount][(int)(totalSweep*relaxationTime)];

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

            Lattice lattice(latSizeRow,latSizeColumn,latNeighbor,interactionJ,temperature[tempindex],externalZMageticField);

            for (int i=0; i<totalSweep*relaxationTime; i++)
            {
                // configEnergy = lattice.GetEnergy();
                // energySequence[simuindex][i] = configEnergy;
                lattice.SweepFlip();
            }

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
}