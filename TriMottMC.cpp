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
        Lattice(int size, int neighborNum, double interactionJ, double temperature);

        void SetNeighborIndex();
        void InitSpin();
        void ShowConfig();
        void SweepFlip();

        double GetEnergy()
        {
            return configE;
        }

        double GetAveMagnetization()
        {
            configM = 0;
            for (int i=0; i<pow(m_size,2); i++)
            {
                configM += spinConfig[i];
            }
            return configM / pow(m_size,2);
        }

        void DebugOutput()
        {
            cout << spinConfig[0] << endl << configE;
        }

        ~Lattice()
        {
            delete temp_stackedIndex;
        }

    private:
        int m_size;
        int m_neighborNum;
        double m_interactionJ;
        double m_temperature;

        double configE;
        double deltaE;
        double configM;

        int spinConfig[LATTICESIZE_MAX * LATTICESIZE_MAX];
        int neighbor[LATTICESIZE_MAX * LATTICESIZE_MAX][NEIGHBOR_MAX];
        int *temp_stackedIndex = new int[2];

        int GetFlattenedCoordinate(int stackedIndex1, int stackedIndex2)
        {
            if (stackedIndex1 < 0)  stackedIndex1 += m_size;
            if (stackedIndex1 >= m_size)  stackedIndex1 -= m_size;
            if (stackedIndex2 < 0)  stackedIndex2 += m_size;
            if (stackedIndex2 >= m_size)  stackedIndex2 -= m_size;

            return stackedIndex1 * m_size + stackedIndex2;
        }

        void GetStackedCoordinate(int flattenedIndex)
        {
            temp_stackedIndex[0] = int(flattenedIndex / m_size);
            temp_stackedIndex[1] = flattenedIndex % m_size;
        }
};

Lattice::Lattice(int size, int neighborNum, double interactionJ, double temperature)
    : m_size(size), m_neighborNum(neighborNum), m_interactionJ(interactionJ), m_temperature(temperature)
{
    SetNeighborIndex();
    InitSpin();
    configE = 0;
    for (int i=0; i<m_size; i++)
    {
        for (int j=0; j<m_size; j++)
        {
            for (int k=0; k<m_neighborNum; k++)
            {
                configE += m_interactionJ * spinConfig[i*m_size+j] * spinConfig[neighbor[i*m_size+j][k]];
            }
        }
    }
    configE /= 2;
}

void Lattice::SetNeighborIndex()
{
    for (int i=0; i<m_size; i++)
    {
        for (int j=0; j<m_size; j++)
        {
            neighbor[i*m_size+j][0] = Lattice::GetFlattenedCoordinate(i,j+1);
            neighbor[i*m_size+j][1] = Lattice::GetFlattenedCoordinate(i-1,j);
            neighbor[i*m_size+j][2] = Lattice::GetFlattenedCoordinate(i,j-1);
            neighbor[i*m_size+j][3] = Lattice::GetFlattenedCoordinate(i+1,j);
        }
    }
}

void Lattice::InitSpin()
{
    for (int i=0; i<m_size; i++)
    {
        for (int j=0; j<m_size; j++)
        {
            spinConfig[i*m_size+j] = (intRand(engineTime) % 2) * 2 - 1;
        }
    }
}

void Lattice::ShowConfig()
{
    cout << "Config:" << endl;
    for (int i=0; i<m_size; i++)
    {
        for (int j=0; j<m_size; j++)
        {
            cout << spinConfig[i*m_size+j] << " ";
        }
        cout << endl;
    }
}

void Lattice::SweepFlip()
{
    for (int i=0; i<m_size; i++)
    {
        for (int j=0; j<m_size; j++)
        {
            deltaE = 0;
            int flipPos = intRand(engineTime) % (int) pow(m_size,2);
            for (int k=0; k<m_neighborNum; k++)
            {
                deltaE += -2 * m_interactionJ * spinConfig[flipPos] * spinConfig[neighbor[flipPos][k]];
            }


            if (realRand(engineTime) < exp(-deltaE / m_temperature))
            {
                spinConfig[flipPos] = -spinConfig[flipPos];
                configE += deltaE;
            }
        }
    }
}


int main()
{
    int latSize = 2;
    int latNeighbor = 4;
    double interactionJ = -1.0;

    double Tmin = 0.5;
    double Tmax = 5;
    double Tstep = 0.05;
    double temperature[(int)((Tmax-Tmin)/Tstep)];

    int tempLength = sizeof(temperature) / sizeof(temperature[0]);
    for (int i=0; i<tempLength; i++)
    {
        temperature[i] = Tmin + Tstep * i;
    }

    int simulationCount = 10;
    int totalSweep = 1E5;
    int sweepPerSample = totalSweep / 50;
    double relaxationTime = 0.1;


    ofstream fout;
    clock_t start_time = clock();

    double resultE[simulationCount][tempLength];
    double resultEsq[simulationCount][tempLength];
    double resultM[simulationCount][tempLength];
    double resultMsq[simulationCount][tempLength];

    for (int simuindex=0; simuindex<simulationCount; simuindex++)
    {
        for (int tempindex=0; tempindex<tempLength; tempindex++)
        {
            double configEnergy = 0;
            double sampledEnergy = 0;
            double sampledEsq = 0;
            double configMagnetization = 0;
            double sampledM = 0;
            double sampledMsq = 0;
            int totalSampleNum = 0;

            Lattice lattice(latSize,latNeighbor,interactionJ,temperature[tempindex]);

            for (int i=0; i<totalSweep*relaxationTime; i++)
            {
                lattice.SweepFlip();
            }

            for (int i=0; i<totalSweep*(1-relaxationTime); i++)
            {
                lattice.SweepFlip();
                if (i % sweepPerSample == 0)
                {
                    totalSampleNum++;
                    configEnergy = lattice.GetEnergy();
                    sampledEnergy += configEnergy;
                    sampledEsq += pow(configEnergy,2);
                    configMagnetization = lattice.GetAveMagnetization();
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