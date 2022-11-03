#include <fstream>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <random>
#include <cmath>
#include <cstring>
// #include "fftw3.h"

#include "lattice.h"

#define TEMPERATURESEQUENCEMAX 1000
#define LATTICESIZEMAX 10000
using namespace std;

void CalculateTheoreticalParams(Lattice lattice, int tempNum, double temperatureSequence[TEMPERATURESEQUENCEMAX], 
    double resultParam[8][TEMPERATURESEQUENCEMAX], double resultNormParam[3][TEMPERATURESEQUENCEMAX][LATTICESIZEMAX])
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

    double MnormZ[tempNum][lattice.m_totalSite] = {0.0};
    double MnormsqZ[tempNum][lattice.m_totalSite] = {0.0};
    double PxnormZ[tempNum][lattice.m_totalSite] = {0.0};
    double PxnormsqZ[tempNum][lattice.m_totalSite] = {0.0};
    double PynormZ[tempNum][lattice.m_totalSite] = {0.0};
    double PynormsqZ[tempNum][lattice.m_totalSite] = {0.0};

    double currentEMP[4];
    double currentNorm[lattice.m_totalSite];

    for (int binConfig=0; binConfig<(1<<(lattice.m_totalSite*2)); binConfig++)
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

            lattice.GetFFTnorm("spin","z",currentNorm);
            for (int i=0; i<lattice.m_totalSite; i++)
            {
                MnormZ[tempindex][i] += currentNorm[i] * exp(-currentEMP[0]/temperatureSequence[tempindex]);
                MnormsqZ[tempindex][i] += pow(currentNorm[i],2) * exp(-currentEMP[0]/temperatureSequence[tempindex]);
            }
            lattice.GetFFTnorm("pseudospin","x",currentNorm);
            for (int i=0; i<lattice.m_totalSite; i++)
            {
                PxnormZ[tempindex][i] += currentNorm[i] * exp(-currentEMP[0]/temperatureSequence[tempindex]);
                PxnormsqZ[tempindex][i] += pow(currentNorm[i],2) * exp(-currentEMP[0]/temperatureSequence[tempindex]);
            }
            lattice.GetFFTnorm("pseudospin","y",currentNorm);
            for (int i=0; i<lattice.m_totalSite; i++)
            {
                PynormZ[tempindex][i] += currentNorm[i] * exp(-currentEMP[0]/temperatureSequence[tempindex]);
                PynormsqZ[tempindex][i] += pow(currentNorm[i],2) * exp(-currentEMP[0]/temperatureSequence[tempindex]);
            }
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
        for (int i=0; i<lattice.m_totalSite; i++)
        {
            resultNormParam[0][tempindex][i] = MnormsqZ[tempindex][i] / pow(MnormZ[tempindex][i],2) * Z[tempindex];
            resultNormParam[1][tempindex][i] = PxnormsqZ[tempindex][i] / pow(PxnormZ[tempindex][i],2) * Z[tempindex];
            resultNormParam[2][tempindex][i] = PynormsqZ[tempindex][i] / pow(PynormZ[tempindex][i],2) * Z[tempindex];
        }
    }
}

int main()
{
    int latSizeXMax = 3;
    int latSizeYMax = 2;
    string latType = "Triangle";
    double interactionJprime = -1.0;
    double interactionJprimeprime = 0.0;
    double externalXElectricField = 0;
    double externalYElectricField = 0;
    double externalZMageticField = 0;
    string interactionType = "Ising";
    // double assignedDirection[4] = {M_PI/4,M_PI/5,M_PI/6,M_PI/7};    // Spherical Coordination
    double assignedDirection[4] = {0,0,0,0};    // Spherical Coordination

    double Tmin = 0.1;
    double Tmax = 5;
    // double Tstep = 0.1;
    double Tstep = 1;
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
    double theoreticalNormResult[3][TEMPERATURESEQUENCEMAX][LATTICESIZEMAX];

    // double energySequence[simulationCount][(int)(totalSweep*relaxationTime)];

    Lattice lattice(latSizeXMax,latSizeYMax,latType,interactionJprime,interactionJprimeprime,externalXElectricField,externalYElectricField,
        externalZMageticField,interactionType,assignedDirection);
    // lattice.DebugInput(17+256);
    // lattice.DebugOutput();
    // double ans[4];
    // lattice.GetBinaryConfigIsingEMP(13,ans);
    // cout << ans[0] << " " << ans[1] << " " << ans[2] << " " << ans[3] << " " << endl;


    if (interactionType == "Ising")
    {
        CalculateTheoreticalParams(lattice, tempLength, temperature, theoreticalResult, theoreticalNormResult);
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
    
    fout.open("theoreticalFFT.csv",ios::out);
    for (int i=0; i<3; i++)
    {
        for (int j=0; j<tempLength; j++)
        {
            for (int k=0; k<lattice.m_totalSite; k++)
            {
                fout << theoreticalNormResult[i][j][k] << ",";
            }
            fout << endl;
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