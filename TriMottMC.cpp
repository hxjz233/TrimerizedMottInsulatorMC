// g++ -o ./TriMottMC.exe ./TriMottMC.cpp -lfftw3 -lm

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <time.h>
#include <random>
#include <cmath>
#include <cstring>
#include <vector>
#include <sstream>
// #include "fftw3.h"

#include "lattice.h"

#define TEMPERATURESEQUENCEMAX 1000
#define LATTICESIZEMAX 10000
using namespace std;


class TheoreticalCalculation
{
public:
    vector<double> resultParam[8];
    vector<vector<double>> resultNormParam[3];
    TheoreticalCalculation(Lattice* lattice, int temperatureLength, double* temperatureInput)
     : m_lattice(lattice), tempNum(temperatureLength), temperatureSequence(temperatureInput)
    {
        for (int i=0; i<8; i++)
        {
            resultParam[i].resize(tempNum);
        }
        
        for (int normindex=0; normindex<3; normindex++)
        {
            resultNormParam[normindex].resize(tempNum);
            for (int tempindex=0; tempindex<tempNum; tempindex++)
            {
                resultNormParam[normindex][tempindex].resize(lattice->m_totalSite);
            }
        }
    }
    void Calculate();
private:
    Lattice* m_lattice;
    int tempNum;
    double* temperatureSequence;
};

void TheoreticalCalculation::Calculate()
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

    double MnormZ[tempNum][m_lattice->m_totalSite] = {0.0};
    double MnormsqZ[tempNum][m_lattice->m_totalSite] = {0.0};
    double PxnormZ[tempNum][m_lattice->m_totalSite] = {0.0};
    double PxnormsqZ[tempNum][m_lattice->m_totalSite] = {0.0};
    double PynormZ[tempNum][m_lattice->m_totalSite] = {0.0};
    double PynormsqZ[tempNum][m_lattice->m_totalSite] = {0.0};

    double currentEMP[4];
    double currentNorm[m_lattice->m_totalSite];

    for (int binConfig=0; binConfig<(1<<(m_lattice->m_totalSite*2)); binConfig++)
    {
        m_lattice->GetBinaryConfigIsingEMP(binConfig,currentEMP);
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

            m_lattice->GetFFTnorm("spin","z",currentNorm);
            for (int i=0; i<m_lattice->m_totalSite; i++)
            {
                MnormZ[tempindex][i] += currentNorm[i] * exp(-currentEMP[0]/temperatureSequence[tempindex]);
                MnormsqZ[tempindex][i] += pow(currentNorm[i],2) * exp(-currentEMP[0]/temperatureSequence[tempindex]);
            }
            m_lattice->GetFFTnorm("pseudospin","x",currentNorm);
            for (int i=0; i<m_lattice->m_totalSite; i++)
            {
                PxnormZ[tempindex][i] += currentNorm[i] * exp(-currentEMP[0]/temperatureSequence[tempindex]);
                PxnormsqZ[tempindex][i] += pow(currentNorm[i],2) * exp(-currentEMP[0]/temperatureSequence[tempindex]);
                // cout << currentNorm[i] << " ";
            }
            // cout << endl;
            m_lattice->GetFFTnorm("pseudospin","y",currentNorm);
            for (int i=0; i<m_lattice->m_totalSite; i++)
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

        for (int i=0; i<m_lattice->m_totalSite; i++)
        {
            resultNormParam[0][tempindex][i] = MnormsqZ[tempindex][i] / pow(MnormZ[tempindex][i],2) * Z[tempindex];
            resultNormParam[1][tempindex][i] = PxnormsqZ[tempindex][i] / pow(PxnormZ[tempindex][i],2) * Z[tempindex];
            resultNormParam[2][tempindex][i] = PynormsqZ[tempindex][i] / pow(PynormZ[tempindex][i],2) * Z[tempindex];
        }
    }
}

int main()
{    
    int latSizeXMax;
    int latSizeYMax;
    string latType;
    double interactionJprime;
    double interactionJprimeprime;
    double externalXElectricField;
    double externalYElectricField;
    double externalZMageticField;
    string interactionType;
    // double assignedDirection[4] = {M_PI/4,M_PI/5,M_PI/6,M_PI/7};    // Spherical Coordination
    double temperature;
    
    int totalSweep = 1E5;
    int sweepPerSample = 50;
    double relaxationTime = 0.1;

    int taskID = 1;
    
    ifstream indata("input.csv",ios::in);
    ofstream fout;
    string outputName;

    if (!indata)
    {
        cout << "file not found!" << endl;
        exit(1);
    }
    string line;
    istringstream sin(line);
    vector<string> fieldTypes;
    string fieldType;

    getline(indata, line);
    sin.clear();
    sin.str(line);
    fieldTypes.clear();

    while (getline(sin, fieldType, ','))
    {
        fieldTypes.push_back(fieldType);
        // cout << fieldType << " " << endl;
    }

    string field;
    while (getline(indata, line))
    {
        clock_t start_time = clock();
        int i=1;
        sin.clear();
        sin.str(line);
        outputName = "simuResult_";
        for (vector<string>::iterator iter=fieldTypes.begin(); iter!=fieldTypes.end(); iter++)
        {
            // cout << i;
            // i++;
            getline(sin,field,',');
            if (*iter == "latSizeXMax")                     latSizeXMax = stoi(field);
            else if (*iter == "latSizeYMax")                latSizeYMax = stoi(field);
            else if (*iter == "latType")                    latType = field;
            else if (*iter == "interactionJprime")          interactionJprime = stod(field);
            else if (*iter == "interactionJprimeprime")     interactionJprimeprime = stod(field);
            else if (*iter == "externalXElectricField")     externalXElectricField = stod(field);
            else if (*iter == "externalYElectricField")     externalYElectricField = stod(field);
            else if (*iter == "externalZMageticField")      externalZMageticField = stod(field);
            else if (*iter == "interactionType")            interactionType = field;
            else if (*iter == "temperature")                temperature = stod(field);
            else
            {
                // cout << i;
                // i++;
                cout << *iter << " ";
                cout << "unknown field!" << endl;
                exit(1);
            }
            outputName += field + "_";
        }
        Lattice lattice(latSizeXMax,latSizeYMax,latType,interactionJprime,interactionJprimeprime,externalXElectricField,externalYElectricField,
            externalZMageticField,interactionType);
        
        lattice.SetTemperature(temperature);

        double configEnergy = 0;
        double sampledEnergy = 0;
        double sampledEsq = 0;
        double configMagnetization = 0;
        double sampledM = 0;
        double sampledMsq = 0;
        double configPolarization[2] = {0};
        double sampledPx = 0;
        double sampledPxsq = 0;
        double sampledPy = 0;
        double sampledPysq = 0;
        double configFFTNorm[3][lattice.m_totalSite] = {0};
        double sampledFFTNorm[6][lattice.m_totalSite] = {0};  // notice only work when concern Mz, Px, Py Binder parameters
        double resultBinder[3][lattice.m_totalSite] = {0};


        int totalSampleNum = 0;
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
                configMagnetization = lattice.GetMagnetization();
                sampledM += configMagnetization;
                sampledMsq += pow(configMagnetization,2);
                configPolarization[0] = lattice.GetPolarizationX();
                sampledPx += configPolarization[0];
                sampledPxsq += pow(configPolarization[0],2);
                configPolarization[1] = lattice.GetPolarizationY();
                sampledPy += configPolarization[1];
                sampledPysq += pow(configPolarization[1],2);
                lattice.GetFFTnorm("spin","z",configFFTNorm[0]);
                lattice.GetFFTnorm("pseudospin","x",configFFTNorm[1]);
                lattice.GetFFTnorm("pseudospin","y",configFFTNorm[2]);
                for (int j=0; j<3; j++)
                {
                    for (int k=0; k<lattice.m_totalSite; k++)
                    {
                        sampledFFTNorm[2*j][k] += configFFTNorm[j][k];
                        sampledFFTNorm[2*j+1][k] += pow(configFFTNorm[j][k],2);
                    }
                }
            }
        }
        sampledEnergy /= totalSampleNum;
        sampledEsq /= totalSampleNum;
        sampledM /= totalSampleNum;
        sampledMsq /= totalSampleNum;
        sampledPx /= totalSampleNum;
        sampledPxsq /= totalSampleNum;
        sampledPy /= totalSampleNum;
        sampledPysq /= totalSampleNum;
        for (int j=0; j<3; j++)
        {
            for (int k=0; k<lattice.m_totalSite; k++)
            {
                resultBinder[j][k] = sampledFFTNorm[2*j+1][k] / pow(sampledFFTNorm[2*j][k],2) * totalSampleNum;
            }
        }
        fout.open(outputName+"Thermo.csv",ios::out);
        fout << sampledEnergy << "," << sampledEsq << "," << sampledM << "," << sampledMsq << "," 
             << sampledPx << "," << sampledPxsq << "," << sampledPy << "," << sampledPysq;
        fout.close();
        fout.open(outputName+"Binder.csv",ios::out);
        for (int i=0; i<3; i++)
        {
            for (int k=0; k<lattice.m_totalSite; k++)
            {
                fout << resultBinder[i][k] << ",";
            }
            fout << endl;
        }
        fout.close();

        clock_t calc_time = clock();
        cout << "Task " << taskID << " completed, elapsed time: " << (double)(calc_time - start_time) / CLOCKS_PER_SEC << "s" << endl;
        taskID++;
    }
}