#ifndef PPGWORKER_H
#define PPGWORKER_H

#include <QObject>
#include <QThread>
#include <QtTest/QtTest>

#include <fstream>
#include "PPG_C++/ppg_analysis.h"

class PPGWorker : public QObject
{
    Q_OBJECT
signals:
    void resultsReady(double *ppg_red, double *ppg_ir, double hr, double rr, double spo2);
    void return_ppg_status(const bool &status);
    void ppg_samples(const int &samples);

public slots:
    void ppg_work();
    void ppg_stop();
    void ppg_status();
    void ppg_save(const int &id);
    void ppg_update(const int &update);

private:
    PPG::PPGAnalysis ppg_analysis;
    bool exec_ppg_thread = false;
    double *linear_interp_10(vector<double> to_interp, int start_index);
    vector<double> ppg_red;
    vector<double> ppg_ir;
    signed int last_index = 0, update=125;

};

#endif // PPGWORKER_H
