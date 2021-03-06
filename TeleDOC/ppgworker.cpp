#include "ppgworker.h"

void PPGWorker::ppg_work()
{
    double tmp1, tmp2;

    double *interp_red, *interp_ir;

    std::ifstream file;
    file.open("E:\\patient_001.txt");

    exec_ppg_thread = true;
    emit return_ppg_status(this->exec_ppg_thread);

    while(exec_ppg_thread) {
        // getting the ppg data from the sensors
        file >> tmp1 >> tmp2;
        this->ppg_red.push_back(tmp1);
        this->ppg_ir.push_back(tmp2);
        emit ppg_samples(ppg_red.size());
        QThread::msleep(40);

        while (this->last_index < (signed) this->ppg_red.size() - 500) {
            interp_red = linear_interp_10(this->ppg_red, this->last_index);
            interp_ir = linear_interp_10(this->ppg_ir, this->last_index);

            ppg_analysis.run(interp_red, interp_ir, 5000, 250);

            emit resultsReady(interp_red, interp_ir, ppg_analysis.get_hr(), ppg_analysis.get_rr(),
                              ppg_analysis.get_rr_std(), ppg_analysis.get_spo2());

            this->last_index += update;
        }

        // dont remove this as this will help calling other functions
        QApplication::processEvents();
    }

}

void PPGWorker::ppg_stop()
{
    this->exec_ppg_thread = false;
    emit return_ppg_status(this->exec_ppg_thread);
}

void PPGWorker::ppg_status()
{
    emit return_ppg_status(this->exec_ppg_thread);
}

void PPGWorker::ppg_save(const int &id)
{
    char tmp[22];
    // saves the patients data in a file with given id
    sprintf(tmp, "./data/patient_%03d.txt", id);
    ofstream file;
    file.open(tmp);

    for (int i = 0; i < this->ppg_ir.size(); i++) {
        file << ppg_red[i] << " " << ppg_ir[i] << endl;
    }
    ppg_red.clear();
    ppg_ir.clear();
    this->last_index = 0;
    file.close();
}

double *PPGWorker::linear_interp_10(vector<double> to_interp, int start_index)
{
    // The formula used is:
    // y = (x - x1) (y2 - y1)
    //              ========= + y1
    //              (x2 - x1)
    // we can assume x starts at zero every time so
    // x - x1 = j
    // x2 -x1 = 10

    // be sure to delete it somewhere
    // it will always return an array of 5000.
    // incase the samples are less than 500 we return an error
    double *to_return = new double[5000];

    for (int i = start_index; i < 500 + start_index; i++) {
        for (int j = 0; j < 10; j++) {
            int new_index = (i - start_index) * 10 + j;

            to_return[new_index] = ((double)j / 10.0) * (to_interp[i + 1] - to_interp[i]) + to_interp[i];
        }
    }

    return to_return;
}


void PPGWorker::ppg_update(const int &update)
{
    this->update = update;
}
