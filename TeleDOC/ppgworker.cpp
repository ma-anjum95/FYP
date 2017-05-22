#include "ppgworker.h"

PPGWorker::PPGWorker()
{
    this->handle = serialOpen ("/dev/ttyAMA0", 9600) ;



}

PPGWorker::~PPGWorker()
{
    serialClose(this->handle);
}

void PPGWorker::ppg_work()
{
    uint32_t tmp1, tmp2;
    double *interp_red, *interp_ir;

    if (!maxim_max30102_init() || this->handle == -1) {
        emit ppg_device_fail();
    } else {
        exec_ppg_thread = true;
        emit return_ppg_status(this->exec_ppg_thread);

    }

    while(exec_ppg_thread) {
        // getting the ppg data from the sensors
        if (maxim_max30102_read_fifo(&tmp1, &tmp2)) {
            this->ppg_red.push_back((double)tmp1);
            this->ppg_ir.push_back((double)tmp2);

            emit ppg_samples(ppg_red.size());


            while (this->last_index < (signed) (this->ppg_red.size() - 500)) {
                interp_red = linear_interp_10(this->ppg_red, this->last_index);
                interp_ir = linear_interp_10(this->ppg_ir, this->last_index);

                ppg_analysis.run(interp_red, interp_ir, 5000, 250);

                emit resultsReady(interp_red, interp_ir, ppg_analysis.get_hr(), ppg_analysis.get_rr(),
                                  ppg_analysis.get_rr_std(), ppg_analysis.get_spo2(),
                                  this->anomaly(ppg_analysis.get_hr(), ppg_analysis.get_rr(), ppg_analysis.get_rr_std()));

                this->sendSerialData(ppg_analysis.get_hr(), ppg_analysis.get_rr(),
                                     ppg_analysis.get_rr_std(), ppg_analysis.get_spo2());

                this->last_index += update;
                QApplication::processEvents();
            }

            // dont remove this as this will help calling other functions
            QApplication::processEvents();
        }
    }
}

void PPGWorker::sendSerialData(double hr, double rr, double rr_dev, double spo2)
{
    serialPrintf(this->handle,"{\"hr\": %f,\"rr\": %f,\"temp\": 0,\"spo2\": %f, \"rr_dev\": %f}", hr, rr, spo2, rr_dev);
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

double PPGWorker::anomaly_gauss(double x, double mu, double sig)
{
    return (1/(sqrt(2 * M_PI) * sig)) * exp(-((x - mu) * (x - mu))/(2.0 * sig * sig));
}

bool PPGWorker::anomaly(double hr, double rr, double rr_dev)
{
    double const HR_MEAN = 124.595617240614;
    double const RR_MEAN = 9.05159415743016;
    double const RR_DEV_MEAN = 2.56199101805027;

    double const HR_STD = 26.3494843832011;
    double const RR_STD = 2.73808772977394;
    double const RR_DEV_STD = 2.11735336729146;

    double prob_HR = this->anomaly_gauss(hr, HR_MEAN, HR_STD) * 10;
    double prob_RR = this->anomaly_gauss(rr, RR_MEAN, RR_STD) * 10;
    double prob_RR_STD = this->anomaly_gauss(rr_dev, RR_DEV_MEAN, RR_DEV_STD) * 10;

    double final_prob = prob_HR * prob_RR * prob_RR_STD;

    // if there is no anomaly return false else return true
    return final_prob > 0.04 ? false : true;
}
