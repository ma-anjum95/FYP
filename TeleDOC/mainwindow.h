#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QThread>
#include <QTimer>
#include <string>
#include <QButtonGroup>
#include <QMessageBox>
#include "ppgworker.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);

    ~MainWindow();

signals:
    void start_ppg();
    void stop_ppg();
    void get_ppg_status();
    void save_ppg(const int &id);
    void ppg_update(const int &update);

private slots:
    //void makePlot();
    void on_start_button_clicked();
    void on_save_button_clicked();
    void radio_toggled();

public slots:
    void handle_ppg_results(double *ppg_red, double *ppg_ir,double hr, double rr, double rr_std, double spo2);
    void ppg_status(const bool &status);
    void ppg_samples(const int &samples);
    void ppg_device_fail();

private:
    Ui::MainWindow *ui;
    PPGWorker *ppg_worker;
    QThread ppg_thread;
    QButtonGroup update_radio_buttons;
    QButtonGroup ppg_signal_buttons;

    // the boolean to keep track if the boolean is working
    bool ppg_working;
    qint64 time_passed;
    bool display_ir = true; // displays ir on graph if true

    void makePlot(double *ppg, int length);
};

#endif // MAINWINDOW_H
