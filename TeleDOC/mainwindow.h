#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QThread>
#include <QTimer>
#include <string>
#include "ppgworker.h"

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = 0);
    void makePlot(double *ppg_ir, int length);
    ~MainWindow();

signals:
    void start_ppg();
    void stop_ppg();
    void get_ppg_status();
    void save_ppg(const int &id);

private slots:
    //void makePlot();
    void on_start_button_clicked();
    void on_save_button_clicked();
    void update_time();

public slots:
    void handle_ppg_results(double *ppg_red, double *ppg_ir,double hr, double rr, double spo2);
    void ppg_status(const bool &status);

private:
    Ui::MainWindow *ui;
    PPGWorker *ppg_worker;
    QThread ppg_thread;

    // the boolean to keep track if the boolean is working
    bool ppg_working;
    qint64 time_passed;
};

#endif // MAINWINDOW_H
