#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{   
    ui->setupUi(this);

    ui->CustomPlot->xAxis->setVisible(false);

    this->ppg_worker = new PPGWorker;
    this->ppg_worker->moveToThread(&this->ppg_thread);

    connect(&this->ppg_thread, &QThread::finished, this->ppg_worker, &QObject::deleteLater);
    connect(this, &MainWindow::start_ppg, this->ppg_worker, &PPGWorker::ppg_work);
    connect(this, &MainWindow::stop_ppg, this->ppg_worker, &PPGWorker::ppg_stop);
    connect(this->ppg_worker, &PPGWorker::resultsReady,
            this, &MainWindow::handle_ppg_results);
    connect(this, &MainWindow::get_ppg_status,
            this->ppg_worker, &PPGWorker::ppg_status);
    connect(this->ppg_worker, &PPGWorker::return_ppg_status,
            this, &MainWindow::ppg_status);
    connect(this, &MainWindow::save_ppg,
            this->ppg_worker, &PPGWorker::ppg_save);
    connect(this->ppg_worker, &PPGWorker::ppg_samples,
            this, &MainWindow::ppg_samples);

    // adding update radio buttons in a qbuttongroup
    this->update_radio_buttons.addButton(ui->radio_1);
    this->update_radio_buttons.addButton(ui->radio_5);
    this->update_radio_buttons.addButton(ui->radio_20);

    connect(ui->radio_1, &QRadioButton::toggled, this, &MainWindow::radio_toggled);
    connect(ui->radio_5, &QRadioButton::toggled, this, &MainWindow::radio_toggled);
    connect(ui->radio_20, &QRadioButton::toggled, this, &MainWindow::radio_toggled);
    connect(this, &MainWindow::ppg_update, this->ppg_worker, &PPGWorker::ppg_update);
    ui->radio_5->setChecked(true);

    this->ppg_signal_buttons.addButton(ui->radio_ir);
    this->ppg_signal_buttons.addButton(ui->radio_red);

    connect(ui->radio_ir, &QRadioButton::toggled, this, &MainWindow::radio_toggled);
    connect(ui->radio_red, &QRadioButton::toggled, this, &MainWindow::radio_toggled);
    ui->radio_ir->setChecked(true);

    connect(this->ppg_worker, &PPGWorker::ppg_device_fail,
            this, &MainWindow::ppg_device_fail);

    this->ppg_thread.start();

    this->ppg_working = false;

}

MainWindow::~MainWindow()
{
    this->ppg_thread.quit();
    this->ppg_thread.wait();
    delete ui;
}


void MainWindow::makePlot(double *ppg, int length)
{
    double max = -1, min = 180000000;
    int *x = new int[length];

    for (int i = 0; i < length; i++)
    {
        x[i] = i + 1;
        if (max < ppg[i])
            max = ppg[i];

        if (min > ppg[i])
            min = ppg[i];
    }

    //converting array X to QVector
    QVector<double> v(length);
    qCopy(x, x+length, v.begin());

    //converting array Y to QVector
    QVector<double> w(length);
    qCopy(ppg, ppg+length, w.begin());

    //setting values to graph
    ui->CustomPlot->addGraph();
    ui->CustomPlot->graph(0)->setData(v,w);

    //set ranges of axis
    ui->CustomPlot->xAxis->setRange(0, length + 1);
    ui->CustomPlot->yAxis->setRange(min, max);

    //plot
    ui->CustomPlot->replot();
}

void MainWindow::on_start_button_clicked()
{
    // if the thread i
    if (!this->ppg_working) {
        emit start_ppg();
    }
    else {
        emit stop_ppg();

    }
}

void MainWindow::handle_ppg_results(double *ppg_red, double *ppg_ir,
                    double hr, double rr, double rr_std, double spo2, bool anomaly)
{
    ui->lcdNumber_HR->display(hr);
    ui->lcdNumber_SPO2->display(spo2);
    ui->lcdNumber_RR->display(rr);
    ui->lcdNumber_RR_STD->display(rr_std);

    if (this->display_ir)
        this->makePlot(ppg_ir, 5000);
    else
        this->makePlot(ppg_red, 5000);

    if (anomaly)
        ui->anomaly_status->setText("Yes");
    else
        ui->anomaly_status->setText("No");
}

void MainWindow::ppg_status(const bool &status)
{ static int i = 10;
    this->ppg_working = status;

    if (status) {
        ui->start_button->setText("            STOP              ");

        QPalette pal = ui->start_button->palette();
        pal.setColor(QPalette::Button, QColor(Qt::red));
        ui->start_button->setAutoFillBackground(true);
        ui->start_button->setPalette(pal);
        ui->start_button->update();

        ui->save_button->setEnabled(false);
    }
     else {
        ui->start_button->setText("            START             ");

        QPalette pal = ui->start_button->palette();
        pal.setColor(QPalette::Button, QColor(Qt::green));
        ui->start_button->setAutoFillBackground(true);
        ui->start_button->setPalette(pal);
        ui->start_button->update();

        ui->save_button->setEnabled(true);
    }
}


void MainWindow::on_save_button_clicked()
{
    char tmp[3];

    static int id = 0;

    emit save_ppg(id++);

    // sets the patients id
    sprintf(tmp, "%03d", id);
    ui->patient_id->setText(tmp);
}


void MainWindow::ppg_samples(const int &samples)
{
    char tmp[10];
    sprintf(tmp, "%d", samples);
    ui->ppg_samples->setText(tmp);
}

void MainWindow::radio_toggled()
{
    if (ui->radio_1->isChecked())
        emit ppg_update(1 * 25);
    else if (ui->radio_5->isChecked())
        emit ppg_update(5 * 25);
    else if (ui->radio_20->isChecked())
        emit ppg_update(20 * 25);

    if (ui->radio_ir->isChecked())
        this->display_ir = true;
    else if (ui->radio_red->isChecked())
        this->display_ir = false;
}


void MainWindow::ppg_device_fail()
{
    QMessageBox::critical(this, tr("MAX117 Initialization failure"),
                       tr("The ppg device MAX117 failed to initialize. Please make sure it is connected."));

}
