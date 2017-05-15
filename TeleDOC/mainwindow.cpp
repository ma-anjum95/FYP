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

    this->ppg_thread.start();

    this->ppg_working = false;
}

MainWindow::~MainWindow()
{
    this->ppg_thread.quit();
    this->ppg_thread.wait();
    delete ui;
}


void MainWindow::makePlot(double *ppg_ir, int length)
{
    double max = -1, min = 180000000;
    int *x = new int[length];

    for (int i = 0; i < length; i++)
    {
        x[i] = i + 1;
        if (max < ppg_ir[i])
            max = ppg_ir[i];

        if (min > ppg_ir[i])
            min = ppg_ir[i];
    }

    //converting array X to QVector
    QVector<double> v(length);
    qCopy(x, x+length, v.begin());

    //converting array Y to QVector
    QVector<double> w(length);
    qCopy(ppg_ir, ppg_ir+length, w.begin());

    //setting values to graph
    ui->CustomPlot->addGraph();
    ui->CustomPlot->graph(0)->setData(v,w);

    //set ranges of axis
    ui->CustomPlot->xAxis->setRange(0, length + 1);
    ui->CustomPlot->yAxis->setRange(min, max);

    //plot
    ui->CustomPlot->replot();
}

void MainWindow::on_pushButton_clicked()
{
    // if the thread i
    if (!this->ppg_working)
        emit start_ppg();
    else
        emit stop_ppg();
}

void MainWindow::handle_ppg_results(double *ppg_red, double *ppg_ir,
                    double hr, double rr, double spo2)
{
    ui->lcdNumber_HR->display(hr);
    ui->lcdNumber_SPO2->display(spo2);
    ui->lcdNumber_RR->display(rr);

    this->makePlot(ppg_ir, 5000);
}

void MainWindow::ppg_status(const bool &status)
{
    this->ppg_working = status;

    if (status == true) {
        ui->pushButton->setText("            STOP              ");
        ui->ppg_status->setText("On");

        QPalette pal = ui->pushButton->palette();
        pal.setColor(QPalette::Button, QColor(Qt::red));
        ui->pushButton->setAutoFillBackground(true);
        ui->pushButton->setPalette(pal);
        ui->pushButton->update();
    }
     else {
        ui->pushButton->setText("            START             ");
        ui->ppg_status->setText("Off");

        QPalette pal = ui->pushButton->palette();
        pal.setColor(QPalette::Button, QColor(Qt::green));
        ui->pushButton->setAutoFillBackground(true);
        ui->pushButton->setPalette(pal);
        ui->pushButton->update();
    }
}
