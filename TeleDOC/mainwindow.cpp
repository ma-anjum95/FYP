#include "mainwindow.h"
#include "ui_mainwindow.h"


MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
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


void MainWindow::makePlot(int y[],int length_y)
{
    //array y for y axis
    //array x for x axis
    int *x = new int[length_y];
//making an array of iterations for the x axis
    for (int i=0; i<length_y; i++)
    {
        x[i] = i;
    }

    //converting array X to QVector
    QVector<double> v(length_y);
    qCopy(x, x+length_y, v.begin());

    //converting array Y to QVector
    QVector<double> w(length_y);
    qCopy(y, y+length_y, w.begin());

    //setting values to graph
    ui->CustomPlot->addGraph();
    ui->CustomPlot->graph(0)->setData(v,w);

    //set names of axis
    ui->CustomPlot->xAxis->setLabel("x");
    ui->CustomPlot->yAxis->setLabel("y");

    //set ranges of axis
    ui->CustomPlot->xAxis->setRange(0,100);
    ui->CustomPlot->yAxis->setRange(0,100);

    //plot
    ui->CustomPlot->replot();
}

void MainWindow::on_pushButton_clicked()
{
    //dummy data for PPG graph
    int x[101];
    for (int i=0; i<101; i++)
    {
        x[i] = i;
    }
    int length_x = sizeof(x)/sizeof(*x);
    //Add ppg value array here with its length
    makePlot(x,length_x);

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
}

void MainWindow::ppg_status(const bool &status)
{
    this->ppg_working = status;

    if (status == true) {
        ui->pushButton->setText("            STOP              ");
        ui->ppg_status->setText("On");
    }
     else {
        ui->pushButton->setText("            START             ");
        ui->ppg_status->setText("Off");

    }
}
