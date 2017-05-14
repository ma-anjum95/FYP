#include "mainwindow.h"
#include "ui_mainwindow.h"

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    //MainWindow::makePlot();
}

MainWindow::~MainWindow()
{
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

    //Add value functions of HR,SPO2 etc here
    ui->lcdNumber_HR->display(71.83);
    ui->lcdNumber_RR->display(6.7);
    ui->lcdNumber_SPO2->display(98.97);
    ui->lcdNumber_TEMP->display(97.34);

}
