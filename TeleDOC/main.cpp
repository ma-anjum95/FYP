#include "mainwindow.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MainWindow w;
   // w.makePlot();
    w.resize(800, 480);
    w.show();

    return a.exec();
}


/*QVector<double> x(101), y(101);
for (int i=0; i<101; i++)
{
    x[i] = i/50.0 -1;
    y[i] = x[i]*x[i];
}*/
