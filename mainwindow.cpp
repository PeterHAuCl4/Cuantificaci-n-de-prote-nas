#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QFile>
#include <QTextStream>
#include <QMessageBox>
#include <math.h>

QString s1;
double beta1 = 0, beta0 = 0 ,  r = 0;
MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
//    MainWindow::makePlot();
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_Solver_clicked()
{
    bool convert;
    int n = ui->data->text().toInt(&convert,10) , i = 0 , inX = ui->interX->text().toInt(&convert,10)
            , inY = ui->interY->text().toInt(&convert,10) , inX_2 = ui->interX_2->text().toInt(&convert,10)
            , inY_2 = ui->interY_2->text().toInt(&convert,10);
    QVector<double> todo, x , y , y1;
    QString archivo;
    QVector<QString> fagocitos;
    QFile file(s1);


    if(!file.open(QIODevice::ReadOnly))
        QMessageBox::information(0,"Info",file.errorString());
    QTextStream in(&file);

    while(!in.atEnd()){
        archivo = in.readLine();
        fagocitos += archivo;
    }
    todo = MainWindow::toFloatVector(fagocitos);

    for(i = 0; i < 2*n; ++i){
        if(i < n){
            x << todo[i];
        }else{
            y << todo[i];
        }
    }
    //===========Se resuelve la ecuación y se grafica===================

for(double j = x[0]; j <= x[n-1]; j+=ui->saltos->text().toDouble(&convert)){
        y1 << interpolate(x,y,j,n);
}
    // create graph and assign data to it:
    QCPGraph *graph2 = ui->customPlot->addGraph();
    graph2->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 1.5), QBrush(Qt::white), 9));
    graph2->setPen(QPen(QColor(0,191,17)));
    ui->customPlot->graph(0)->setData(x, y);
    QCPGraph *graph1 = ui->customPlot->addGraph();
    graph1->setData(x, y1);
    graph1->setPen(QPen(QColor(255, 243, 0)));
    // give the axes some labels:
    ui->customPlot->xAxis->setLabel("x");
    ui->customPlot->yAxis->setLabel("y");
    ui->customPlot->xAxis->setBasePen(QPen(Qt::white, 1));
    ui->customPlot->yAxis->setBasePen(QPen(Qt::white, 1));
    ui->customPlot->xAxis->setTickPen(QPen(Qt::white, 1));
    ui->customPlot->yAxis->setTickPen(QPen(Qt::white, 1));
    ui->customPlot->xAxis->setSubTickPen(QPen(Qt::white, 1));
    ui->customPlot->yAxis->setSubTickPen(QPen(Qt::white, 1));
    ui->customPlot->xAxis->setTickLabelColor(Qt::white);
    ui->customPlot->yAxis->setTickLabelColor(Qt::white);
    ui->customPlot->xAxis->grid()->setPen(QPen(QColor(140, 140, 140), 1, Qt::DotLine));
    ui->customPlot->yAxis->grid()->setPen(QPen(QColor(140, 140, 140), 1, Qt::DotLine));
    ui->customPlot->xAxis->grid()->setSubGridPen(QPen(QColor(80, 80, 80), 1, Qt::DotLine));
    ui->customPlot->yAxis->grid()->setSubGridPen(QPen(QColor(80, 80, 80), 1, Qt::DotLine));
    ui->customPlot->xAxis->grid()->setSubGridVisible(true);
    ui->customPlot->yAxis->grid()->setSubGridVisible(true);
    ui->customPlot->xAxis->grid()->setZeroLinePen(Qt::NoPen);
    ui->customPlot->yAxis->grid()->setZeroLinePen(Qt::NoPen);
    ui->customPlot->xAxis->setUpperEnding(QCPLineEnding::esSpikeArrow);
    ui->customPlot->yAxis->setUpperEnding(QCPLineEnding::esSpikeArrow);
    QLinearGradient plotGradient;
    plotGradient.setStart(0, 0);
    plotGradient.setFinalStop(0, 350);
    plotGradient.setColorAt(0, QColor(80, 80, 80));
    plotGradient.setColorAt(1, QColor(50, 50, 50));
    ui->customPlot->setBackground(plotGradient);
    QLinearGradient axisRectGradient;
    axisRectGradient.setStart(0, 0);
    axisRectGradient.setFinalStop(0, 350);
    axisRectGradient.setColorAt(0, QColor(80, 80, 80));
    axisRectGradient.setColorAt(1, QColor(30, 30, 30));
    ui->customPlot->axisRect()->setBackground(axisRectGradient);
    // set axes ranges, so we see all data:
    ui->customPlot->xAxis->setRange(inX, inX_2);
    ui->customPlot->yAxis->setRange(inY, inY_2);
    ui->customPlot->replot();

}
double MainWindow::interpolate(QVector<double> x, QVector<double> y, double xi, int n)
{
    double result = 0; // Initialize result

    for (int i=0; i<n; i++)
    {
        // Compute individual terms of above formula
        double term = y[i];
        for (int j=0;j<n;j++)
        {
            if (j!=i)
                term = term*(xi - x[j])/(x[i] - x[j]);
        }

        // Add current term to result
        result += term;
    }

    return result;
}
void MainWindow::solveEc(QVector<double> x, QVector<double> y, int n){

    int i = 0;
    double sumx = 0, sumy = 0 , sumxy = 0, sumx2 = 0, sumy2 = 0 , sxx = 0 , sxy = 0;
    for(i = 0; i < n; ++i)
        sumx += x[i];
    for(i = 0; i < n; ++i)
        sumy += y[i];
    for(i = 0; i < n; ++i)
        sumxy += (y[i]*x[i]);
    for(i = 0; i < n; ++i)
        sumx2 += pow(x[i],2);
    for(i = 0; i < n;++i)
        sumy2 += pow(y[i],2);
    sxx = sumx2 - (pow(sumx,2))/n;
    sxy = sumxy - (sumx*sumy)/n;

    beta1 = sxy/sxx;
    beta0 = (sumy/n) - (sxy/sxx)*(sumx/n);
    r = ((n*sumxy)-(sumx*sumy))/sqrt((n*sumx2-(pow(sumx,2)))*(n*sumy2-(pow(sumy,2))));
}
void MainWindow::on_inData_clicked()
{
    s1 =
            QFileDialog::getOpenFileName(this, "Open a file", "directoryToOpen",
                "Text files (*.txt);;Images (*.png *.xpm *.jpg);;XML files (*.xml)");

    ui->Output->setText(s1);
}
QVector<double> MainWindow::toFloatVector(QVector<QString> &aVector){
    QVector<double> vector;

    for(const auto& item : aVector){
        bool ok = true;
        const double value = item.toDouble(&ok);
        if(ok)
            vector << value;
    }
    return vector;
}
void MainWindow::on_pushButton_clicked()
{
    bool convert;
    int n = ui->data->text().toInt(&convert,10) , i = 0;
    double inX = ui->interX->text().toDouble(&convert) , inY = ui->interY->text().toDouble(&convert)
            , inX_2 = ui->interX_2->text().toDouble(&convert)
            , inY_2 = ui->interY_2->text().toDouble(&convert);
    QVector<double> todo, x , y , y1;
    QString archivo;
    QVector<QString> fagocitos;
    QFile file(s1);

    if(!file.open(QIODevice::ReadOnly))
        QMessageBox::information(0,"Info",file.errorString());
    QTextStream in(&file);

    while(!in.atEnd()){
        archivo = in.readLine();
        fagocitos += archivo;
    }
    todo = MainWindow::toFloatVector(fagocitos);

    for(i = 0; i < 2*n; ++i){
        if(i < n){
            x << todo[i];
        }else{
            y << todo[i];
        }
    }
    solveEc(x,y,n);

    for(double j = x[0]; j <= x[n-1]; j+=ui->saltos->text().toDouble(&convert)){
            y1 << beta1*j+beta0;
    }

    // create graph and assign data to it:
    QCPGraph *graph2 = ui->customPlot->addGraph();
    graph2->setScatterStyle(QCPScatterStyle(QCPScatterStyle::ssCircle, QPen(Qt::black, 1.5), QBrush(Qt::white), 9));
    graph2->setPen(QPen(QColor(0,191,17)));
    ui->customPlot->graph(0)->setData(x, y);
    QCPGraph *graph1 = ui->customPlot->addGraph();
    graph1->setData(x, y1);
    graph1->setPen(QPen(QColor(255, 243, 0)));
    // give the axes some labels:
    ui->customPlot->xAxis->setLabel("x");
    ui->customPlot->yAxis->setLabel("y");
    ui->customPlot->xAxis->setBasePen(QPen(Qt::white, 1));
    ui->customPlot->yAxis->setBasePen(QPen(Qt::white, 1));
    ui->customPlot->xAxis->setTickPen(QPen(Qt::white, 1));
    ui->customPlot->yAxis->setTickPen(QPen(Qt::white, 1));
    ui->customPlot->xAxis->setSubTickPen(QPen(Qt::white, 1));
    ui->customPlot->yAxis->setSubTickPen(QPen(Qt::white, 1));
    ui->customPlot->xAxis->setTickLabelColor(Qt::white);
    ui->customPlot->yAxis->setTickLabelColor(Qt::white);
    ui->customPlot->xAxis->grid()->setPen(QPen(QColor(140, 140, 140), 1, Qt::DotLine));
    ui->customPlot->yAxis->grid()->setPen(QPen(QColor(140, 140, 140), 1, Qt::DotLine));
    ui->customPlot->xAxis->grid()->setSubGridPen(QPen(QColor(80, 80, 80), 1, Qt::DotLine));
    ui->customPlot->yAxis->grid()->setSubGridPen(QPen(QColor(80, 80, 80), 1, Qt::DotLine));
    ui->customPlot->xAxis->grid()->setSubGridVisible(true);
    ui->customPlot->yAxis->grid()->setSubGridVisible(true);
    ui->customPlot->xAxis->grid()->setZeroLinePen(Qt::NoPen);
    ui->customPlot->yAxis->grid()->setZeroLinePen(Qt::NoPen);
    ui->customPlot->xAxis->setUpperEnding(QCPLineEnding::esSpikeArrow);
    ui->customPlot->yAxis->setUpperEnding(QCPLineEnding::esSpikeArrow);
    QLinearGradient plotGradient;
    plotGradient.setStart(0, 0);
    plotGradient.setFinalStop(0, 350);
    plotGradient.setColorAt(0, QColor(80, 80, 80));
    plotGradient.setColorAt(1, QColor(50, 50, 50));
    ui->customPlot->setBackground(plotGradient);
    QLinearGradient axisRectGradient;
    axisRectGradient.setStart(0, 0);
    axisRectGradient.setFinalStop(0, 350);
    axisRectGradient.setColorAt(0, QColor(80, 80, 80));
    axisRectGradient.setColorAt(1, QColor(30, 30, 30));
    ui->customPlot->axisRect()->setBackground(axisRectGradient);
    // set axes ranges, so we see all data:
    ui->customPlot->xAxis->setRange(inX, inX_2);
    ui->customPlot->yAxis->setRange(inY, inY_2);
    ui->customPlot->replot();
    QString ssa = "La ecuación de la recta es: y = " + QString::number(beta1)  + "x + " + QString::number(beta0)
            + "           |         Coeficiente de correlación: " + QString::number(r);
    ui->Output->setText(ssa);
}
void MainWindow::on_BtnAbs_clicked()
{
    bool ok;
    if(beta0 == 0.0 && beta1 == 0.0){
        QMessageBox::warning(0,"¡Advertencia!", "No se puede realizar ese cálculo");
    }else{
        double resul = (ui->Abs->text().toDouble(&ok)-beta0)/beta1;
        QString sa = "La concentración es: " + QString::number(resul);
        QMessageBox::information(0,"Info",sa);
    }
}
