#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>

namespace Ui {
class MainWindow;
}

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    explicit MainWindow(QWidget *parent = nullptr);
    ~MainWindow();

private slots:
    void makePlot();
    void on_Solver_clicked();
    void solveEc(QVector<double> x, QVector<double> y, int n);
    void on_inData_clicked();
    QVector<double> toFloatVector(QVector<QString> &aVector);
    double interpolate(QVector<double> x, QVector<double> y, double xi, int n);

    void on_interX_textChanged(const QString &arg1);

    void on_pushButton_clicked();

    void on_BtnAbs_clicked();

private:
    Ui::MainWindow *ui;
};

#endif // MAINWINDOW_H
