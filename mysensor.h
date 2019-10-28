#ifndef MYSENSOR_H
#define MYSENSOR_H

#include <QMainWindow>

#include <QGraphicsView>
#include <QWidget>
#include <QGraphicsScene>
#include <QtGui>
#include <QtCore>
#include <QDialog>
#include <QString>
#include <string>
#include <gps_receiver.h>
#include <Kriging_Lib.h>

#define LDR_PATH "/sys/bus/iio/devices/iio:device0/in_voltage"

namespace Ui {
class MySensor;
}

class MySensor : public QMainWindow
{
    Q_OBJECT

public:
    explicit MySensor(QWidget *parent = 0);
    ~MySensor();     
    bool pointInPolygon(double xx, double yy);
    int readAnalog(int number);
    double readingECa();

    QGraphicsScene *scene;
    QGraphicsScene *sceneZmap;
    QGraphicsScene *sceneEmap;
    Kriging_Lib *startKrig;
    Kriging_Lib *Exporting;
    QPixmap Zpixmap, Epixmap;
    double zmax, zmin, ezmin, ezmax;

    GPS_Receiver *gps;

    int ii, NC, NContour, NG,typePoint,sample;
    double pixelSize;
    double newxUTM, newyUTM;
    double sumsqECa;
    double sumECa;
    int nECa;
    int cECa;
    double varECa;
    double dxECa;
    double dyECa;
    double antECa;
    double dsumECa;
    double dsumsqECa;
    int dnECa;
    double dvarECa;
    int dcECa;

    QString statusGPS;    
    double reqDist;
    double distFactor;
    double signalF;
    double distanceToNP;
    double ecaValue;
    double maxxUTM;
    double minxUTM;
    double maxyUTM;
    double minyUTM;
    int zone;
    int iSize;
    int jSize;
    int iLoc;
    int jLoc;

    double CxPos[1000];
    double CyPos[1000];
    double CzPos[1000];
    bool ICPos[1000];
    double GCxPos[1000];
    double GCyPos[1000];

    bool already_clicked;
    double x1, x2, x3;
    double y1, y2, y3;
    double dx,dy;
    double d;
    double A[2][2];
    double xUTM, x;
    double yUTM, y;
    double alt;

    //QString timeString;
    double lat, longi, data,hora,pdop;
    int ihmed;   // mean value of current read values that are higher than the mean value
    int ilmed;   // mean value of current read values that are lower than the mean value
    int vhmed;   // mean value of voltage read values that are higher than the mean value
    int vlmed;   // mean value of voltage read values that are lower than the mean value
    bool contourIsDone;



private slots:
    void on_pushButton_clicked();
    void on_pushButton_2_clicked();
    void on_pushButton_3_clicked();
    void on_pushButton_4_clicked();
    void on_pushButton_5_clicked();
    void on_pushButton_6_clicked();
    void show_picture();

private:
    Ui::MySensor *ui;
};

#endif // MYSENSOR_H
