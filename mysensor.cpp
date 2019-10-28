#include "mysensor.h"
#include "ui_mysensor.h"
#include <gps_receiver.h>
#include <Kriging_Lib.h>
#include "BlackLib.h"
#include "BlackCore.h"
#include "BlackDef.h"
#include "BlackErr.h"
#include "BlackPWM.h"
#include <QTimer>
#include <QGraphicsTextItem>
#include <QGraphicsScene>
#include <QFile>
#include <QString>
#include <QStringList>
#include <QTextStream>
#include <QtCore>
#include <QtGui>
#include <QDebug>
#include <QMessageBox>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include "time.h"
using namespace std;
using namespace BlackLib;
MySensor::MySensor(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MySensor)
{
    ui->setupUi(this);
    NC=0;                   // total number of collected points 
    NContour=0;             // number of contour points collected
    NG=0;                   // number of points inside the polygon collected
    pixelSize=0.25;         // size of the pixel, m
    newxUTM=-9999.0;        // coordinates easting of next point
    newyUTM=-9999.0;        // coordinates northing of next point
    sumsqECa=0;             // sum of squared values of ECa
    sumECa=0;               // sum of ECa values
    nECa=0;                 // number of ECa values measured
    cECa=0;                 //
    varECa=0;               // variance of ECa values
    distFactor=1.0;
    contourIsDone=false;
    statusGPS="Void";
//
    if(contourIsDone==false){typePoint=1;}
//
    ui->lineEdit->setText(statusGPS);               // Show GPS Status on screen
    distanceToNP=ui->lineEdit_3->text().toDouble(); // Distance to the next point, m
    ecaValue=ui->lineEdit_4->text().toDouble();     // Measured ECa value
//
    gps->GetGPSPosition(statusGPS, lat, longi, x, y, alt, hora, data, pdop, zone);
//
    maxxUTM=x+50.0;
    minxUTM=x-50.0;
    maxyUTM=y+50.0;
    minyUTM=y-50.0;
//
    iSize=(int)((maxxUTM-minxUTM)/pixelSize);
    jSize=(int)((maxyUTM-minyUTM)/pixelSize);
//
    scene = new QGraphicsScene(this);     
    sceneZmap = new QGraphicsScene(this); 
    sceneEmap = new QGraphicsScene(this); 
    ui->graphicsView->setScene(scene);    
//
    QTimer *timer = new QTimer();
    QObject::connect(timer,SIGNAL(timeout()),this,SLOT(show_picture())); 
    timer->start(2000);
}
MySensor::~MySensor()
{
    delete ui;
}
void MySensor::on_pushButton_3_clicked() //SETUP CONFIGURATION
{
    // 
    reqDist=ui->lineEdit_2->text().toDouble();   
    signalF=ui->lineEdit_5->text().toDouble();  
    BlackPWM myPWM(P8_13);
    myPWM.setDutyPercent(100);
    myPWM.setPeriodTime(1000.0/signalF,milisecond);
    myPWM.setDutyPercent(50);
    sleep(1);
    QMessageBox::information(this,"Setup Configuration","The distance between points and frenquency \n of the signal were configured successfully!!!");
}
// function that shows the sampled points in the screen of the graphical interface
void MySensor::show_picture(){
   int ixx, iyy;
   int ixx1, iyy1;
   int ixx2, iyy2;
   double myslider;
   scene->clear();
//
   gps->GetGPSPosition(statusGPS, lat, longi, x, y, alt, hora, data, pdop, zone);
   ui->lineEdit->setText(statusGPS);
//
   xUTM=x;
   yUTM=y;
//
   myslider=ui->horizontalSlider->value();
   pixelSize=0.1+myslider*5/99;
//
   if (minxUTM-x<2*reqDist) {
       minxUTM=x-5*reqDist;
   }
   if (x-maxxUTM<2*reqDist) {
       maxxUTM=x+5*reqDist;
   }
   if (minyUTM-y<2*reqDist) {
       minyUTM=y-5*reqDist;
   }
   if (y-maxyUTM<2*reqDist) {
       maxyUTM=y+5*reqDist;
   }
   iSize=(int)((maxxUTM-minxUTM)/pixelSize);
   jSize=(int)((maxyUTM-minyUTM)/pixelSize);
   ixx1=(int) ((xUTM-minxUTM)/pixelSize);
   iyy1=(int) ((yUTM-minyUTM)/pixelSize);
   ixx2=0;
   iyy2=0;
//
   if(newxUTM!=-9999.0) {
       ixx2=(int) ((newxUTM-minxUTM)/pixelSize);
       iyy2=(int) ((newyUTM-minyUTM)/pixelSize);
   }
//
   QPixmap pix(iSize,jSize);
   pix.fill(Qt::white);
   QPixmap rotate(pix.size());
   QPainter *paint = new QPainter(&rotate);
   paint->rotate(-90);
// 
   paint->setPen(QColor(255,0,0,255));
   paint->drawPoint(ixx1,iyy1);
   paint->drawPoint(ixx1-1,iyy1);
   paint->drawPoint(ixx1+1,iyy1);
   paint->drawPoint(ixx1-2,iyy1);
   paint->drawPoint(ixx1+2,iyy1);
   paint->drawPoint(ixx1,iyy1-1);
   paint->drawPoint(ixx1,iyy1+1);
   paint->drawPoint(ixx1,iyy1-2);
   paint->drawPoint(ixx1,iyy1+2);
//
   if(newxUTM!=-9999.0){
       paint->setPen(QColor(0,255,0,255));
       paint->drawPoint(ixx2,iyy2);
       paint->drawPoint(ixx2-1,iyy2);
       paint->drawPoint(ixx2+1,iyy2);
       paint->drawPoint(ixx2-2,iyy2);
       paint->drawPoint(ixx2+2,iyy2);
       paint->drawPoint(ixx2,iyy2-1);
       paint->drawPoint(ixx2,iyy2+1);
       paint->drawPoint(ixx2-1,iyy2);
       paint->drawPoint(ixx2+1,iyy2);
       paint->drawPoint(ixx2-2,iyy2);
       paint->drawPoint(ixx2+2,iyy2);
    }
    paint->setPen(QColor(0,0,255,255));
    for (int i=0;i<NC;i++){
        ixx=(int) ((CxPos[i]-minxUTM)/pixelSize);
        iyy=(int) ((CyPos[i]-minyUTM)/pixelSize);
        paint->drawPoint(ixx,iyy);
        paint->drawPoint(ixx-1,iyy);
        paint->drawPoint(ixx+1,iyy);
        paint->drawPoint(ixx-2,iyy);
        paint->drawPoint(ixx+2,iyy);
        paint->drawPoint(ixx,iyy-1);
        paint->drawPoint(ixx,iyy+1);
        paint->drawPoint(ixx,iyy-2);
        paint->drawPoint(ixx,iyy+2);
    }
    paint->setPen(QColor(0,0,0,255));
    for (int i=0;i<NG;i++){
        ixx=(int) ((GCxPos[i]-minxUTM)/pixelSize);
        iyy=(int) ((GCyPos[i]-minyUTM)/pixelSize);
        paint->drawPoint(ixx,iyy);
        paint->drawPoint(ixx-1,iyy);
        paint->drawPoint(ixx+1,iyy);
        paint->drawPoint(ixx,iyy-1);
        paint->drawPoint(ixx,iyy+1);
     }
     delete paint;
     scene->addPixmap(pix);
     ui->graphicsView->centerOn(QPointF(ixx2, iyy2));
     QCoreApplication::processEvents();
}
double MySensor::readingECa(){
    // Read ECa
    // P9-35 AIN6 will be used for reading the electrical current applied
    // P9-36 AIN5 will be used for reading the voltage difference
    // do not forget to setup beaglebone
    // export SLOTS=/sys/devices/bone_capemgr.9/slots
    // sudo sh -c "echo BB-UART1 > $SLOTS"
    // sudo sh -c "echo BB-ADC > $SLOTS"
    double il;  // current low
    double ih;  // current high
    double vl;  // voltage low
    double vh;  // voltage high
    double a=0.3; // distance between electrodes, in meters
    double ECa;
    int j1, j2, j3, j4;
    int i_value[100];
    int v_value[100];
    int imed;    // mean value read for current
    int vmed;    // mean value read for voltage
    double R1, R2, R3, R4, R0;
    bool okFlag;
    QString myMessage;
    R1=9.2;   //kohms
    R2=1.5;   //kohms
    R3=4.7;   //kohms
    R4=7.1;   //kohms
    R0 = 180.0;  // ohms, the resistance used for current measurement;
    //calculating the current and voltage
    imed=0;
    vmed=0;
    j1=0;
    while(j1<100){
        usleep(1000);
        i_value[j1] = readAnalog(6)*1000*1.8/4095;
        usleep(1000);
        v_value[j1] = readAnalog(5)*1000*1.8/4095;
        imed=imed+i_value[j1];        
        vmed=vmed+v_value[j1];        
        j1=j1+1;
    }
    imed=imed/100;
    vmed=vmed/100;
    j1=0;  // number of point with current greater than the mean value
    j2=0;  // number of point with current lower than the mean value
    j3=0;  // number of point with voltage greater than the mean value
    j4=0;  // number of point with voltage lower than the mean value
    ihmed=0;
    ilmed=0;
    vhmed=0;
    vlmed=0;
    for(int j=0;j<100;j++){
       if(i_value[j]>imed){
          ihmed=ihmed+i_value[j];
          j1=j1+1;
       }
       else {
           ilmed=ilmed+i_value[j];
           j2=j2+1;
       }
       if(v_value[j]>vmed){
           vhmed=vhmed+v_value[j];
           j3=j3+1;
       }
       else {
           vlmed=vlmed+v_value[j];
           j4=j4+1;
       }
   }
   okFlag=true;
   if(j1==0 || j2==0 || j3==0 || j4==0) {okFlag=false;}
   if(okFlag==true){
       ihmed=ihmed/j1;
       ilmed=ilmed/j2;
       vhmed=vhmed/j3;
       vlmed=vlmed/j4;
       ECa=0.000;
       vh=1.8*(R3/(R3+R4))*((R1+R2)/R2)-((vhmed/1000.0)*(R1/R2));            
       vl=1.8*(R3/(R3+R4))*((R1+R2)/R2)-((vlmed/1000.0)*(R1/R2));      
       ih=1.8*(R3/(R3+R4))*((R1+R2)/R2)-((ihmed/1000.0)*(R1/R2));      
       il=1.8*(R3/(R3+R4))*((R1+R2)/R2)-((ilmed/1000.0)*(R1/R2));
       if((ih<0)&(il<0)){
         QMessageBox::information(this,"ECa Measurement","A problem has occurred. ih<0 and il<0 !!!");
         ECa=-9999.0;
        }
       else if((ih>0)&(il>0)){
             QMessageBox::information(this,"ECa Measurement","A problem has occurred. ih>0 and il>0 !!!");
             ECa=-9999.0;
             }
       else {
       ih=ih/R0;
       il=il/R0;
         if(ECa>=0.000){
           ECa=(il/(2*3.141592653589793238*a*vl));
           ECa=ECa+(ih/(2*3.141592653589793238*a*vh));
           ECa=1000.0*ECa/2.0;
           myMessage="ECa was measured = "+QString::number(ECa,'f',3)+" mS. You can save it !!!";
           QMessageBox::information(this,"ECa Measurement",myMessage);
         }
       }
       QFile file1("ECa_Signal.txt");
       file1.open(QFile::WriteOnly |QFile::Append | QFile::Text);
       QTextStream out1(&file1);
       out1<<"sample,j1,i_value[j1],v_value[j1],il,ih,vl,vh,signalF"<<endl;
	  for(int j1=0;j1<100;j1++){
          out1<<sample<<","<<j1<<","<<i_value[j1]<<","<<v_value[j1]<<","<<il<<","<<ih<<","<<vl<<","<<vh<<","<<signalF<<endl;
       }
       file1.flush();
       file1.close();
   }
   else {
        QMessageBox::information(this,"ECa Measurement","A problem has occurred. okFlag=false !!!");
   }
   return ECa;
}
//
int MySensor::readAnalog(int number)
{
    //reading analog ports of beaglebone
    stringstream ss;
    ss << LDR_PATH << number << "_raw";
    fstream fs;
    number=0;
    try{
        fs.open(ss.str().c_str(),fstream::in);
        fs >> number;
        fs.close();
    }
    catch(...)
    {
        QMessageBox::information(this, "Error","Error has occurred when reading ADC!!");
    }
    return number;
}
//
void MySensor::on_pushButton_clicked() // Save coordinates & ECa
{
        double ECa;
        QString myString;
        double dd;
        double dECadx;
        gps->GetGPSPosition(statusGPS, lat, longi, x, y, alt, hora, data, pdop, zone);
        xUTM=x;
        yUTM=y;
        CxPos[NC]=x;
        CyPos[NC]=y;
        CzPos[NC]=alt;
        ICPos[NC]='C';
        ECa=readingECa();
        if(ECa>=0.000){
	sample++;
	ui->lineEdit_8->setText(QString::number(sample));
	sumECa=sumECa+ECa;
        sumsqECa=sumsqECa+ECa*ECa;
        nECa=nECa+1;
        if (nECa>1) {
            varECa=(sumsqECa-(sumECa*sumECa)/nECa)/(nECa-1);
            if (varECa<0) varECa=0;
            cECa=(int) ((ECa-sumECa/nECa)/pow(varECa,0.5));
            ui->lineEdit_6->setText(QString::number(cECa));
            dd=pow((dxECa-xUTM)*(dxECa-xUTM)+(dyECa-yUTM)*(dyECa-yUTM),0.5);
            dECadx=(ECa-antECa)/dd;
            dsumECa=dsumECa+dECadx;
            dsumsqECa=dsumsqECa+dECadx*dECadx;
            dnECa=dnECa+1;
            if(dnECa>1){
                dvarECa=(dsumsqECa-(dsumECa*dsumECa)/dnECa)/(dnECa-1);
                if (dvarECa<0) dvarECa=0;
                dcECa=(int) ((dECadx-dsumECa/dnECa)/pow(dvarECa,0.5));
                ui->lineEdit_7->setText(QString::number(dcECa));
            }
        }
        dxECa=xUTM;
        dyECa=yUTM;
        antECa=ECa;
        myString=QString::number(typePoint)+",";
        myString=myString+QString::number(sample)+",";
        myString=myString+QString::number(data)+",";
        myString=myString+QString::number(hora)+",";
        myString=myString+QString::number(lat, 'f', 6)+",";
        myString=myString+QString::number(longi, 'f', 6)+",";
        myString=myString+QString::number(x, 'f', 3)+",";
        myString=myString+QString::number(y, 'f', 3)+",";
        myString=myString+QString::number(alt, 'f', 3)+",";
        myString=myString+QString::number(zone)+",";
        myString=myString+QString::number(pdop, 'f', 2)+",";
        myString=myString+ui->lineEdit_5->text()+",";
        myString=myString+QString::number(ECa, 'f', 3);
        ui->lineEdit_4->setText(QString::number(ECa, 'f', 3));
        already_clicked=true;
        NC=NC+1;
        NContour=NC;
        // Calculating the coordinate of the next point
        if(NC>1){
                x1=CxPos[NC-2];
                y1=CyPos[NC-2];
                x2=CxPos[NC-1];
                y2=CyPos[NC-1];
                dx=x2-x1;
                dy=y2-y1;
                d=pow(dx*dx+dy*dy,0.5);
                A[0][0]=dx/d;
                A[0][1]=-dy/d;
                A[1][0]=-A[0][1];
                A[1][1]=A[0][0];
                distFactor=1.0;
                if (ui->pushButton_2->isChecked()){
                    switch(abs(dcECa)) {
                        case 0:
                            break;
                        case 1:
                            distFactor=distFactor-0.10;
                            break;
                        case 2:
                            distFactor=distFactor-0.20;
                            break;
                        case 3:
                            distFactor=distFactor-0.30;
                            break;
                        default:
                            distFactor=distFactor-0.40;
                            break;
                    }
                    switch(abs(cECa)) {
                        case 0:
                            break;
                        case 1:
                            distFactor=distFactor-0.10;
                            break;
                        case 2:
                            distFactor=distFactor-0.20;
                            break;
                        case 3:
                            distFactor=distFactor-0.30;
                            break;
                        default:
                            distFactor=distFactor-0.40;
                            break;
                    }
                }
                y3=d+distFactor*reqDist;
                y3=0;
                newxUTM=x1+A[0][0]*x3+A[0][1]*y3;
                newyUTM=y1+A[1][0]*x3+A[1][1]*y3;
        }
        // Save the ECa readings and coordinates of the point in a file called contour.txt
        QFile file1("contour.txt");
        file1.open(QFile::WriteOnly | QFile::Append | QFile::Text);
        QTextStream out1(&file1);
        out1 << myString << endl;
        file1.flush();
        file1.close();
    // Increase 10Hz in the current frequency signal, up to 40Hz
	signalF=signalF+10;
	if(signalF>40){signalF=10;}
	ui->lineEdit_5->setText(QString::number(signalF));
	BlackPWM myPWM(P8_13);
	myPWM.setDutyPercent(100);
	myPWM.setPeriodTime(1000.0/signalF,milisecond);
	myPWM.setDutyPercent(50);
	sleep(1);
	}
        else {
             QMessageBox::information(this,"ECa Measurement","Measurement error. Please repeat ECa measurement !!!");
             return;
	}
}
void MySensor::on_pushButton_2_clicked() // Contour is done
{
    // Create the grid inside the polygon
    bool inside;
    double dx,dy;
    y=minyUTM;
    contourIsDone=true;
//
    if(contourIsDone==true){typePoint=0;}
//
    while(y<maxyUTM){
       x=minxUTM;
       while(x<maxxUTM){
          inside=pointInPolygon(x, y);
          if(inside==true){
             for(int kk=0;kk<NContour;kk++){
                dx=x-CxPos[kk];
                dy=y-CyPos[kk];
                d=pow((dx*dx+dy*dy),0.5);
                if(d<(0.5*reqDist)) {inside=false;}
             }
          }
//
          if(inside==true){
             GCxPos[NG]=x;
             GCyPos[NG]=y;
             NG=NG+1;
          }
          x=x+reqDist;
       }
       y=y+reqDist;
    }
    QMessageBox::information(this,"Contour Point","Grid points were generated successfully. Start with internal points reading !!!");
}
//
void MySensor::on_pushButton_4_clicked() // Data Acquisition Ended
{
    // Process the collected data by performing the semivariance analysis and the kriging interpolation
    QString MyFileName="contour.txt";
    startKrig->main_function(MyFileName, Zpixmap, Epixmap, zmin, zmax, ezmin, ezmax);
}
void MySensor::on_pushButton_5_clicked() // Show the Soil Apparente Electrical Conductivity Map
{
    // Show the ECa map
    scene->clear();
    sceneZmap->clear();
    ui->graphicsView->setScene(sceneZmap);
    sceneZmap->addPixmap(Zpixmap);
}
void MySensor::on_pushButton_6_clicked() // Show Error Map
{
    // Show the ECa error map
    scene->clear();
    sceneEmap->clear();
    ui->graphicsView->setScene(sceneEmap);
    sceneEmap->addPixmap(Epixmap);
}
bool MySensor::pointInPolygon(double xx, double yy) {
    // Verify if a point is inside a polygon (that defines the area contour)
    int   i;
    int j=NContour-1;
    bool  oddNodes=false;
    for (i=0; i<NContour; i++) {
        if ((((CyPos[i]< yy) && (CyPos[j]>=yy)) || (CyPos[j]< yy && CyPos[i]>=yy)) && (CxPos[i]<=xx || CxPos[j]<=xx)) {
            if (CxPos[i]+(yy-CyPos[i])/(CyPos[j]-CyPos[i])*(CxPos[j]-CxPos[i])<xx) {
                oddNodes=!oddNodes;
            }
        }
        j=i;
    }
  return oddNodes;
}