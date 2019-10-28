#include "Kriging_Lib.h"
#include <QVector>
#include <QDebug>
#include <QFileDialog>
#include <QFile>
#include <QMessageBox>
#include <QStringList>
#include "math.h"
#include <iostream>
using namespace std;
Kriging_Lib::Kriging_Lib()
{
}
// main function that process the kriging interpolation and generates the map to be displayed in the graphical interface
void Kriging_Lib::main_function(QString &MyFileName, QPixmap &Zpixmap, QPixmap &Epixmap,double &zMin, double &zMax, double &ezMin, double &ezMax)
{
    QVector <int> it;
    QVector <int> pID;
    QVector <double> x;
    QVector <double> y;
    QVector <double> z;    
    it.clear();
    pID.clear();
    x.clear();
    y.clear();
    z.clear();
    bool inside[300][300];
    double zValues[300][300];
    double eValues[300][300];
    int polyCorners=0;
    int sWindow=300;
    int iaux;
    double dtopix;
    double xMax,xMin,yMax,yMin;    
    int neigID[16];
    double neigDist[16];
    double neigZ[16];
    int nNeig=16;
    int k;
    double maxD;
    double minD;
    double h[20];
    double gamma[20];
    QString modelo;
    double efeitoPepita;
    double alcance;
    double patamar;
    QString inLine;
    QStringList splittedLine;
    double sqdesvio;
    double semiModel[3][4];
    double xx,yy;
    double zEstimated,eEstimated;
    int myzcolor;
    double covMa[17][17];
    double xDca[17];
    double xDa[17];
    double readingECa;
    QFile file1;
    file1.setFileName(MyFileName);
    file1.open(QIODevice::ReadOnly | QIODevice::Text);
    QTextStream in(&file1);
    k=0;
    while (!in.atEnd()) {
      inLine = in.readLine();
      splittedLine=inLine.split(',');
      iaux=splittedLine.at(0).toInt();
      readingECa=splittedLine.at(11).toDouble();
      if(readingECa>0){
      	if (iaux==1) polyCorners++;
      	it << iaux;
      	pID << k;
        x << splittedLine.at(6).toDouble();
        y << splittedLine.at(7).toDouble();
        z << splittedLine.at(12).toDouble();
      	k++;
      }
    }
    file1.close();
    for(int j=0;j<sWindow;j++){
        for(int i=0;i<sWindow;i++){
            zValues[i][j]=0.0;
            eValues[i][j]=0.0;
        }
    }
    maxminDist(x,y,maxD,minD);
    semiVariogram(x,y,z,maxD,minD,h,gamma);
    modelo="exponencial";
    ajustaModelo(modelo, efeitoPepita, alcance, patamar, h, gamma);
    sqdesvio=somaQuadradoDesvio(modelo, efeitoPepita, alcance, patamar, h, gamma);
    semiModel[0][0]=sqdesvio;
    semiModel[0][1]=efeitoPepita;
    semiModel[0][2]=alcance;
    semiModel[0][3]=patamar;
    modelo="esferico";
    ajustaModelo(modelo, efeitoPepita, alcance, patamar, h, gamma);
    sqdesvio=somaQuadradoDesvio(modelo, efeitoPepita, alcance, patamar, h, gamma);
    semiModel[1][0]=sqdesvio;
    semiModel[1][1]=efeitoPepita;
    semiModel[1][2]=alcance;
    semiModel[1][3]=patamar;
    modelo="gaussiano";
    ajustaModelo(modelo, efeitoPepita, alcance, patamar, h, gamma);
    sqdesvio=somaQuadradoDesvio(modelo, efeitoPepita, alcance, patamar, h, gamma);
    semiModel[2][0]=sqdesvio;
    semiModel[2][1]=efeitoPepita;
    semiModel[2][2]=alcance;
    semiModel[2][3]=patamar;
    modelo="exponencial";
    sqdesvio=semiModel[0][0];
    efeitoPepita=semiModel[0][1];
    alcance=semiModel[0][2];
    patamar=semiModel[0][3];
    for (int i=1;i<3;i++){
        if (sqdesvio> semiModel[i][0]){
            sqdesvio=semiModel[i][0];
            efeitoPepita=semiModel[i][1];
            alcance=semiModel[i][2];
            patamar=semiModel[i][3];
            if (i==1) modelo="esferico";
            if (i==2) modelo="gaussiano";
        }
    }
    xMin=x[0];
    xMax=x[0];
    yMin=y[0];
    yMax=y[0];
    for(int i=1;i<x.size();i++){
        if (xMax<x[i]) xMax=x[i];
        if (xMin>x[i]) xMin=x[i];
        if (yMax<y[i]) yMax=y[i];
        if (yMin>y[i]) yMin=y[i];
    }
    if((xMax-xMin)>=(yMax-yMin)){
        dtopix=(xMax-xMin)/sWindow;
        yMin=(yMax+yMin)/2.0-sWindow*dtopix/2.0+dtopix/2;
        xMin=xMin+dtopix/2;
    }
    else {
        dtopix=(yMax-yMin)/sWindow;
        xMin=(xMax+xMin)/2.0-sWindow*dtopix/2.0+dtopix/2;
        yMin=yMin+dtopix/2;
    }
    creating_matriz_io(sWindow,polyCorners,inside,x,y,dtopix,xMin,yMin);
    zMax=0.0;
    zMin=1.0e30;
    ezMax=0.0;
    ezMin=1.0e30;
    for(int i=0;i<sWindow;i++){
        for(int j=0;j<sWindow;j++){
            if(inside[i][j]==true){
                yy=yMin+i*dtopix;
                xx=xMin+j*dtopix;
                find_PNeigborns(nNeig,xx,yy,pID,x,y,z,neigID,neigZ,neigDist);
                kriging(xx,yy,nNeig,modelo,efeitoPepita,alcance,patamar,neigID,neigZ,neigDist,x,y,zEstimated,eEstimated);
                zValues[i][j]=zEstimated;
                eValues[i][j]=eEstimated;
                if (zMax<zEstimated) zMax=zEstimated;
                if (zMin>zEstimated) zMin=zEstimated;
                if (ezMax<eEstimated) ezMax=eEstimated;
                if (ezMin>eEstimated) ezMin=eEstimated;
            }
         }
     }
     QImage Zimage(sWindow, sWindow, QImage::Format_RGB888);
     QImage Eimage(sWindow, sWindow, QImage::Format_RGB888);
     for(int i=0;i<sWindow;i++){         
         for(int j=0;j<sWindow;j++){             
             if(inside[i][j]==true){
                    myzcolor=(15*(zValues[i][j]-zMin)/(zMax-zMin));
                    switch (myzcolor){
                        case 14:
                            Zimage.setPixel(j, 300-i, qRgb(200, 0, 200));
                            break;
                        case 13:
                            Zimage.setPixel(j, 300-i, qRgb(200, 0, 0));
                            break;
                        case 12:
                            Zimage.setPixel(j, 300-i, qRgb(255, 0, 0));
                            break;
                        case 11:
                            Zimage.setPixel(j, 300-i, qRgb(255, 150, 0));
                            break;
                        case 10:
                            Zimage.setPixel(j, 300-i, qRgb(255, 200, 0));
                            break;
                        case 9:
                            Zimage.setPixel(j, 300-i, qRgb(255, 255, 0));
                            break;
                        case 8:
                            Zimage.setPixel(j, 300-i, qRgb(0, 100, 0));
                            break;
                        case 7:
                            Zimage.setPixel(j, 300-i, qRgb(0, 150, 0));
                            break;
                        case 6:
                            Zimage.setPixel(j, 300-i, qRgb(0, 200, 0));
                            break;
                        case 5:
                            Zimage.setPixel(j, 300-i, qRgb(0, 255, 50));
                            break;
                        case 4:
                            Zimage.setPixel(j, 300-i, qRgb(0, 255, 150));
                            break;
                        case 3:
                            Zimage.setPixel(j, 300-i, qRgb(0, 0, 150));
                            break;
                        case 2:
                            Zimage.setPixel(j, 300-i, qRgb(0, 0, 255));
                            break;
                        case 1:
                            Zimage.setPixel(j, 300-i, qRgb(0, 100, 255));
                            break;
                        case 0:
                            Zimage.setPixel(j, 300-i, qRgb(50, 175, 255));
                            break;
                    }
             } else{Zimage.setPixel(j, 300-i, qRgb(255, 255, 255));}
         }
     }
     for(int mzcolor=0;mzcolor<=14;mzcolor++){
        switch (mzcolor){
            case 14:
                for(int j=160;j<170;j++){
                    for(int i=280;i<290;i++){
                        Zimage.setPixel(j, i, qRgb(200, 0, 200));
                    }
                }
                break;
            case 13:
                for(int j=150;j<160;j++){
                    for(int i=280;i<290;i++){
                        Zimage.setPixel(j, i, qRgb(200, 0, 0));
                    }
                }
                break;
            case 12:
                for(int j=140;j<150;j++){
                    for(int i=280;i<290;i++){
                        Zimage.setPixel(j, i, qRgb(255, 0, 0));
                    }
                }
                break;
            case 11:
                for(int j=130;j<140;j++){
                    for(int i=280;i<290;i++){
                        Zimage.setPixel(j, i, qRgb(255, 150, 0));
                    }
                }
                break;
            case 10:
                for(int j=120;j<130;j++){
                    for(int i=280;i<290;i++){
                        Zimage.setPixel(j, i, qRgb(255, 200, 0));
                    }
                }
                break;
            case 9:
                for(int j=110;j<120;j++){
                    for(int i=280;i<290;i++){
                        Zimage.setPixel(j, i, qRgb(255, 255, 0));
                    }
                }
                break;
            case 8:
                for(int j=100;j<110;j++){
                    for(int i=280;i<290;i++){
                        Zimage.setPixel(j, i, qRgb(0, 100, 0));
                    }
                }
                break;
            case 7:
                for(int j=90;j<100;j++){
                    for(int i=280;i<290;i++){
                        Zimage.setPixel(j, i, qRgb(0, 150, 0));
                    }
                }
                break;
            case 6:
                for(int j=80;j<90;j++){
                    for(int i=280;i<290;i++){
                        Zimage.setPixel(j, i, qRgb(0, 200, 0));
                    }
                }
                break;
            case 5:
                for(int j=70;j<80;j++){
                    for(int i=280;i<290;i++){
                        Zimage.setPixel(j, i, qRgb(0, 255, 50));
                    }
                }
                break;
            case 4:
                for(int j=60;j<70;j++){
                    for(int i=280;i<290;i++){
                        Zimage.setPixel(j, i, qRgb(0, 255, 150));
                    }
                }
                break;
            case 3:
                for(int j=50;j<60;j++){
                    for(int i=280;i<290;i++){
                        Zimage.setPixel(j, i, qRgb(0, 0, 150));
                        }
                    }
                break;
            case 2:
                for(int j=40;j<50;j++){
                    for(int i=280;i<290;i++){
                        Zimage.setPixel(j, i, qRgb(0, 0, 255));
                    }
                }
                break;
            case 1:
                for(int j=30;j<40;j++){
                    for(int i=280;i<290;i++){
                        Zimage.setPixel(j, i, qRgb(0, 100, 255));
                    }
                }
                break;
            case 0:
                for(int j=20;j<30;j++){
                    for(int i=280;i<290;i++){
                        Zimage.setPixel(j, i, qRgb(50, 175, 255));
                    }
                }
                break;
        }
     }
     Zpixmap = QPixmap::fromImage(Zimage);
     for(int i=0;i<sWindow;i++){
         for(int j=0;j<sWindow;j++){
             if(inside[i][j]==true){
                    myzcolor=(15*(eValues[i][j]-ezMin)/(ezMax-ezMin));
                    switch (myzcolor){
                        case 14:
                            Eimage.setPixel(j, 300-i, qRgb(200, 0, 200));
                            break;
                        case 13:
                            Eimage.setPixel(j, 300-i, qRgb(200, 0, 0));
                            break;
                        case 12:
                            Eimage.setPixel(j, 300-i, qRgb(255, 0, 0));
                            break;
                        case 11:
                            Eimage.setPixel(j, 300-i, qRgb(255, 150, 0));
                            break;
                        case 10:
                            Eimage.setPixel(j, 300-i, qRgb(255, 200, 0));
                            break;
                        case 9:
                            Eimage.setPixel(j, 300-i, qRgb(255, 255, 0));
                            break;
                        case 8:
                            Eimage.setPixel(j, 300-i, qRgb(0, 100, 0));
                            break;
                        case 7:
                            Eimage.setPixel(j, 300-i, qRgb(0, 150, 0));
                            break;
                        case 6:
                            Eimage.setPixel(j, 300-i, qRgb(0, 200, 0));
                            break;
                        case 5:
                            Eimage.setPixel(j, 300-i, qRgb(0, 255, 50));
                            break;
                        case 4:
                            Eimage.setPixel(j, 300-i, qRgb(0, 255, 150));
                            break;
                        case 3:
                            Eimage.setPixel(j, 300-i, qRgb(0, 0, 150));
                            break;
                        case 2:
                            Eimage.setPixel(j, 300-i, qRgb(0, 0, 255));
                            break;
                        case 1:
                            Eimage.setPixel(j, 300-i, qRgb(0, 100, 255));
                            break;
                        case 0:
                            Eimage.setPixel(j, 300-i, qRgb(50, 175, 255));
                            break;
                    }
             } else{Eimage.setPixel(j, 300-i, qRgb(255, 255, 255));}
         }
     }
     for(int mzcolor=0;mzcolor<=14;mzcolor++){
        switch (mzcolor){
            case 14:
                for(int j=160;j<170;j++){
                    for(int i=280;i<290;i++){
                        Eimage.setPixel(j, i, qRgb(200, 0, 200));
                    }
                }
                break;
            case 13:
                for(int j=150;j<160;j++){
                    for(int i=280;i<290;i++){
                        Eimage.setPixel(j, i, qRgb(200, 0, 0));
                    }
                }
                break;
            case 12:
                for(int j=140;j<150;j++){
                    for(int i=280;i<290;i++){
                        Eimage.setPixel(j, i, qRgb(255, 0, 0));
                    }
                }
                break;
            case 11:
                for(int j=130;j<140;j++){
                    for(int i=280;i<290;i++){
                        Eimage.setPixel(j, i, qRgb(255, 150, 0));
                    }
                }
                break;
            case 10:
                for(int j=120;j<130;j++){
                    for(int i=280;i<290;i++){
                        Eimage.setPixel(j, i, qRgb(255, 200, 0));
                    }
                }
                break;
            case 9:
                for(int j=110;j<120;j++){
                    for(int i=280;i<290;i++){
                        Eimage.setPixel(j, i, qRgb(255, 255, 0));
                    }
                }
                break;
            case 8:
                for(int j=100;j<110;j++){
                    for(int i=280;i<290;i++){
                        Eimage.setPixel(j, i, qRgb(0, 100, 0));
                    }
                }
                break;
            case 7:
                for(int j=90;j<100;j++){
                    for(int i=280;i<290;i++){
                        Eimage.setPixel(j, i, qRgb(0, 150, 0));
                    }
                }
                break;
            case 6:
                for(int j=80;j<90;j++){
                    for(int i=280;i<290;i++){
                        Eimage.setPixel(j, i, qRgb(0, 200, 0));
                    }
                }
                break;
            case 5:
                for(int j=70;j<80;j++){
                    for(int i=280;i<290;i++){
                        Eimage.setPixel(j, i, qRgb(0, 255, 50));
                    }
                }
                break;
            case 4:
                for(int j=60;j<70;j++){
                    for(int i=280;i<290;i++){
                        Eimage.setPixel(j, i, qRgb(0, 255, 150));
                    }
                }
                break;
            case 3:
                for(int j=50;j<60;j++){
                    for(int i=280;i<290;i++){
                        Eimage.setPixel(j, i, qRgb(0, 0, 150));
                        }
                    }
                break;
            case 2:
                for(int j=40;j<50;j++){
                    for(int i=280;i<290;i++){
                        Eimage.setPixel(j, i, qRgb(0, 0, 255));
                    }
                }
                break;
            case 1:
                for(int j=30;j<40;j++){
                    for(int i=280;i<290;i++){
                        Eimage.setPixel(j, i, qRgb(0, 100, 255));
                    }
                }
                break;
            case 0:
                for(int j=20;j<30;j++){
                    for(int i=280;i<290;i++){
                        Eimage.setPixel(j, i, qRgb(50, 175, 255));
                    }
                }
                break;
        }
     }
     Epixmap = QPixmap::fromImage(Eimage);
}
// checking the points that inside the polygon to process the kriging interpolation
void Kriging_Lib::creating_matriz_io(int sWindow,int polyCorners,bool inside[300][300],const QVector <double> &coordX,const QVector <double> &coordY,double dtopix,double xMin,double yMin){
    double xx;
    double yy;
    int kk=0;
    QVector <double> polyX;
    QVector <double> polyY;
    for (kk=0;kk<polyCorners;kk++){
            polyX << coordX[kk];
            polyY << coordY[kk];
    }
    yy=yMin;
    for(int i=0;i<sWindow;i++){
        xx=xMin;
        for(int j=0;j<sWindow;j++){
            inside[i][j]=pointInPolygon(polyCorners,polyX,polyY,xx, yy);
            xx=xx+dtopix;
            kk=kk+1;
        }
        yy=yy+dtopix;
    }
}
// function to verify if a point is inside a polygon
bool Kriging_Lib::pointInPolygon(int polyCorners,const QVector <double> &polyX,const QVector <double> &polyY,double xx, double yy) {
    int   i;
    int j=polyCorners-1;
    bool  oddNodes=false;
    for (i=0; i<polyCorners; i++) {
        if ((((polyY[i]< yy) && (polyY[j]>=yy))  ||   (polyY[j]< yy && polyY[i]>=yy)) &&  (polyX[i]<=xx || polyX[j]<=xx)) {
            if (polyX[i]+(yy-polyY[i])/(polyY[j]-polyY[i])*(polyX[j]-polyX[i])<xx) {
                oddNodes=!oddNodes;
            }
        }
        j=i;
    }
  return oddNodes;
}
//
// Covariance function evaluation based on semivariogram parameters
double Kriging_Lib::covar_Funct(QString modelo,double efeitoPepita,double alcance, double patamar, double hh){
      double output=0.0;
      if (modelo=="exponencial"){
              if(hh>0.0){
                  output=(patamar-efeitoPepita)*exp(-3.0*hh/alcance);
              }
              else {
                  output=patamar;
              }
      }
      if (modelo=="gaussiano"){
              if(hh>0.0){
                  output=(patamar-efeitoPepita)*(1.0-exp(-3.0*hh*hh/(alcance*alcance)));
              }
              else {
                  output=patamar;
              }
      }
      if (modelo=="esferico"){
              if(hh>=alcance){
                  output=0.0;
              }
              else {
                  if(hh>0.0 && hh<alcance){
                      output=(patamar-efeitoPepita)*(3*hh/(2.0*alcance)-hh*hh*hh/(2.0*alcance*alcance*alcance));
                  }
                  else {
                      output=patamar;
                  }
              }
      }
      return output;
}
// Pivot function, to prevent zeros in the main diagonal, built base on code presented on 
// Chapra, Steven C., and Raymond P. Canale. 2014. Numerical methods for engineers. 
// 7nd Edition. Boston: McGraw-Hill Higher Education.
void Kriging_Lib::Pivot(double a[17][17], double b[17], double s[17], int k, int NROWS){
    int n=NROWS;
    int p=k;
    double big=fabs(a[k][k]/s[k]);
    double dummy;
    for(int ii=k+1;ii<n;ii++){
        dummy=fabs(a[ii][k]/s[ii]);
        if(dummy>big) {
            big=dummy;
            p=ii;
        }
    }
    if(p!=k){
        for(int jj=k;jj<n;jj++){
            dummy=a[p][jj];
            a[p][jj]=a[k][jj];
            a[k][jj]=dummy;
        }
        dummy=b[p];
        b[p]=b[k];
        b[k]=dummy;
        dummy=s[p];
        s[p]=s[k];
        s[k]=dummy;
     }
}
// Eliminate function, processing Gauss elimination, built base on code presented on 
// Chapra, Steven C., and Raymond P. Canale. 2014. Numerical methods for engineers. 
// 7nd Edition. Boston: McGraw-Hill Higher Education.
void Kriging_Lib::Eliminate(double a[17][17], double s[17], double b[17], double tol, bool &er, int NROWS){
    int n=NROWS;
    double factor;
    for(int k=0;k<n-1;k++){
        Pivot(a, b, s, k, n);
        if(fabs(a[k][k]/s[k])<tol){
                er=true;
                break;
        }
        for(int i=k+1;i<n;i++){
            factor=a[i][k]/a[k][k];
            for(int j=k+1;j<n;j++){
                a[i][j]=a[i][j]-factor*a[k][j];
            }
            b[i]=b[i]-factor*b[k];
        }
    }
    if(fabs(a[n-1][n-1]/s[n-1])<tol) er=true;
}
// Substitute function, to process the final step of the Gauss elimination, built base on code presented on 
// Chapra, Steven C., and Raymond P. Canale. 2014. Numerical methods for engineers. 
// 7nd Edition. Boston: McGraw-Hill Higher Education.
void Kriging_Lib::Substitute(double a[17][17],double b[17], double x[17], int NROWS){
    double sum;
    int n=NROWS;
    x[n-1]=b[n-1]/a[n-1][n-1];
    for(int i=n-2;i>=0;i--){
            sum=0.0;
            for (int j=i+1;j<n;j++){
                sum=sum+a[i][j]*x[j];
            }
            x[i]=(b[i]-sum)/a[i][i];
    }
}
// Gauss function, to solve a linear system of equations using Gauss elimination method, built base on code presented on 
// Chapra, Steven C., and Raymond P. Canale. 2014. Numerical methods for engineers. 
// 7nd Edition. Boston: McGraw-Hill Higher Education.
void Kriging_Lib::Gauss(double a[17][17], double b[17], double x[17], int NROWS){
    bool er=false;
    int n=NROWS;
    double tol=0.00001;
    double s[17];
    for(int i=0;i<n;i++){
        s[i]=fabs(a[i][0]);
        for(int j=1;j<n;j++){
            if(fabs(a[i][j])>s[i]) s[i]=fabs(a[i][j]);
        }
    }
    Eliminate(a,s,b,tol,er,n);
    if(er!=true) Substitute(a,b,x,n);
}
// Calculating the covariance matrix for the kriging interpolation
void Kriging_Lib::covarianceMatriz(int nNeig, QString modelo,double efeitoPepita,double alcance,double patamar,int neigID[],double neigZ[],double neigDist[],const QVector <double> &x,const QVector <double> &y, double covM[17][17],double xDc[17]){
    double xx, gamma;
    int ik,jk;
    for(int j=0;j<nNeig;j++){
        for(int i=j+1;i<nNeig;i++){
            ik=neigID[i];
            jk=neigID[j];
            xx=sqrt((x[ik]-x[jk])*(x[ik]-x[jk])+(y[ik]-y[jk])*(y[ik]-y[jk]));
            gamma=covar_Funct(modelo,efeitoPepita,alcance,patamar,xx);
            covM[i][j]=gamma;
            covM[j][i]=gamma;
        }
    }
    for(int i=0;i<nNeig;i++){
        covM[i][nNeig]=1.0;
        covM[nNeig][i]=1.0;
        xx=0;
        gamma=covar_Funct(modelo,efeitoPepita,alcance,patamar,xx);
        covM[i][i]=gamma;
        xx=neigDist[i];
        gamma=covar_Funct(modelo,efeitoPepita,alcance,patamar,xx);
        xDc[i]=gamma;
    }
    covM[nNeig][nNeig]=0.0;
    xDc[nNeig]=1.0;
}
// Running the kriging interpolation, built based on system of equations presented by
// Isaaks, E.H. & R.M. Srivastava. 1989. An introduction to applied geostatistics. Oxford University Press. New York.
void Kriging_Lib::kriging(double xx, double yy, int nNeig, QString modelo,double efeitoPepita,double alcance,double patamar,int neigID[],double neigZ[],double neigDist[],const QVector <double> &x,const QVector <double> &y,double &zEstimated,double &eEstimated){
    int n;
    double covM[17][17];
    double xD[17];
    double xDc[17];
    double b[17];
    int dummyID;
    double dummyZ;
    double dummyDist;
    int dummynNeig;
    bool negativeW;
    int k;
    negativeW=true;
    k=0;
    n=nNeig+1-k;
    while((negativeW==true) && (nNeig+1-k>=4)){
        n=nNeig+1-k;
        dummynNeig=nNeig-k;
        if (xx>1120 && xx<1130 && yy>2995 && yy<3000){
            for (int i=0;i<dummynNeig;i++){
            }
        }
        covarianceMatriz(dummynNeig,modelo,efeitoPepita,alcance,patamar,neigID,neigZ,neigDist,x,y,covM,xDc);
        for (int i=0;i<n;i++){
            b[i]=xDc[i];
        }
        Gauss(covM,b,xD,n);
        negativeW=false;
        for(int i=0;i<n-1;i++){
            if(xD[i]<0) {
                negativeW=true;
            }
        }
        if(negativeW==true) {
           k=k+4;
        }
    }
    zEstimated=0.0;
    eEstimated=patamar;
    for(int i=0;i<n-1;i++){
        zEstimated=zEstimated+xD[i]*neigZ[i];    //troquei i+1 por i
        eEstimated=eEstimated-xD[i]*xDc[i];
     }
     eEstimated=eEstimated-xD[n-1];
     if (eEstimated<0) {
         for(int i=0;i<n-1;i++){
         }
     eEstimated=0.0;
     }
     if(eEstimated>0.0){
        eEstimated=sqrt(eEstimated);
     }
}
// Function to find neighbors to process the kriging interpolation. 
void Kriging_Lib::find_PNeigborns(int nNeig,double xx,double yy,const QVector <int> &pID,const QVector <double> &x,const QVector <double> &y,const QVector <double> &z,int neigID[],double neigZ[],double neigDist[]){
    int k=x.size();
    double maxDist;
    double deltaX;
    double deltaY;
    double distancia;
    int maxJ;
    double minDist;
    int minNum;
    double minZ;
    int minID;
    for (int i=0;i<nNeig;i++){
           neigDist[i]=1E30;
    }
    for (int i=0;i<k;i++){
        maxDist=neigDist[0];
        maxJ=0;
        for(int j=1;j<nNeig;j++){
            if (maxDist<neigDist[j]){
                maxDist=neigDist[j];
                maxJ=j;
            }
        }
        deltaX=fabs(xx-x[i]);
        deltaY=fabs(yy-y[i]);
        distancia=sqrt(deltaX*deltaX+deltaY*deltaY);
        if(distancia<=maxDist){
            neigDist[maxJ]=distancia;
            neigID[maxJ]=pID[i];
            neigZ[maxJ]=z[i];
        }
    }
    for (int i=0;i<nNeig-1;i++){
        minDist=neigDist[i];
        minNum=i;
        minID=neigID[i];
        minZ=neigZ[i];
        for (int j=i+1;j<nNeig;j++){
            if(minDist>neigDist[j]){
                minDist=neigDist[j];
                minNum=j;
                minID=neigID[j];
                minZ=neigZ[j];
            }
        }
        if(i!=minNum){
            neigDist[minNum]=neigDist[i];
            neigID[minNum]=neigID[i];
            neigZ[minNum]=neigZ[i];
            neigDist[i]=minDist;
            neigID[i]=minID;
            neigZ[i]=minZ;
        }
    }
}
//    
// function to calcule the maximimun and minimum distances between points in the data set to define the semivariogram active lag distance   
void Kriging_Lib::maxminDist(const QVector <double> &x,const QVector <double> &y,double &maxD,double &minD){
    int k=x.size();
    double dist;
    maxD=sqrt((x[k-1]-x[k-2])*(x[k-1]-x[k-2])+(y[k-1]-y[k-2])*(y[k-1]-y[k-2]));
    minD=maxD;
    for(int i=0;i<k-1;i++){
        for (int j=i+1;j<k;j++){
           dist=sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]));
           if (dist>maxD) maxD=dist;
           if (dist<minD) minD=dist;
        }
    }
}
 //
 // calculating the experimental semivariogram with 20 data points (lags)
void Kriging_Lib::semiVariogram(const QVector <double> &x,const QVector <double> &y,const QVector <double> &z,double maxD,double minD,double h[], double gamma[]){
    double dist;
    double delta=(maxD-minD)*0.6/20;
    h[0]=minD+delta/2;
    int num[20];
    gamma[0]=0.0;
    num[0]=0;
    int k;
    int n=x.size();
    for (int i=1;i<20;i++) {
        h[i]=h[i-1]+delta;
        gamma[i]=0.0;
        num[i]=0;
    }
    for(int i=0;i<n-1;i++){
        for (int j=i+1;j<n;j++){
           dist=sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j]));
           k=0;
           while(dist>=(h[k]+delta/2)){
              k=k+1;
              if(k>=20) break;
           }
           if(k<20){
              gamma[k]=gamma[k]+(z[i]-z[j])*(z[i]-z[j]);
              num[k]=num[k]+1;
           }
        }
    }
    for (int i=0;i<20;i++){
        if(num[i]!=0){
            gamma[i]=gamma[i]/(2*num[i]);
        }
        else{
            gamma[i]=-99999.0;
        }
    }
}

// calculating the sum of squares error between the experimental and the theoretical semivariogram
double Kriging_Lib::somaQuadradoDesvio(QString modelo, double efeitoPepita, double alcance, double patamar, double h[], double gamma[]){
    int k=19;
    int imod;
    double somaQ=0.0;
    double diferenca=0.0;
    imod=3;
    if (modelo == "esferico") imod=1;
    if (modelo == "exponencial") imod=2;
    if (modelo == "gaussiano") imod=3;
//
    switch (imod){
        case 1:
            for(int i=0;i<k;i++){
                diferenca=0.0;
                if(h[i]<alcance){
                    if (gamma[i]>0) diferenca=gamma[i]-(efeitoPepita+(patamar-efeitoPepita)*(1.5*h[i]/alcance-0.5*h[i]*h[i]*h[i]/(alcance*alcance*alcance)));
                }
                else {
                    if (gamma[i]>0) diferenca=gamma[i]-patamar;
                }
                somaQ=somaQ+diferenca*diferenca;
            }
            break;
        case 2:
            for(int i=0;i<k;i++){
                diferenca=0.0;
                if (gamma[i]>0) diferenca=gamma[i]-(efeitoPepita+(patamar-efeitoPepita)*(1-exp(-3.0*h[i]/alcance)));
                somaQ=somaQ+diferenca*diferenca;
            }
            break;
        case 3:
            for(int i=0;i<k;i++){
                diferenca=0.0;
                if (gamma[i]>0) diferenca=gamma[i]-(efeitoPepita+(patamar-efeitoPepita)*(1-exp(-3.0*h[i]*h[i]/(alcance*alcance))));
                somaQ=somaQ+diferenca*diferenca;
            }
            break;
    }
    return somaQ;
}
// using the golden rule to find the semivariogram parameters that minimize the sum of squares error. Developed based on code presented by
// Chapra, Steven C., and Raymond P. Canale. 2014. Numerical methods for engineers. 
// 7nd Edition. Boston: McGraw-Hill Higher Education.
double Kriging_Lib::gold(QString problem,char ivar,double xlow, double xhigh, QString modelo, double x, double y, double z, int maxIt, double es, double h[], double gamma[]){
    int fp,iter;
    double r, xl, xu, d, x1, x2, f1, f2, ea;
    double xopt, fx, xint;
    if (problem.toLower()=="minimizar"){
        fp = -1;
    }
    else {
        fp = 1;
    }
   //
    r = 0.618033989;
    xl = xlow;
    xu = xhigh;
    iter = 1;
    d = r*(xu-xl);
    //
    x1 = xl+d;
    x2 = xu-d;
    fx=0.0;
    switch (ivar){
        case 'x':
            f1 = somaQuadradoDesvio(modelo, x1, y, z, h, gamma);
            f2 = somaQuadradoDesvio(modelo, x2, y, z, h, gamma);
            break;
        case 'y':
            f1 = somaQuadradoDesvio(modelo, x, x1, z, h, gamma);
            f2 = somaQuadradoDesvio(modelo, x, x2, z, h, gamma);
            break;
        case 'z':
            f1 = somaQuadradoDesvio(modelo, x, y, x1, h, gamma);
            f2 = somaQuadradoDesvio(modelo, x, y, x2, h, gamma);
            break;
    }
    if (f1*fp > f2*fp) {
        xopt = x1;
        fx = f1;
    }
    else{
        xopt = x2;
        fx = f2;
    }
    //
    while (true){
        d = r*d;
        xint = xu-xl;
        //
        if (f1*fp > f2*fp) {
            xl = x2;
            x2 = x1;
            x1 = xl+d;
            f2 = f1;
            switch (ivar){
                case 'x':
                    f1 = somaQuadradoDesvio(modelo, x1, y, z, h, gamma);
                    break;
                case 'y':
                    f1 = somaQuadradoDesvio(modelo, x, x1, z, h, gamma);
                    break;
                case 'z':
                    f1 = somaQuadradoDesvio(modelo, x, y, x1, h, gamma);
                    break;
            }
        }
        else{
            xu = x1;
            x1 = x2;
            x2 = xu-d;
            f1 = f2;
            switch (ivar){
                case 'x':
                    f2 = somaQuadradoDesvio(modelo, x2, y, z, h, gamma);
                    break;
                case 'y':
                    f2 = somaQuadradoDesvio(modelo, x, x2, z, h, gamma);
                    break;
                case 'z':
                    f2 = somaQuadradoDesvio(modelo, x, y, x2, h, gamma);
                    break;
            }
        }
        iter = iter+1;
        //
        if (f1*fp > f2*fp){
            xopt = x1;
            fx = f1;
        }
        else{
            xopt = x2;
            fx = f2;
        }
        //
        if (xopt != 0.0) ea = (1-r)*fabs(xint/xopt)*100;
        if (ea <= es or iter >= maxIt) break;
    }
    return xopt;
}
//
// calculating the parameters of the theoretical semivariogram
void Kriging_Lib::ajustaModelo(QString &modelo, double &efeitoPepita, double &alcance, double &patamar, double h[], double gamma[]){
        QString problem="minimizar";
        double fant,error,fxyz;
        double xL, xU;
        int j,imaxit;
        int maxIt=25;
        char ivar;
        double es=0.01;
        efeitoPepita=(gamma[1]*h[0]-gamma[0]*h[1])/(h[0]-h[1]);     // the first guess for the semivariogram nugget effect
        patamar=(gamma[16]+gamma[17]+gamma[18]+gamma[19])/4.0;      // the first guess for the semivariogram sill
        alcance=h[10];                                              // the first guess for the semivariogram range
        fant=somaQuadradoDesvio(modelo, efeitoPepita, alcance, patamar, h, gamma); 
        imaxit = 25;                                                
        j = 1;                                                      
        //
        // the main loop to find the nugget effect, the sill an the range of the semivariogram
        while (true){
            //
            // first find the optimum value of the nugget efect
            ivar = 'x';
            xL=0.00001;
            xU= (gamma[16]+gamma[17]+gamma[18]+gamma[19])/4.0;
            efeitoPepita=gold(problem, ivar, xL, xU, modelo, efeitoPepita, alcance, patamar, maxIt, es, h, gamma);
            //
            //
            // then find the best value for the range
            ivar = 'y';
            xL=0.0001;
            xU=h[19];
            alcance=gold(problem, ivar, xL, xU, modelo, efeitoPepita, alcance, patamar, maxIt, es, h, gamma);
            //
            //
            // then find the best value for the sill
            ivar = 'z';
            xL=(gamma[0]+gamma[1])/2;
            xU=1.5*(gamma[16]+gamma[17]+gamma[18]+gamma[19])/4.0;
            patamar=gold(problem, ivar, xL, xU, modelo, efeitoPepita, alcance, patamar, maxIt, es, h, gamma);
            //
            if (fant != 0) {
                j = j + 1;
                //
                // calculate the sum of squares error
                fxyz=somaQuadradoDesvio(modelo, efeitoPepita, alcance, patamar, h, gamma);
                error = 100 * fabs((fant - fxyz) / fant);
                fant = fxyz;
                //
                // check if the solution is OK
                if ((j >= imaxit) or (error < es)) break;
            }
            else {
            }
         }
}