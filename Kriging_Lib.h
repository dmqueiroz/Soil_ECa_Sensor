#ifndef KRIGING_LIB_H
#define KRIGING_LIB_H

#include <QDialog>
#include <QGraphicsScene>
#include <QVector>
#include <QFile>
#include <QString>

class Kriging_Lib
{
public:
    Kriging_Lib();
    void main_function(QString &MyFileName, QPixmap &Zpixmap, QPixmap &Epixmap, double &zMin, double &zMax, double &ezMin, double &ezMax);
    void maxminDist(const QVector <double> &x, const QVector <double> &y, double &maxD, double &minD);
    void semiVariogram(const QVector <double> &x, const QVector <double> &y, const QVector <double> &z, double maxD, double minD, double h[], double gamma[]);
    void ajustaModelo(QString &modelo, double &efeitoPepita, double &alcance, double &patamar, double h[], double gamma[]);
    double somaQuadradoDesvio(QString modelo, double efeitoPepita, double alcance, double patamar, double h[], double gamma[]);
    double gold(QString problem,char ivar,double xlow, double xhigh, QString modelo, double x, double y, double z,int maxIt, double es, double h[], double gamma[]);

    void creating_matriz_io(int sWindow,int polyCorners,bool inside[300][300],const QVector <double> &coordX,const QVector <double> &coordY,double dtopix,double xMin,double yMin);
    bool pointInPolygon(int polyCorners,const QVector <double> &polyX,const QVector <double> &polyY,double xx, double yy);
    void find_PNeigborns(int nNeig,double xx,double yy,const QVector <int> &pID,const QVector <double> &x,const QVector <double> &y,const QVector <double> &z,int neigID[],double neigZ[],double neigDist[]);
    void kriging(double xx, double yy, int nNeig,QString modelo,double efeitoPepita,double alcance,double patamar,int neigID[], double neigZ[],double neigDist[],const QVector <double> &x,const QVector <double> &y,double &zEstimated,double &eEstimated);
    double covar_Funct(QString modelo,double efeitoPepita,double alcance, double patamar, double hh);
    void Gauss(double a[17][17], double b[17], double x[17], int n);
    void Substitute(double a[17][17],double b[17], double x[17], int n);
    void Eliminate(double a[17][17], double s[17], double b[17], double tol, bool &er, int n);
    void Pivot(double a[17][17], double b[17], double s[17], int k, int n);
    void covarianceMatriz(int nNeig, QString modelo,double efeitoPepita,double alcance,double patamar,int neigID[],double neigZ[],double neigDist[],const QVector <double> &x,const QVector <double> &y, double covM[17][17],double xDc[17]);

};

#endif // KRIGING_LIB_H
