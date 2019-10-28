#include <gps_receiver.h>
#include <stdlib.h>
#include <QDebug>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <termios.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <QString>
#include <QStringList>
#include <math.h>
#include <cmath>


GPS_Receiver::GPS_Receiver()
{
}

void GPS_Receiver::GetGPSPosition(QString &statusGPS, double &lat, double &longi, double &x, double &y, double &alt, double &hora, double &data, double &pdop, int &zone)
{
    int fd;
    statusGPS="Void";   
    int res;
    int ret;
    int kk;
    bool gpsS;
    bool gpsFlag;
    QString myData;
    QStringList splittedData1;
    QStringList splittedData2;
    QString minutesdecimal;    
    double degrees;
    double minutes;
    struct termios oldtio, newtio;    
    char buf[255];
    ret=0;
        // Load the pin configuration
    ret = system("echo UART1 > /sys/devices/bone_capemgr.9/slots");
        // Open modem device for reading and writing and not as controlling tty
        // because we don't want to get killed if linenoise sends CTRL-C.
    fd = open(MODEMDEVICE, O_RDWR | O_NOCTTY );
    if (fd < 0) { perror(MODEMDEVICE); exit(-1); }

    bzero(&newtio, sizeof(newtio)); // clear struct for new port settings

        //   BAUDRATE: Set bps rate. You could also use cfsetispeed and cfsetospeed.
        //   CRTSCTS : output hardware flow control (only used if the cable has
        //             all necessary lines. See sect. 7 of Serial-HOWTO)
        //   CS8     : 8n1 (8bit,no parity,1 stopbit)
        //   CLOCAL  : local connection, no modem contol
        //   CREAD   : enable receiving characters
    newtio.c_cflag = BAUDRATE | CRTSCTS | CS8 | CLOCAL | CREAD;

        // IGNPAR  : ignore bytes with parity errors
        //   otherwise make device raw (no other input processing)
    newtio.c_iflag = IGNPAR;

        //  Raw output
    newtio.c_oflag = 0;

        // ICANON  : enable canonical input
        // disable all echo functionality, and don't send signals to calling program
    newtio.c_lflag = ICANON;
    // now clean the modem line and activate the settings for the port
    tcflush(fd, TCIFLUSH);
    tcsetattr(fd,TCSANOW,&newtio);
    // NMEA command to ouput all sentences
    // Note that this code & format values in manual are hexadecimal
    write(fd, "$PTNLSNM,273F,01*27\r\n", 21);
    // terminal settings done, now handle input
    kk=0;
    gpsFlag=true;
    lat=0.0;
    longi=0.0;
    while (gpsFlag)
    {   // loop continuously
        // read blocks program execution until a line terminating character is
        // input, even if more than 255 chars are input. If the number
        // of characters read is smaller than the number of chars available,
        // subsequent reads will return the remaining chars. res will be set
        // to the actual number of characters actually read
        res = read(fd, buf, 255);
        buf[res] = 0;             
        if(kk>6)
        {
            myData="";
            for (int ii=0;ii<=res;ii++)
            {
                myData=myData+buf[ii];
            }

            splittedData1 = myData.split(',');
            if(splittedData1.at(0)=="$GPRMC")
            {
                if(splittedData1.at(2)=="A")
                {
                    hora=splittedData1.at(1).toDouble();
                    gpsS=true;
                    splittedData2=splittedData1.at(3).split('.');
                    degrees=(double)(splittedData2.at(0).toInt()/100);
                    minutes=(double)(splittedData2.at(0).toInt()%100);
                    minutesdecimal="0."+splittedData2.at(1);
                    minutes=minutes+minutesdecimal.toDouble();
                    lat=degrees+minutes/60;
                    if(splittedData1.at(4)=="S"){lat=-lat;}
                    splittedData2=splittedData1.at(5).split('.');
                    degrees=(double)(splittedData2.at(0).toInt()/100);
                    minutes=(double)(splittedData2.at(0).toInt()%100);
                    minutesdecimal="0."+splittedData2.at(1);
                    minutes=minutes+minutesdecimal.toDouble();
                    longi=degrees+minutes/60;
                    if(splittedData1.at(6)=="W"){longi=-longi;}
                    data=splittedData1.at(9).toDouble();
                    statusGPS="Active";
                }
                else
                {
                    gpsS=false;
                }
            }
            if(splittedData1.at(0)=="$GPGSA")
            {
                if(splittedData1.at(1)=="A")
                {
                    pdop=splittedData1.at(15).toDouble();
                    statusGPS="Active";
                }
                else
                {
                    gpsS=false;
                }
            }
            if(splittedData1.at(0)=="$GPGGA")
            {
                if(splittedData1.at(6)=="1")
                {
                    alt=splittedData1.at(9).toDouble();
                    statusGPS="Active";

                }
                else
                {
                    gpsS=false;
                }
            }
        }
        kk=kk+1;
        if(kk>16){gpsFlag=false;};
    }
    tcsetattr(fd, TCSANOW, &oldtio);
    LatLongToUTM(lat,longi, x, y, zone);
}
//
// Converting geographical coordinates to UTM
void GPS_Receiver::LatLongToUTM (double Lat, double Longi, double &x, double &y, int &zone)
{
    double a=6378137.0;                  // equatorial radius
    double b=6356752.314;                // polar radius
    double k0=0.9996;                    // scale factor
    double e=pow((1.0-(b*b)/(a*a)),0.5); // eccentricity    
    double n=(a-b)/(a+b);                // 3d flattening   
    double n2=n*n;
    double n3=n2*n;
    double n4=n2*n2;
    double n5=n4*n;
    double n6=n4*n2;
    double n7=n6*n;
    double n8=n4*n4;
    double n9=n8*n;
    double n10=n8*n2;

    double AA=(a/( 1.0 + n ) )*(1.0 + (1.0/4.0)* n2 + (1.0/64.0)*n4 + (1.0/256.0)*n6 + (25.0/16384.0)*n8 + (49.0/65536.0)*n10);  // Meridian Radius

    // Kruger Series
    double alpha1 = (1.0/2.0)*n - (2.0/3.0)*n2 + (5.0/16.0)*n3 + (41.0/180.0)*n4 - (127.0/288.0)*n5;
    alpha1 = alpha1 + (7891.0/37800.0)*n6 + (72161.0/387072.0)*n7 - (18975107.0/50803200.0)*n8;
    alpha1 = alpha1 + (60193001.0/290304000.0)*n9 + (134592031.0/1026432000.0)*n10;

    double alpha2 = (13.0/48.0)*n2 - (3.0/5.0)*n3 + (557.0/1440.0)*n4 + (281.0/630.0)*n5 - (1983433.0/1935360.0)*n6;
    alpha2 = alpha2 + (13769.0/28800.0)*n7 + (148003883.0/174182400.0)*n8 - (705286231.0/465696000.0)*n9;
    alpha2 = alpha2 + (1703267974087.0/3218890752000.0)*n10;

    double alpha3 = (61.0/240.0)*n3 - (103.0/140.0)*n4 + (15061.0/26880.0)*n5 + (167603.0/181440.0)*n6;
    alpha3 = alpha3 - (67102379.0/29030400.0)*n7 + (79682431.0/79833600.0)*n8 + (6304945039.0/2128896000.0)*n9;
    alpha3 = alpha3 -  (6601904925257.0/1307674368000.0)*n10;

    double alpha4 = (49561.0/161280.0)*n4 - (179.0/168.0)*n5 + (6601661.0/7257600.0)*n6;
    alpha4 = alpha4 + (97445.0/49896.0)*n7 - (40176129013.0/7664025600.0)*n8;
    alpha4= alpha4 + (138471097.0/66528000.0)*n9 + (48087451385201.0/5230697472000.0)*n10;

    double alpha5 = (34729.0/80640.0)*n5 - (3418889.0/1995840.0)*n6 + (14644087.0/9123840.0)*n7;
    alpha5 = alpha5 +   (2605413599.0/622702080.0)*n8 - (31015475399.0/2583060480.0)*n9;
    alpha5 = alpha5 +  (5820486440369.0/1307674368000.0)*n10;

    double alpha6 = (212378941.0/319334400.0)*n6 - (30705481.0/10378368.0)*n7 + (175214326799.0/58118860800.0)*n8;
    alpha6 = alpha6 + (870492877.0/96096000.0)*n9 - (1328004581729000.0/47823519744000.0)*n10;

    double alpha7 = (1522256789.0/1383782400.0)*n7 - (16759934899.0/3113510400.0)*n8 + (1315149374443.0/221405184000.0)*n9;
    alpha7 = alpha7 + (71809987837451.0/3629463552000.0)*n10;

    double alpha8 = (1424729850961.0/743921418240.0)*n8 -   (256783708069.0/25204608000.0)*n9;
    alpha8 = alpha8 + (2468749292989890.0/203249958912000.0)*n10;

    double beta1 = (1.0/2.0)*n-(2.0/3.0)*n2 + (37.0/96.0)*n3 - (1.0/360.0)*n4 - (81.0/512.0)*n5 + (96199.0/604800.0)*n6;
    beta1 = beta1 - (5406467.0/38707200.0)*n7 + (7944359.0/67737600.0)*n8 - (7378753979.0/97542144000.0)*n9;
    beta1 = beta1 + (25123531261.0/804722688000.0)*n10;

    double beta2 = (1.0/48.0)*n2 + (1.0/15.0)*n3 - (437.0/1440.0)*n4 + (46.0/105.0)*n5 - (1118711.0/3870720.0)*n6;
    beta2 = beta2 + (51841.0/1209600.0)*n7 + (24749483.0/348364800.0)*n8 - (115295683.0/1397088000.0)*n9;
    beta2 = beta2 + (5487737251099.0/51502252032000.0)*n10;

    double beta3 = (17.0/480.0)*n3 - (37.0/840.0)*n4 - (209.0/4480.0)*n5 + (5569.0/90720.0)*n6 + (9261899.0/58060800.0)*n7;
    beta3 = beta3 - (6457463.0/17740800.0)*n8 + (2473691167.0/9289728000.0)*n9 - (852549456029.0/20922789888000.0)*n10;

    double beta4 = (4397.0/161280.0)*n4 - (11.0/504.0)*n5 - (830251.0/7257600.0)*n6 + (466511.0/2494800.0)*n7;
    beta4 = beta4 + (324154477.0/7664025600.0)*n8 - (937932223.0/3891888000.0)*n9 - (89112264211.0/5230697472000.0)*n10;

    double beta5 = (4583.0/161280.0)*n5 - (108847.0/3991680.0)*n6 - (8005831.0/63866880.0)*n7 + (22894433.0/124540416.0)*n8;
    beta5 = beta5 + (112731569449.0/557941063680.0)*n9 - (5391039814733.0/10461394944000.0)*n10;

    double beta6 = (20648693.0/638668800.0)*n6 -  (16363163.0/518918400.0)*n7 - (2204645983.0/12915302400.0)*n8;
    beta6 = beta6 + (4543317553.0/18162144000.0)*n9 + (54894890298749.0/167382319104000.0)*n10;

    double beta7 = (219941297.0/5535129600.0)*n7 - (497323811.0/12454041600.0)*n8 - (79431132943.0/332107776000.0)*n9;
    beta7 = beta7 + (4346429528407.0/12703122432000.0)*n10;

    double beta8 = (191773887257.0/3719607091200.0)*n8 -  (17822319343.0/336825216000.0)*n9;
    beta8 = beta8 - (497155444501631.0/1422749712384000.0)*n10;

    double abslat = fabs(Lat);

    double latr = abslat*3.14159265358979/180.0;

    int zone1;
    if(Longi<0) {
            zone1 = 1+(int)((180+Longi)/6);
    }
    else{
        zone1 = 31+(int)(Longi/6);
    }
    int zoneCM = 6*zone1-183;
    int Meridian = zoneCM;
    int EastCM = 1;
    if(Longi<Meridian){EastCM=-1;}
    double DlongD = fabs(Longi-Meridian);
    double dlong = DlongD;
    double dlongr = dlong*3.14159265358979/180;
    double conflat = atan(sinh(asinh(tan(latr))-e*atanh(e*sin(latr))));
    double sigma = sinh(e*atanh(e*tan(latr)/pow((1.0+(tan(latr)*tan(latr))),0.5)));
    conflat = atan(tan(latr)*pow((1.0+sigma*sigma),0.5)-sigma*pow((1+tan(latr)*tan(latr)),0.5));
    double tauprime = tan(conflat);
    double xiprime = atan(tauprime/cos(dlongr));
    double etaprime = asinh(sin(dlongr)/pow((tauprime*tauprime+(cos(dlongr)*cos(dlongr))),0.5));
    double xi=xiprime+alpha1*sin(2*xiprime)*cosh(2*etaprime)+alpha2*sin(4*xiprime)*cosh(4*etaprime);
    xi = xi+alpha3*sin(6*xiprime)*cosh(6*etaprime)+alpha4*sin(8*xiprime)*cosh(8*etaprime);
    xi = xi+alpha5*sin(10*xiprime)*cosh(10*etaprime)+alpha6*sin(12*xiprime)*cosh(12*etaprime);
    xi = xi+alpha7*sin(14*xiprime)*cosh(14*etaprime);

    double eta = etaprime+alpha1*cos(2*xiprime)*sinh(2*etaprime)+alpha2*cos(4*xiprime)*sinh(4*etaprime);
    eta = eta+alpha3*cos(6*xiprime)*sinh(6*etaprime)+alpha4*cos(8*xiprime)*sinh(8*etaprime);
    eta = eta+alpha5*cos(10*xiprime)*sinh(10*etaprime)+alpha6*cos(12*xiprime)*sinh(12*etaprime);
    eta = eta+alpha6*cos(14*xiprime)*sinh(14*etaprime);

    double easting = k0*AA*eta;
    double northing = k0*AA*xi;

    if(Lat<0){northing=10000000-northing;}
    easting=500000.0+EastCM*easting;

    x=easting;
    y=northing;
    zone=zone1;
}

