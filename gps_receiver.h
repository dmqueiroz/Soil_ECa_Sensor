#ifndef GPS_RECEIVER_H
#define GPS_RECEIVER_H
#include <QString>

/* baudrate settings are defined in <asm/termbits.h>, which is
   included by <termios.h> */
#define BAUDRATE B9600   // Change as needed, keep B

/* change this definition for the correct port */
#define MODEMDEVICE "/dev/ttyO1" //Beaglebone Black serial port

#define _POSIX_SOURCE 1 /* POSIX compliant source */

/*
#define FALSE 0
#define TRUE 1
*/

class GPS_Receiver
{
public:
    GPS_Receiver();
    void GetGPSPosition(QString &statusGPS, double &lat, double &longi, double &x, double &y, double &alt, double &hora, double &data, double &pdop, int &zone);
    void LatLongToUTM (double Lat, double Longi, double &x, double &y, int &zone);

};

#endif // GPS_RECEIVER_H
