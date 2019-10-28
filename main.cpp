#include "mysensor.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    MySensor w;
    w.show();

    return a.exec();
}
