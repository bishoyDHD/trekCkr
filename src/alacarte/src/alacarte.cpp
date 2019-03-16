#include <QApplication>
#include "mainwindow.h"

int main (int argc, char ** argv)
{

  QApplication app(argc,argv);
  app.setOrganizationName(.cooker");
  app.setApplicationName("A la carte");

  app.setStyleSheet("QLineEdit::disabled {background: darkgray;}");

  MainWindow mainWin;
  mainWin.show();
  return app.exec();
}

