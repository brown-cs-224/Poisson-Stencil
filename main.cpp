#include <QCoreApplication>
#include <QCommandLineParser>

#include <iostream>
#include <chrono>

#include "poisson.h"

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[])
{
    QCoreApplication a(argc, argv);

    ////////////////////////////////////////////////////////////////////////////////
    QCommandLineParser parser;
    parser.addHelpOption();
    parser.addPositionalArgument("infile", "Input .ply file path");
    parser.addPositionalArgument("outfile", "Output .obj file path");

    parser.process(a);

    const QStringList args = parser.positionalArguments();
    if(args.size() < 2) {
        cerr << "Error: Wrong number of arguments" << endl;
        a.exit(1);
        return 1;
    }
    QString infile = args[0];
    QString outfile = args[1];

    Poisson reconstructor;
    reconstructor.reconstruct(infile.toStdString(), outfile.toStdString());


    ////////////////////////////////////////////////////////////////////////////////



    a.exit();
}
