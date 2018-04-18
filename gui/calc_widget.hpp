#ifndef CALCWIDGET_HPP
#define CALCWIDGET_HPP

// #include <iostream>
// #include <fstream>

#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QDoubleSpinBox>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QProgressBar>
#include <QtWidgets/QMessageBox>
#include <QtConcurrent/QtConcurrentRun>
#include <QtCore/QTextStream>

#include "constants.hpp"
#include "calculation.hpp"

namespace GUI {

class CalcWidget : public QWidget {
    Q_OBJECT

    public:
        CalcWidget(const Arguments& args);

    signals:
        void StartingCalculation();
        void OverwriteFile();

    public slots:
        void Go();
        void ChangeOutput(QTextStream* newOutStream);
        void GiveOverwriteWarnings(const bool newValue);

    private:
        enum OverwriteWarning { WARNING_ON  = 0b10, WARNING_ACTIVE = 0b01 };
        void Calculate();

        QTextStream* console;
        QTextStream* outStream;
        int warningStatus;
        // number fields for N, L, P; I think that's the following 4 things:
        QSpinBox* nBox;
        QSpinBox* lBox;
        QSpinBox* pBox;
        QDoubleSpinBox* msqBox;
        QCheckBox* testCheckBox;
        QPushButton* goButton;
        QProgressBar* progressBar;
};

}

#endif
