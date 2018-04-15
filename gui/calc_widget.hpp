#ifndef CALCWIDGET_HPP
#define CALCWIDGET_HPP

#include <iostream>
// #include <fstream>

#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QSpinBox>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QProgressBar>

#include "constants.hpp"
#include "calculation.hpp"

namespace GUI {

class CalcWidget : public QWidget {
    Q_OBJECT

    public:
        CalcWidget(const Arguments& args);

    public slots:
        void ChangeOutput(std::ostream* newOutStream);

    private:
        std::ostream* outStream;
        // number fields for N, L, P; I think that's the following 4 things:
        QVBoxLayout* layout;
        QFrame* inputBoxes;
        QSpinBox* nBox;
        QSpinBox* lBox;
        QSpinBox* pBox;
        QCheckBox* testCheckBox;
        QPushButton* goButton;
        QProgressBar* progressBar;

        void Calculate();
};

}

#endif
