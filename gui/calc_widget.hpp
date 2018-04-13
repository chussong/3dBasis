#ifndef CALCWIDGET_HPP
#define CALCWIDGET_HPP

#include <iostream>

#include <QtWidgets/QCheckBox>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QProgressBar>

#include "constants.hpp"
#include "calculation.hpp"

namespace GUI {

class CalcWidget : public QVBoxLayout {
    Q_OBJECT

    public:
        CalcWidget(const Arguments& args);

    private:
        std::ostream& outStream;
        // number fields for N, L, P; I think that's the following 4 things:
        QHBoxLayout* inputBoxGrid;
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
