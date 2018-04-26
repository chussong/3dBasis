#include "calc_widget.hpp"

namespace GUI {

CalcWidget::CalcWidget(const Arguments& args): console(args.console), 
        outStream(args.outStream), warningStatus(WARNING_ON), 
        nBox(new QSpinBox), lBox(new QSpinBox), dBox(new QDoubleSpinBox),
        pBox(new QSpinBox), msqBox(new QDoubleSpinBox), 
        lambdaBox(new QDoubleSpinBox), freeButton(new QRadioButton("&Free")),
        interactingButton(new QRadioButton("&Interacting")),
        testButton(new QRadioButton("&Tests only")), 
        goButton(new QPushButton("&go")), progressBar(new QProgressBar) {
    QVBoxLayout* layout = new QVBoxLayout;
    SetupInputBoxes(layout, args);
    SetupParameterBoxes(layout, args);
    SetupButtons(layout);
    layout->addWidget(progressBar);
    setLayout(layout);
}

void CalcWidget::SetupInputBoxes(QLayout* layout, const Arguments& args) {
    QFrame* inputBoxes = new QFrame;
    QGridLayout* inputLayout = new QGridLayout;

    QRadioButton* oneLevel = new QRadioButton(tr("&One (n,l) level"));
    inputLayout->addWidget(oneLevel, 0, 0);
    QFrame* oneLevelBoxes = SetupOneLevelFrame(args);
    inputLayout->addWidget(oneLevelBoxes, 0, 3, 1, 4);
    connect(oneLevel, &QAbstractButton::toggled, 
            oneLevelBoxes, &QWidget::setEnabled);

    QRadioButton* allLevels = new QRadioButton(tr("&All (n,l) levels"));
    inputLayout->addWidget(allLevels, 1, 0);
    QFrame* allLevelsBoxes = SetupAllLevelsFrame(args);
    inputLayout->addWidget(allLevelsBoxes, 1, 3, 1, 4);
    connect(allLevels, &QAbstractButton::toggled, 
            allLevelsBoxes, &QWidget::setEnabled);

    oneLevel->setChecked(true);
    allLevelsBoxes->setEnabled(false);

    inputBoxes->setLayout(inputLayout);
    layout->addWidget(inputBoxes);
}

QFrame* CalcWidget::SetupOneLevelFrame(const Arguments& args) {
    QFrame* frame = new QFrame;
    QHBoxLayout* layout = new QHBoxLayout;

    QLabel* nLabel = new QLabel("n");
    nLabel->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
    nLabel->setBuddy(nBox);
    nBox->setRange(2, 9);
    nBox->setValue(args.numP);
    nBox->setStatusTip(tr("Number of particles"));
    layout->addWidget(nLabel);
    layout->addWidget(nBox);

    QLabel* lLabel = new QLabel("l");
    lLabel->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
    lLabel->setBuddy(lBox);
    lBox->setRange(2, 10);
    lBox->setValue(args.degree);
    lBox->setSingleStep(2);
    lBox->setStatusTip(tr("Max number of derivatives above Dirichlet"));
    layout->addWidget(lLabel);
    layout->addWidget(lBox);

    frame->setLayout(layout);
    frame->setFrameStyle(QFrame::StyledPanel);
    return frame;
}

QFrame* CalcWidget::SetupAllLevelsFrame(const Arguments& args) {
    QFrame* frame = new QFrame;
    QHBoxLayout* layout = new QHBoxLayout;

    QLabel* dLabel = new QLabel(QString(0x0394));
    dLabel->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
    dLabel->setBuddy(dBox);
    layout->addWidget(dLabel);
    layout->addWidget(dBox);
    dBox->setRange(3.0, 20.0);
    dBox->setValue(args.delta);
    dBox->setSingleStep(0.5);
    dBox->setStatusTip(tr("Maximum total dimension"));

    frame->setLayout(layout);
    frame->setFrameStyle(QFrame::StyledPanel);
    return frame;
}

void CalcWidget::SetupParameterBoxes(QLayout* layout, const Arguments& args) {
    QFrame* paramBoxes = new QFrame;
    QHBoxLayout* paramBoxGrid = new QHBoxLayout;

    QLabel* msqLabel = new QLabel("m" + QString(0x00b2));
    msqLabel->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
    msqLabel->setBuddy(msqBox);
    msqBox->setRange(0.0, 1.0);
    msqBox->setValue(args.msq);
    msqBox->setSingleStep(0.05);
    msqBox->setStatusTip(tr("Coefficient of mass term in Hamiltonian"));

    paramBoxGrid->addWidget(msqLabel);
    paramBoxGrid->addWidget(msqBox);

    QLabel* lambdaLabel = new QLabel(QString(0x03bb));
    lambdaLabel->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
    lambdaLabel->setBuddy(lambdaBox);
    lambdaBox->setRange(0.0, 1.0);
    lambdaBox->setValue(args.msq);
    lambdaBox->setSingleStep(0.05);
    lambdaBox->setStatusTip(tr("Coefficient of interaction term in Hamiltonian"));

    paramBoxGrid->addWidget(lambdaLabel);
    paramBoxGrid->addWidget(lambdaBox);

    QLabel* pLabel = new QLabel("p");
    pLabel->setAlignment(Qt::AlignVCenter | Qt::AlignRight);
    pLabel->setBuddy(pBox);
    pBox->setRange(1, 100);
    pBox->setValue(args.partitions);
    pBox->setStatusTip(tr("Number of mu^2 partitions per operator"));

    paramBoxGrid->addWidget(pLabel);
    paramBoxGrid->addWidget(pBox);

    paramBoxes->setLayout(paramBoxGrid);
    layout->addWidget(paramBoxes);
}

void CalcWidget::SetupButtons(QLayout* layout) {
    connect(goButton, &QAbstractButton::clicked, this, &CalcWidget::Go);

    freeButton->setChecked(true);

    QButtonGroup* buttonGroup = new QButtonGroup;
    buttonGroup->addButton(freeButton);
    buttonGroup->addButton(interactingButton);
    buttonGroup->addButton(testButton);
    // layout->addWidget(buttonGroup);
    layout->addWidget(freeButton);
    layout->addWidget(interactingButton);
    layout->addWidget(testButton);
    layout->addWidget(goButton);
}

void CalcWidget::Go() {
    // if warning is both on and active, pop up a confirmation box
    if (warningStatus == (WARNING_ON | WARNING_ACTIVE)) {
        QMessageBox confirm;
        confirm.setText("This file has already been written to.");
        confirm.setInformativeText("Do you want to overwrite it?");
        confirm.setStandardButtons(QMessageBox::Yes | QMessageBox::Cancel);
        confirm.setDefaultButton(QMessageBox::Cancel);
        if (confirm.exec() == QMessageBox::Cancel) {
            return;
        } else {
            emit OverwriteFile();
        }
    }

    if (outStream != console) warningStatus |= WARNING_ACTIVE;
    goButton->setEnabled(false);

    emit StartingCalculation();
    QtConcurrent::run(this, &CalcWidget::Calculate);
}

void CalcWidget::Calculate() {
    Arguments args;
    args.numP = nBox->value();
    args.degree = lBox->value();
    args.delta = dBox->isEnabled() ? dBox->value() : 0.0;
    args.partitions = pBox->value();
    args.msq = msqBox->value();
    args.lambda = lambdaBox->value();
    args.outStream = outStream;
    args.console = console;

    if (outStream != console) args.options |= OPT_MATHEMATICA;
    if (interactingButton->isChecked()) {
        args.options |= OPT_INTERACTING;
    } else if (testButton->isChecked()) {
        args.options |= OPT_TEST;
    }

    ::Calculate(args);

    *console << "***** Calculation Complete *****" << endl;
    goButton->setEnabled(true);
}

void CalcWidget::ChangeOutput(QTextStream* newOutStream) {
    // console << "outStream updated" << std::endl;
    outStream = newOutStream;
    goButton->setEnabled(outStream != nullptr);
    warningStatus &= ~WARNING_ACTIVE;
}

void CalcWidget::GiveOverwriteWarnings(const bool newValue) {
    if (newValue == true) {
        warningStatus |= WARNING_ON;
    } else {
        warningStatus &= ~WARNING_ON;
    }
}

} // namespace GUI
