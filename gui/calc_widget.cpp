#include "calc_widget.hpp"

namespace GUI {

CalcWidget::CalcWidget(const Arguments& args): console(args.console), 
        outStream(args.outStream), warningStatus(WARNING_ON), 
        nBox(new QSpinBox), lBox(new QSpinBox), pBox(new QSpinBox), 
        msqBox(new QDoubleSpinBox), 
        testCheckBox(new QCheckBox("Run &tests only")), 
        goButton(new QPushButton("&go")), progressBar(new QProgressBar) {
    // put the boxes into inputBoxGrid
    QFrame* inputBoxes = new QFrame;
    QHBoxLayout* inputBoxGrid = new QHBoxLayout;
    inputBoxGrid->addWidget(new QLabel("n"));
    inputBoxGrid->addWidget(nBox);
    inputBoxGrid->addWidget(new QLabel("l"));
    inputBoxGrid->addWidget(lBox);
    inputBoxGrid->addWidget(new QLabel("p"));
    inputBoxGrid->addWidget(pBox);
    inputBoxGrid->addWidget(new QLabel("m^2")); // FIXME: unicode exponent
    inputBoxGrid->addWidget(msqBox);
    inputBoxes->setLayout(inputBoxGrid);

    nBox->setRange(2, 9);
    nBox->setValue(args.numP);
    nBox->setStatusTip(tr("Number of particles"));

    lBox->setRange(2, 10);
    lBox->setValue(args.degree);
    lBox->setSingleStep(2);
    lBox->setStatusTip(tr("Max number of derivatives above Dirichlet"));

    pBox->setRange(1, 100);
    pBox->setValue(args.partitions);
    pBox->setStatusTip(tr("Number of mu^2 partitions per operator"));

    msqBox->setRange(0.0, 1.0);
    msqBox->setValue(args.msq);
    msqBox->setSingleStep(0.05);
    msqBox->setStatusTip(tr("Coefficient of mass operator in Hamiltonian"));

    connect(goButton, &QAbstractButton::clicked, this, &CalcWidget::Go);

    QVBoxLayout* layout = new QVBoxLayout;
    layout->addWidget(inputBoxes);
    layout->addWidget(testCheckBox);
    layout->addWidget(goButton);
    layout->addWidget(progressBar);
    setLayout(layout);
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
    args.partitions = pBox->value();
    args.msq = msqBox->value();
    args.outStream = outStream;
    args.console = console;

    if (outStream != console) args.options |= OPT_MATHEMATICA;
    if (testCheckBox->isChecked()) args.options |= OPT_TEST;

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
