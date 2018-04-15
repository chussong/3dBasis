#include "calc_widget.hpp"

namespace GUI {

CalcWidget::CalcWidget(const Arguments& args): outStream(args.outputStream),
        layout(new QVBoxLayout), inputBoxes(new QFrame), 
        nBox(new QSpinBox), lBox(new QSpinBox), pBox(new QSpinBox), 
        testCheckBox(new QCheckBox("Run &tests only")), 
        goButton(new QPushButton("&go")), progressBar(new QProgressBar) {
    // put the boxes into inputBoxGrid
    QHBoxLayout* inputBoxGrid = new QHBoxLayout;
    inputBoxGrid->addWidget(new QLabel("n"));
    inputBoxGrid->addWidget(nBox);
    inputBoxGrid->addWidget(new QLabel("l"));
    inputBoxGrid->addWidget(lBox);
    inputBoxGrid->addWidget(new QLabel("p"));
    inputBoxGrid->addWidget(pBox);
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

    connect(goButton, &QAbstractButton::clicked, this, &CalcWidget::Calculate);

    layout->addWidget(inputBoxes);
    layout->addWidget(testCheckBox);
    layout->addWidget(goButton);
    layout->addWidget(progressBar);
    setLayout(layout);
}

void CalcWidget::Calculate() {
    // FIXME: mutex lock or something to make sure this never gets run twice
    goButton->setEnabled(false);

    Arguments args;
    args.numP = nBox->value();
    args.degree = lBox->value();
    args.partitions = pBox->value();
    args.outputStream = outStream;

    if (testCheckBox->isChecked()) args.options |= OPT_TEST;

    ::Calculate(args);
    std::cout << "***** Calculation Complete *****" << std::endl;

    goButton->setEnabled(true);
    // FIXME: release mutex (or whatever else we end up using)
}

void CalcWidget::ChangeOutput(std::ostream* newOutStream) {
    std::cout << "outStream updated" << std::endl;
    outStream = newOutStream;
}

} // namespace GUI
