#include "file_widget.hpp"

namespace GUI {

FileWidget::FileWidget(std::ostream* outStream): outStream(outStream), 
    outPath(new QLineEdit), dontSave(new QCheckBox),
    suppressOverwriteWarning(new QCheckBox), appendContents(new QCheckBox) {
    QVBoxLayout* layout = new QVBoxLayout;
    setLayout(layout);

    outPath->setText(tr("output file"));
    dontSave->setText(tr("don't save output"));
    suppressOverwriteWarning->setText(tr("don't warn about overwriting files"));
    appendContents->setText(tr("append to file contents (don't overwrite)"));

    layout->addWidget(outPath);
    layout->addWidget(dontSave);
    layout->addWidget(suppressOverwriteWarning);
    layout->addWidget(appendContents);

    // FIXME: need validator for outPath
    QRegExp pathRegex();
    QRegExpValidator* pathValidator = new QRegExpValidator(pathRegex);
    outPath->setValidator(pathValidator);

    connect(outPath, &QLineEdit::editingFinished, 
            this, &FileWidget::ChangeOutputStream);
    connect(dontSave, &QCheckBox::stateChanged, 
            this, &FileWidget::ChangeOutputStream);
    connect(appendContents, &QCheckBox::stateChanged, 
            this, &FileWidget::ChangeOutputStream);

    connect(dontSave, &QCheckBox::stateChanged,
            [=](bool checked){outPath->setEnabled(!checked);});
}

void FileWidget::ChangeOutputStream() {
    if (outStream != nullptr) {
        outStream->flush(); // FIXME?? is this a bad idea?
        if (outStream->rdbuf() != std::cout.rdbuf()
            && outStream->rdbuf() != std::cerr.rdbuf()) delete outStream;
    }

    if (dontSave->isChecked()) {
        // std::cout << "outStream set to std::cout" << std::endl;
        outStream = &std::cout;
        emit OutputChanged(outStream);
        return;
    }

    auto writeMode = std::ios_base::out;
    if (appendContents->isChecked()) {
        writeMode |= std::ios_base::app;
    } else {
        writeMode |= std::ios_base::trunc;
    }
    // std::cout << "outStream set to \"" << 
        // (std::string)outPath->text().toLocal8Bit() << "\"" << std::endl;
    outStream = new std::ofstream(outPath->text().toLocal8Bit(), writeMode);

    emit OutputChanged(outStream);
}

} // namespace GUI
