#include "file_widget.hpp"

namespace GUI {

namespace {
char outPathDefaultText[] = "output file (*.txt)";
} // anonymous namespace

FileWidget::FileWidget(): outStream(), 
    outPath(new QLineEdit), dontSave(new QCheckBox), 
    suppressOverwriteWarning(new QCheckBox), appendContents(new QCheckBox) {
    QVBoxLayout* layout = new QVBoxLayout;
    setLayout(layout);

    outPath->setText(tr(outPathDefaultText));
    outPath->setEnabled(false);
    dontSave->setText(tr("don't save output"));
    dontSave->setChecked(true);
    suppressOverwriteWarning->setText(tr("don't warn about overwriting files"));
    appendContents->setText(tr("append to file contents (don't overwrite)"));

    layout->addWidget(outPath);
    layout->addWidget(dontSave);
    layout->addWidget(suppressOverwriteWarning);
    layout->addWidget(appendContents);

    QRegExp pathRegex("\\w+\\.txt");
    QRegExpValidator* pathValidator = new QRegExpValidator(pathRegex);
    outPath->setValidator(pathValidator);

    connect(outPath, &QLineEdit::textChanged, 
            this, &FileWidget::ChangeOutputFileName);

    connect(dontSave, &QCheckBox::stateChanged, 
            this, &FileWidget::ChangeOutputStream);
    connect(appendContents, &QCheckBox::stateChanged, 
            this, &FileWidget::ChangeOutputStream);

    connect(dontSave, &QCheckBox::stateChanged,
            [=](bool checked){outPath->setEnabled(!checked);});

    connect(suppressOverwriteWarning, &QCheckBox::stateChanged,
            this, &FileWidget::OverwriteWarningSlot);
}

// this is for when the FILE needs to be changed; it's not called if only the
// stream is different
void FileWidget::ChangeOutputFileName() {
    if (!outPath->hasAcceptableInput()) {
        DisableOutput();
        return;
    }

    QFileInfo fileInfo(outPath->text().toLocal8Bit());
    if (fileInfo.exists()) {
        if (fileInfo.isFile() && fileInfo.isWritable()) {
            QMessageBox confirm;
            confirm.setText("This file already exists.");
            confirm.setInformativeText("Do you want to overwrite it?");
            confirm.setStandardButtons(QMessageBox::Yes | QMessageBox::Cancel);
            confirm.setDefaultButton(QMessageBox::Cancel);
            if (confirm.exec() == QMessageBox::Yes) {
                ChangeOutputStream();
            } else {
                outPath->setText(tr(outPathDefaultText));
                return;
            }
        } else {
            QMessageBox badName;
            badName.setText("This name already exists and is not a writable "
                    "file.");
            badName.exec();
            outPath->setText(tr(outPathDefaultText));
            return;
        }
    } else {
        // file doesn't exist, go ahead and make it
        ChangeOutputStream();
    }
}

// is is called if the stream needs to change for any reason, such as a new
// filename or different write mode
//
// if this function is called, we assume that the contents of *outPath form a 
// valid file name
void FileWidget::ChangeOutputStream() {
    if (!dontSave->isChecked() && !outPath->hasAcceptableInput()) {
        DisableOutput();
        return;
    }

    if (dontSave->isChecked()) {
        // if (outStream != nullptr) {
        if (outStream.rdbuf() != std::cout.rdbuf()) {
            std::cout << "outStream set to std::cout" << std::endl;
            outStream.close();
            // fileBuffer = outStream.rdbuf(std::cout.rdbuf());
            // emit OutputChanged(outStream.get());
            // outStream->close();
            // outStream = nullptr;
            emit OutputChanged(&std::cout);
        }
        return;
    }

    auto writeMode = std::ios_base::out;
    if (appendContents->isChecked()) {
        writeMode |= std::ios_base::app;
    } else {
        writeMode |= std::ios_base::trunc;
    }
    std::cout << "outStream set to \"" << 
        (std::string)outPath->text().toLocal8Bit() << "\"" << std::endl;
    if (outStream.rdbuf() != std::cout.rdbuf()) outStream.close();
    // outStream.rdbuf(fileBuffer);
    outStream.open(outPath->text().toLocal8Bit(), writeMode);
    // std::unique_ptr<std::ofstream> newOutStream = std::make_unique<std::ofstream>(
                                                // outPath->text().toLocal8Bit(), 
                                                // writeMode);
    emit OutputChanged(&outStream);

    // if (outStream != nullptr) outStream->close();
    // outStream = std::move(newOutStream);
}

void FileWidget::OverwriteWarningSlot() {
    emit OverwriteWarningSignal(!suppressOverwriteWarning->isChecked());
}

void FileWidget::DisableOutput() {
    emit OutputChanged(nullptr);
}

} // namespace GUI
