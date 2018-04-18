#include "file_widget.hpp"

namespace GUI {

namespace {
char outPathDefaultText[] = "output file (*.txt)";
} // anonymous namespace

FileWidget::FileWidget(QTextStream* consoleStream): consoleStream(consoleStream),
    fileStream(nullptr), openFile(nullptr), outPath(new QLineEdit), 
    outPathButton(new QPushButton), dontSave(new QCheckBox), 
    suppressOverwriteWarning(new QCheckBox), appendContents(new QCheckBox) {
    QVBoxLayout* layout = new QVBoxLayout;
    setLayout(layout);

    outPath->setText(tr(outPathDefaultText));
    outPath->setEnabled(false);
    dontSave->setText(tr("don't save output"));
    dontSave->setChecked(true);
    suppressOverwriteWarning->setText(tr("don't warn about overwriting files"));
    appendContents->setText(tr("append to file contents (don't overwrite)"));

    // FIXME: WIP of moving to dialog-only file selection
    outPath->setReadOnly(true);
    outPathButton->setText("choose output file");
    outPathButton->setEnabled(false);
    connect(outPathButton, &QAbstractButton::clicked,
            this, &FileWidget::ChooseOutputFile);

    layout->addWidget(outPath);
    layout->addWidget(outPathButton);
    layout->addWidget(dontSave);
    layout->addWidget(suppressOverwriteWarning);
    layout->addWidget(appendContents);

    QRegExp pathRegex("\\w+\\.txt");
    QRegExpValidator* pathValidator = new QRegExpValidator(pathRegex);
    outPath->setValidator(pathValidator);

    // connect(outPath, &QLineEdit::textChanged, 
            // this, &FileWidget::ChangeOutputFileName);

    connect(dontSave, &QCheckBox::stateChanged, 
            this, &FileWidget::ChangeOutputStream);
    connect(appendContents, &QCheckBox::stateChanged, 
            this, &FileWidget::ChangeOutputStream);

    connect(dontSave, &QCheckBox::stateChanged,
            [=](bool checked){outPath->setEnabled(!checked);});
    connect(dontSave, &QCheckBox::stateChanged,
            [=](bool checked){outPathButton->setEnabled(!checked);});

    connect(suppressOverwriteWarning, &QCheckBox::stateChanged,
            this, &FileWidget::OverwriteWarningSlot);
}

QTextStream* FileWidget::OutStream() {
    return openFile != nullptr ? fileStream : consoleStream;
}

void FileWidget::ChooseOutputFile() {
    // QFileDialog dialog(this);
    // dialog.setAcceptMode(QFileDialog::AcceptSave);
    // dialog.setNameFilter(tr("Text files (*.txt)"));
    // dialog.setDefaultSuffix(".txt");
    // dialog.setViewMode(QFileDialog::List);
    // if (!dialog.exec() || dialog.selectedFile() == outPath->text()) return;
// 
    // outPath->setText(dialog.selectedFile());
    // ChangeOutputFileName();
    QString fileName = QFileDialog::getSaveFileName(this, 
            tr("choose output file"), "" /* FIXME: current directory */,
            tr("Text files (*.txt)"));
    if (fileName.isEmpty() || fileName == outPath->text()) return;

    outPath->setText(fileName);
    ChangeOutputStream();
}

// this is for when the FILE needs to be changed; it's not called if only the
// stream is different
/*void FileWidget::ChangeOutputFileName() {
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
}*/

// this is called if the stream needs to change for any reason, such as a new
// filename or different write mode
//
// if this function is called, we assume that the contents of *outPath form a 
// valid file name
//
// some variation of OutputChanged MUST be emitted before this exits
void FileWidget::ChangeOutputStream() {
    if (!dontSave->isChecked() && outPath->text() == outPathDefaultText) {
        DisableOutput();
        return;
    }

    if (dontSave->isChecked()) {
        if (openFile != nullptr) {
            *consoleStream << "outStream set to stdout" << endl;
            emit OutputChanged(consoleStream);

            CloseFileStream();
        } else {
            emit OutputChanged(consoleStream);
        }
        return;
    }

    QFile::OpenMode writeMode = QFile::WriteOnly;
    if (appendContents->isChecked()) {
        writeMode = writeMode | QFile::Append;
    } else {
        writeMode |= QFile::Truncate;
    }

    OpenFileStream(outPath->text().toLocal8Bit(), writeMode);
    emit OutputChanged(fileStream);
}

void FileWidget::OverwriteWarningSlot() {
    emit OverwriteWarningSignal(!suppressOverwriteWarning->isChecked());
}

void FileWidget::DisableOutput() {
    emit OutputChanged(nullptr);
}

void FileWidget::ReopenFileStream() {
    // looks like this is actually sufficient
    ChangeOutputStream();
}

void FileWidget::OpenFileStream(const QString& fileName, 
                                const QFile::OpenMode writeMode) {
    *consoleStream << "outStream set to \"" << fileName << "\"" << endl;
    if (openFile != nullptr) CloseFileStream();
    openFile = new QFile(fileName);
    if (!openFile->open(writeMode)) {
        throw std::runtime_error("FileWidget failed to open a file.");
    }
    fileStream = new QTextStream(openFile);
}

void FileWidget::CloseFileStream() {
    fileStream->flush();
    delete fileStream;
    fileStream = nullptr;
    openFile->close();
    delete openFile;
    openFile = nullptr;
}

} // namespace GUI
