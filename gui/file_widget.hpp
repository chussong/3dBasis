#ifndef FILEWIDGET_HPP
#define FILEWIDGET_HPP

#include <iostream>
#include <fstream>
#include <memory>

#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QFileDialog>
#include <QtCore/QFileInfo>
#include <QtCore/QTextStream>
#include <QtGui/QRegExpValidator>

namespace GUI {

class FileWidget : public QWidget {
    Q_OBJECT

    public:
        FileWidget(QTextStream* consoleStream);
        QTextStream* OutStream();

    signals:
        void OutputChanged(QTextStream* newOutStream);
        void OverwriteWarningSignal(const bool newStatus);

    public slots:
        void ReopenFileStream();

    private slots:
        void ChooseOutputFile();
        // void ChangeOutputFileName();
        void ChangeOutputStream();
        void OverwriteWarningSlot();

    private:
        void DisableOutput();
        void OpenFileStream(const QString& fileName, 
                            const QFile::OpenMode writeMode);
        void CloseFileStream();
        QTextStream* consoleStream;
        QTextStream* fileStream;
        QFile* openFile;

        // QVBoxLayout* layout;
        QLineEdit* outPath;
        QPushButton* outPathButton;
        QCheckBox* dontSave; // checking this disables the text field
        QCheckBox* suppressOverwriteWarning;
        QCheckBox* appendContents; // i.e. don't overwrite the file
};

} // namespace GUI

#endif
