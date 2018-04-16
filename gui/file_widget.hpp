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
#include <QtCore/QFileInfo>
#include <QtGui/QRegExpValidator>

namespace GUI {

class FileWidget : public QWidget {
    Q_OBJECT

    public:
        FileWidget();
        std::ostream* OutStream() { return &outStream; }

    signals:
        void OutputChanged(std::ostream* newOutStream);
        void OverwriteWarningSignal(const bool newStatus);

    private slots:
        void ChangeOutputFileName();
        void ChangeOutputStream();
        void OverwriteWarningSlot();

    private:
        void DisableOutput();
        std::ofstream outStream;

        // QVBoxLayout* layout;
        QLineEdit* outPath;
        QCheckBox* dontSave; // checking this disables the text field
        QCheckBox* suppressOverwriteWarning;
        QCheckBox* appendContents; // i.e. don't overwrite the file
};

} // namespace GUI

#endif
