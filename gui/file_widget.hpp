#ifndef FILEWIDGET_HPP
#define FILEWIDGET_HPP

#include <iostream>
#include <fstream>

#include <QtWidgets/QVBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QCheckBox>

namespace GUI {

class FileWidget : public QWidget {
    Q_OBJECT

    public:
        FileWidget(std::ostream* outStream);

    signals:
        void OutputChanged(std::ostream* newOutStream);

    private slots:
        void ChangeOutputStream();

    private:
        std::ostream* outStream;

        // QVBoxLayout* layout;
        QLineEdit* outPath;
        QCheckBox* dontSave; // checking this disables the text field
        QCheckBox* suppressOverwriteWarning;
        QCheckBox* appendContents; // i.e. don't overwrite the file
};

} // namespace GUI

#endif
