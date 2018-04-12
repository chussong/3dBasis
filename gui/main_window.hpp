#ifndef MAIN_WINDOW_HPP
#define MAIN_WINDOW_HPP

#include "QtWidgets"

// these might be in QtWidgets already?
#include <QApplication>
#include <QCommandLineParser>
#include <QCommandLineOption>

#include "constants.hpp"

namespace GUI {

class MainWindow : public QMainWindow {
    Q_OBJECT // this is a macro :(

    public:
        MainWindow();

        void LoadFile(const QString& filename);

    protected:
        void closeEvent(QCloseEvent* event) override;

    private slots:
        void NewFiles();
        void Open();
        void Save();
        void SaveAs();
        void About();
        void DocumentWasModified();
#ifndef QT_NO_SESSIONMANAGER
        void CommitData(QSessionManager&);
#endif

    private:
        void CreateActions();
        void CreateStatusBar();
        void ReadSettings();
        void WriteSettings();
        bool MaybeSave();
        bool SaveFile(const QString& fileName);
        void SetCurrentFile(const QString& fileName);
        QString StrippedName(const QString& fullFileName);

        QPlainTextEdit* textEdit;
        QString curFile;
};

bool StartGUI(const Arguments& args);

} // namespace GUI

#endif
