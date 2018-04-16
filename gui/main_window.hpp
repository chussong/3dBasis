#ifndef MAIN_WINDOW_HPP
#define MAIN_WINDOW_HPP

#include <QtCore/QSettings>
#include <QtCore/QResource>
#include <QtCore/QCommandLineParser>
#include <QtCore/QCommandLineOption>
#include <QtGui/QSessionManager>
#include <QtWidgets/QDesktopWidget>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QPlainTextEdit>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QAction>
#include <QtWidgets/QDockWidget>
#include <QtWidgets/QApplication>
#include <QtWidgets/QCheckBox>

#include "constants.hpp"
#include "calc_widget.hpp"
#include "file_widget.hpp"

namespace GUI {

class MainWindow : public QMainWindow {
    Q_OBJECT

    public:
        MainWindow(const Arguments& args);

        void LoadFile(const QString& filename);

    protected:
        void closeEvent(QCloseEvent* event) override;

    private slots:
        void NewFile();
        void Open();
        bool Save();
        bool SaveAs();
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
        bool overwriteWithoutPrompting = false;
};

bool StartGUI(int argc, char** argv, const Arguments& args);

} // namespace GUI

#endif
