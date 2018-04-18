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
#include "console_widget.hpp"

namespace GUI {

class MainWindow : public QWidget {
    Q_OBJECT

    public:
        MainWindow(Arguments args);

    protected:
        void closeEvent(QCloseEvent* event) override;

    private slots:
        void About();

    private:
        void CreateActions();
        void CreateStatusBar();
        void ReadSettings();
        void WriteSettings();
};

bool StartGUI(int argc, char** argv, const Arguments& args);

} // namespace GUI

#endif
