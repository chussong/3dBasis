#include "main_window.hpp"

namespace GUI {

MainWindow::MainWindow(Arguments args) {
    // CreateActions();
    // CreateStatusBar();

    QGridLayout* layout = new QGridLayout;

    ConsoleWidget* consoleWidget = new ConsoleWidget();
    FileWidget* fileWidget = new FileWidget(consoleWidget->OutStream());
    args.outStream = consoleWidget->OutStream();
    args.console = consoleWidget->OutStream();
    CalcWidget* calcWidget = new CalcWidget(args);
    layout->addWidget(fileWidget, 0, 0, 10, 6);
    layout->addWidget(calcWidget, 0, 6, 10, 6);
    layout->addWidget(consoleWidget, 10, 0, 20, 12);
    setLayout(layout);

    connect(fileWidget, &FileWidget::OutputChanged,
            calcWidget, &CalcWidget::ChangeOutput);
    connect(fileWidget, &FileWidget::OverwriteWarningSignal,
            calcWidget, &CalcWidget::GiveOverwriteWarnings);
    connect(calcWidget, &CalcWidget::StartingCalculation,
            consoleWidget, &ConsoleWidget::clear);
    connect(calcWidget, &CalcWidget::OverwriteFile,
            fileWidget, &FileWidget::ReopenFileStream);

    ReadSettings();
}

void MainWindow::closeEvent(QCloseEvent* event) {
    WriteSettings();
    event->accept();
}

void MainWindow::About() {
    QMessageBox::about(this, tr("About 3dBasis"),
            tr("<b>3dBasis</b> is a program for numerically computing a "
                "truncated polynomial basis of a conformal field theory, as "
                "well as various properties of the basis. Find updates "
                "<a href=\"https://www.github.com/chussong/3dBasis\">here</a>.")
            );
}

/*void MainWindow::CreateActions() {
    QMenu* fileMenu = menuBar()->addMenu(tr("&File"));

    // fileMenu->addSeparator(); FIXME: add other "file" actions above this

    const QIcon exitIcon = QIcon::fromTheme("document-exit");
    QAction* exitAct = fileMenu->addAction(exitIcon, tr("E&xit..."), this,
            &QWidget::close);
    exitAct->setShortcuts(QKeySequence::Quit);
    exitAct->setStatusTip(tr("Exit the application"));

    QMenu* helpMenu = menuBar()->addMenu(tr("&Help"));
    QAction* aboutAct = helpMenu->addAction(tr("&About"), this, 
            &MainWindow::About);
    aboutAct->setStatusTip(tr("Show 3dBasis's About box"));

    QAction* aboutQtAct = helpMenu->addAction(tr("About &Qt"), qApp,
            &QApplication::aboutQt);
    aboutQtAct->setStatusTip(tr("Show the Qt library's About box"));
}*/

// void MainWindow::CreateStatusBar() {
    // statusBar()->showMessage(tr("Ready"));
// }

void MainWindow::ReadSettings() {
    QSettings settings(QCoreApplication::organizationName(),
            QCoreApplication::applicationName());
    const QByteArray geometry = 
        settings.value("geometry", QByteArray()).toByteArray();
    if (geometry.isEmpty()) {
        const QRect availableGeometry = 
            QApplication::desktop()->availableGeometry(this);
        resize(availableGeometry.width()/3, availableGeometry.height()/2);
        move((availableGeometry.width() - width())/2,
                (availableGeometry.height() - height())/2);
    } else {
        restoreGeometry(geometry);
    }
}

void MainWindow::WriteSettings() {
    QSettings settings(QCoreApplication::organizationName(),
            QCoreApplication::applicationName());
    settings.setValue("geometry", saveGeometry());
}

bool StartGUI(int argc, char** argv, const Arguments& args) {
    // Q_INIT_RESOURCE(application);
    // QResource::registerResource("/home/charles/Dropbox/Research/3dBasis/gui/resources.rcc");
    QApplication app(argc, argv);
    QCoreApplication::setOrganizationName("CharlesHussong");
    QCoreApplication::setApplicationName("3dBasis");
    QCoreApplication::setApplicationVersion("{VERSION PLACEHOLDER}");

    app.setStyleSheet("QTextEdit { background-color: black; color: white; }");

    // QCommandLineParser parser;
    // parser.setApplicationDescription(QCoreApplication::applicationName());
    // parser.addHelpOption();
    // parser.addVersionOption();
    // parser.addPositionalArgument("file", "The file to open.");
    // parser.process(app);

    MainWindow mainWin(args);
    // if (!parser.positionalArguments().isEmpty()) {
        // mainWin.LoadFile(parser.positionalArguments().first());
    // }

    mainWin.show();
    return app.exec();
}

} // namespace GUI
