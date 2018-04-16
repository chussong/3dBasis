#include "main_window.hpp"

namespace GUI {

MainWindow::MainWindow(const Arguments& args): textEdit(new QPlainTextEdit) {
    setCentralWidget(textEdit);

    CreateActions();
    CreateStatusBar();

    CalcWidget* calcWidget = new CalcWidget(args);
    QDockWidget* calcDock = new QDockWidget(tr("Calc Widget"), this);
    calcDock->setAllowedAreas(Qt::LeftDockWidgetArea
                                | Qt::RightDockWidgetArea);
    calcDock->setWidget(calcWidget);
    addDockWidget(Qt::RightDockWidgetArea, calcDock);

    FileWidget* fileWidget = new FileWidget();
    QDockWidget* fileDock = new QDockWidget(tr("File Widget"), this);
    fileDock->setAllowedAreas(Qt::LeftDockWidgetArea
                                | Qt::RightDockWidgetArea);
    fileDock->setWidget(fileWidget);
    addDockWidget(Qt::LeftDockWidgetArea, fileDock);

    connect(fileWidget, &FileWidget::OutputChanged,
            calcWidget, &CalcWidget::ChangeOutput);
    connect(fileWidget, &FileWidget::OverwriteWarningSignal,
            calcWidget, &CalcWidget::GiveOverwriteWarnings);

    ReadSettings();

    textEdit->setReadOnly(true);
    connect(textEdit->document(), &QTextDocument::contentsChanged, this,
            &MainWindow::DocumentWasModified);

#ifndef QT_NO_SESSIONMANAGER
    QGuiApplication::setFallbackSessionManagementEnabled(false);
    connect(qApp, &QGuiApplication::commitDataRequest, this,
            &MainWindow::CommitData);
#endif

    SetCurrentFile(QString());
    setUnifiedTitleAndToolBarOnMac(true);
}

void MainWindow::closeEvent(QCloseEvent* event) {
    if (MaybeSave()) {
        WriteSettings();
        event->accept();
    } else {
        event->ignore();
    }
}

void MainWindow::NewFile() {
    if (MaybeSave()) {
        textEdit->clear();
        SetCurrentFile(QString());
    }
}

void MainWindow::Open() {
    if (MaybeSave()) {
        QString fileName = QFileDialog::getOpenFileName(this);
        if (!fileName.isEmpty()) LoadFile(fileName);
    }
}

bool MainWindow::Save() {
    if (curFile.isEmpty()) {
        return SaveAs();
    } else {
        return SaveFile(curFile);
    }
}

bool MainWindow::SaveAs() {
    QFileDialog dialog(this);
    dialog.setWindowModality(Qt::WindowModal);
    dialog.setAcceptMode(QFileDialog::AcceptSave);
    if (dialog.exec() != QDialog::Accepted) return false;
    return SaveFile(dialog.selectedFiles().first());
}

void MainWindow::About() {
    QMessageBox::about(this, tr("About 3dBasis"),
            tr("<b>3dBasis</b> is a program for numerically computing a "
                "truncated polynomial basis of a conformal field theory, as "
                "well as various properties of the basis. Find updates "
                "<a href=\"https://www.github.com/chussong/3dBasis\">here</a>.")
            );
}

void MainWindow::DocumentWasModified() {
    setWindowModified(textEdit->document()->isModified());
}

void MainWindow::CreateActions() {
    QMenu* fileMenu = menuBar()->addMenu(tr("&File"));
    QToolBar* fileToolBar = addToolBar(tr("File"));

    const QIcon newIcon = QIcon::fromTheme("document-new", 
            QIcon(":/images/new.png"));
    QAction* newAct = new QAction(newIcon, tr("&New"), this);
    newAct->setShortcuts(QKeySequence::New);
    newAct->setStatusTip(tr("Create a new file"));
    connect(newAct, &QAction::triggered, this, &MainWindow::NewFile);
    fileMenu->addAction(newAct);
    fileToolBar->addAction(newAct);

    const QIcon openIcon = QIcon::fromTheme("document-open",
            QIcon(":/images/open.png"));
    QAction* openAct = new QAction(openIcon, tr("&Open..."), this);
    openAct->setShortcuts(QKeySequence::Open);
    openAct->setStatusTip(tr("Open an existing file"));
    connect(openAct, &QAction::triggered, this, &MainWindow::Open);
    fileMenu->addAction(openAct);
    fileToolBar->addAction(openAct);

    const QIcon saveIcon = QIcon::fromTheme("document-save",
            QIcon(":/images/save.png"));
    QAction* saveAct = new QAction(saveIcon, tr("&Save..."), this);
    saveAct->setShortcuts(QKeySequence::Save);
    saveAct->setStatusTip(tr("Save current file to disk"));
    connect(saveAct, &QAction::triggered, this, &MainWindow::Save);
    fileMenu->addAction(saveAct);
    fileToolBar->addAction(saveAct);

    const QIcon saveAsIcon = QIcon::fromTheme("document-save-as");
    QAction* saveAsAct = fileMenu->addAction(saveAsIcon, tr("Save &As..."), 
            this, &MainWindow::SaveAs);
    saveAsAct->setShortcuts(QKeySequence::SaveAs);
    saveAsAct->setStatusTip(tr("Save the document under a new name"));
    
    fileMenu->addSeparator();

    const QIcon exitIcon = QIcon::fromTheme("document-exit");
    QAction* exitAct = fileMenu->addAction(exitIcon, tr("E&xit..."), this,
            &QWidget::close);
    exitAct->setShortcuts(QKeySequence::Quit);
    exitAct->setStatusTip(tr("Exit the application"));

    QMenu* editMenu = menuBar()->addMenu(tr("&Edit"));
    QToolBar* editToolBar = addToolBar(tr("Edit"));

#ifndef QT_NO_CLIPBOARD
    const QIcon cutIcon = QIcon::fromTheme("edit-cut", 
            QIcon(":/images/cut.png"));
    QAction* cutAct = new QAction(cutIcon, tr("Cu&t"), this);
    cutAct->setShortcuts(QKeySequence::Cut);
    cutAct->setStatusTip(tr("Cut the current selection's contents to the "
                "system clipboard."));
    connect(cutAct, &QAction::triggered, textEdit, &QPlainTextEdit::cut);
    editMenu->addAction(cutAct);
    editToolBar->addAction(cutAct);

    const QIcon copyIcon = QIcon::fromTheme("edit-copy", 
            QIcon(":/images/copy.png"));
    QAction* copyAct = new QAction(copyIcon, tr("&Copy"), this);
    copyAct->setShortcuts(QKeySequence::Copy);
    copyAct->setStatusTip(tr("Copy the current selection's contents to the "
                "system clipboard."));
    connect(copyAct, &QAction::triggered, textEdit, &QPlainTextEdit::copy);
    editMenu->addAction(copyAct);
    editToolBar->addAction(copyAct);

    const QIcon pasteIcon = QIcon::fromTheme("edit-paste", 
            QIcon(":/images/paste.png"));
    QAction* pasteAct = new QAction(pasteIcon, tr("&Paste"), this);
    pasteAct->setShortcuts(QKeySequence::Paste);
    pasteAct->setStatusTip(tr("Paste the system clipboards's contents into the "
                "curent selection."));
    connect(pasteAct, &QAction::triggered, textEdit, &QPlainTextEdit::paste);
    editMenu->addAction(pasteAct);
    editToolBar->addAction(pasteAct);

    menuBar()->addSeparator();
#endif // !QT_NO_CLIPBOARD

    QMenu* helpMenu = menuBar()->addMenu(tr("&Help"));
    QAction* aboutAct = helpMenu->addAction(tr("&About"), this, 
            &MainWindow::About);
    aboutAct->setStatusTip(tr("Show 3dBasis's About box"));

    QAction* aboutQtAct = helpMenu->addAction(tr("About &Qt"), qApp,
            &QApplication::aboutQt);
    aboutQtAct->setStatusTip(tr("Show the Qt library's About box"));

#ifndef QT_NO_CLIPBOARD
    cutAct->setEnabled(false);
    copyAct->setEnabled(false);
    connect(textEdit, &QPlainTextEdit::copyAvailable, cutAct, 
            &QAction::setEnabled);
    connect(textEdit, &QPlainTextEdit::copyAvailable, copyAct, 
            &QAction::setEnabled);
#endif // !QT_NO_CLIPBOARD
}

void MainWindow::CreateStatusBar() {
    statusBar()->showMessage(tr("Ready"));
}

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

bool MainWindow::MaybeSave() {
    if (!textEdit->document()->isModified()) return true;
    const QMessageBox::StandardButton ret
        = QMessageBox::warning(this, tr("Application"),
                              tr("The document has been modified.\n"
                                  "Do you want to save your changes?"),
                              QMessageBox::Save | QMessageBox::Discard 
                              | QMessageBox::Cancel);
    switch (ret) {
        case QMessageBox::Save:
            return Save();
        case QMessageBox::Cancel:
            return false;
        default:
            break;
    }
    return true;
}

void MainWindow::LoadFile(const QString& fileName) {
    QFile file(fileName);
    if (!file.open(QFile::ReadOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot read file %1:\n%2.")
                             .arg(QDir::toNativeSeparators(fileName),
                                  file.errorString()));
        return;
    }

    QTextStream in(&file);
#ifndef QT_NO_CURSOR
    QApplication::setOverrideCursor(Qt::WaitCursor);
#endif
    textEdit->setPlainText(in.readAll());
#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
#endif

    SetCurrentFile(fileName);
    statusBar()->showMessage(tr("File loaded"), 2000);
}

bool MainWindow::SaveFile(const QString& fileName) {
    QFile file(fileName);
    if (!file.open(QFile::WriteOnly | QFile::Text)) {
        QMessageBox::warning(this, tr("Application"),
                             tr("Cannot write file %1:\n%2.")
                             .arg(QDir::toNativeSeparators(fileName),
                                  file.errorString()));
        return false;
    }

    QTextStream out(&file);
#ifndef QT_NO_CURSOR
    QApplication::setOverrideCursor(Qt::WaitCursor);
#endif
    out << textEdit->toPlainText();
#ifndef QT_NO_CURSOR
    QApplication::restoreOverrideCursor();
#endif

    SetCurrentFile(fileName);
    statusBar()->showMessage(tr("File saved"), 2000);
    return true;
}

void MainWindow::SetCurrentFile(const QString& fileName) {
    curFile = fileName;
    textEdit->document()->setModified(false);
    setWindowModified(false);

    QString shownName = curFile;
    if (curFile.isEmpty()) shownName = "untitled.txt";
    setWindowFilePath(shownName);
}

QString MainWindow::StrippedName(const QString& fullFileName) {
    return QFileInfo(fullFileName).fileName();
}

#ifndef QT_NO_SESSIONMANAGER
void MainWindow::CommitData(QSessionManager& manager) {
    if (manager.allowsInteraction()) {
        if (!MaybeSave()) {
            manager.cancel();
        } else {
            // save without asking
            if (textEdit->document()->isModified()) Save();
        }
    }
}
#endif

bool StartGUI(int argc, char** argv, const Arguments& args) {
    // Q_INIT_RESOURCE(application);
    QResource::registerResource("/home/charles/Dropbox/Research/3dBasis/gui/resources.rcc");
    // QApplication app(argc, argv);
    QApplication app(argc, argv);
    QCoreApplication::setOrganizationName("CharlesHussong");
    QCoreApplication::setApplicationName("3dBasis");
    QCoreApplication::setApplicationVersion("{VERSION PLACEHOLDER}");

    QCommandLineParser parser;
    parser.setApplicationDescription(QCoreApplication::applicationName());
    parser.addHelpOption();
    parser.addVersionOption();
    parser.addPositionalArgument("file", "The file to open.");
    parser.process(app);

    MainWindow mainWin(args);
    if (!parser.positionalArguments().isEmpty()) {
        mainWin.LoadFile(parser.positionalArguments().first());
    }

    mainWin.show();
    return app.exec();
}

} // namespace GUI
