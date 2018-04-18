#ifndef CONSOLE_WIDGET_HPP
#define CONSOLE_WIDGET_HPP

#include <iostream>

#include <QtWidgets/QTextEdit>
#include <QtCore/QTextStream>

namespace GUI {

class ConsoleListener : public QIODevice {
    Q_OBJECT

    public:
        ConsoleListener(QObject* const parent):
            QIODevice(parent) {
                open(QIODevice::WriteOnly | QIODevice::Text);
            }

    signals:
        void SendText(QString text);

    protected:
        qint64 readData(char*, qint64) { return 0; }
        qint64 writeData(const char* data, qint64 maxSize) {
            // if (console != nullptr) {
                // console->append(data);
            // }
            // return maxSize;
            emit SendText(QString(data));
            return maxSize;
        }

    private:
        // ConsoleWidget* console;
};

class ConsoleWidget : public QTextEdit {
    Q_OBJECT

    public:
        ConsoleWidget();
        QTextStream* OutStream() { return &outStream; }

    public slots:
        void WriteText(QString text);

    private:
        ConsoleListener* listener;
        QTextStream outStream;
};

} // namespace GUI

#endif
