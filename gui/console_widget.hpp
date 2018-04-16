#ifndef CONSOLE_WIDGET_HPP
#define CONSOLE_WIDGET_HPP

#include <iostream>

#include <QtWidgets/QTextEdit>
#include <QtWidgets/QTextStream>

namespace GUI {

class ConsoleWidget : public QTextEdit {
    Q_OBJECT

    public:
        ConsoleWidget();

    private:
        QTextStream* outStream;
};

} // namespace GUI

#endif
