#ifndef CONSOLE_HPP
#define CONSOLE_HPP

#include <iostream>

#include <QtWidgets/QTextEdit>
#include <QtWidgets/QTextStream>

namespace GUI {

class Console : public QTextEdit {
    Q_OBJECT

    public:
        Console();

    private:
        QTextEdit* display;
        QTextStream* outStream;
};

} // namespace GUI

#endif
