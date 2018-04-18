#include "console_widget.hpp"

namespace GUI {

ConsoleWidget::ConsoleWidget(): listener(new ConsoleListener(this)), 
    outStream(listener) {
    connect(listener, &ConsoleListener::SendText,
            this, &ConsoleWidget::WriteText);
    setReadOnly(true);
    // setBackgroundRole(QPalette::Shadow);
    // setTextBackgroundColor(Qt::black);
    // setTextColor(Qt::darkMagenta);
    // setText("Ready to calculate.");

    // for (int i = 0; i < 100; ++i) {
    outStream << "Ready to calculate." << endl;
    // }
}

void ConsoleWidget::WriteText(QString text) {
    append(text);
}

} // namespace GUI
