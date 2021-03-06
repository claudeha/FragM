#pragma once

#ifndef TEXTEDIT_H
#define TEXTEDIT_H

#include "MainWindow.h"
#include <QAbstractScrollArea>
#include <QAction>
#include <QCheckBox>
#include <QComboBox>
#include <QLabel>
#include <QMainWindow>
#include <QMenu>
#include <QPlainTextEdit>
#include <QScrollBar>
#include <QSpinBox>
#include <QStackedWidget>
#include <QTabBar>
#include <QTextBlock>
#include <QTextEdit>

class MainWindow;
class FragmentHighlighter;

namespace Fragmentarium
{
namespace GUI
{

// A modified QTextEdit with an extended context menu, syntax highlighting and
// line numbers
class TextEdit : public QPlainTextEdit
{
    Q_OBJECT

public:
    TextEdit(MainWindow *parent = 0);

    void lineNumberAreaPaintEvent(QPaintEvent *event);
    int lineNumberAreaWidth();
    void contextMenuEvent(QContextMenuEvent *event);
    void insertFromMimeData (const QMimeData * source );
    void saveSettings ( QString ss )
    {
        savedSettings = ss;
    }
    QString lastSettings()
    {
        return savedSettings;
    }
    void parmsToTest ( QString ss )
    {
        testSettings = ss;
    }
    QString testParms()
    {
        return testSettings;
    }
    FragmentHighlighter *fh;
    QWidget *lineNumberArea;
public slots:
    void insertText();

protected:
    void resizeEvent ( QResizeEvent *ev );

private slots:
    void updateLineNumberAreaWidth(int newBlockCount);
    void highlightCurrentLine();
    void updateLineNumberArea(const QRect &, int);
    void matchParentheses();

private:
    MainWindow* mainWindow;
    QString savedSettings;
    QString testSettings;
    bool matchLeftParenthesis ( QTextBlock currentBlock, int i, int numLeftParentheses );
    bool matchRightParenthesis ( QTextBlock currentBlock, int i, int numRightParentheses );
    void createParenthesisSelection(int pos);
};


class LineNumberArea : public QWidget
{
public:
    LineNumberArea ( TextEdit *editor ) : QWidget ( editor )
    {
        codeEditor = editor;
    }

    QSize sizeHint() const
    {
        return QSize(codeEditor->lineNumberAreaWidth(), 0);
    }

protected:
    void paintEvent ( QPaintEvent *event )
    {
        codeEditor->lineNumberAreaPaintEvent(event);
    }

private:
    TextEdit *codeEditor;
};
} // namespace GUI
} // namespace Fragmentarium
#endif // TEXTEDIT_H
