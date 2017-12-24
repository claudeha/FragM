/*
 * TimeLineDialog : gui representation of easingcurves for Fragmentarium
 * Copyright (C) 2016  R P <digilanti@hotmail.com>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef TIMELINE_H
#define TIMELINE_H

#include <QColor>
#include <QtGui>
#include <QtWidgets/QDialog>
#include <QtWidgets/QDialogButtonBox>
#include <QtWidgets/QGraphicsView>
#include <QtWidgets/QGraphicsScene>
#include <QtWidgets/QGraphicsItem>
#include <QtWidgets/QGraphicsTextItem>
#include <QtWidgets/QGraphicsRectItem>
#include <QtWidgets/QGraphicsEllipseItem>
#include <QtWidgets/QGraphicsItemGroup>

#include "ui_TimeLineDialog.h"
#include "TextEdit.h"

class MainWindow;

namespace Fragmentarium {
  namespace GUI {
    
//   class EasingItem : public QGraphicsItem
//   {
//   public:
//       EasingItem(EasingInfo ei) : data(ei);
//       ~EasingItem();
//       QRectF boundingRect() const
//       {
//           qreal penWidth = 1;
//           return QRectF(-10 - penWidth / 2, -10 - penWidth / 2,
//                         20 + penWidth, 20 + penWidth);
//       }
// 
//       void paint(QPainter *painter, const QStyleOptionGraphicsItem *option,
//                  QWidget *widget)
//       {
//           painter->drawRoundedRect(-10, -10, 20, 20, 5, 5);
//       }
//       EasingInfo data;
//   };

 class TimeLineDialog : public QDialog
    {
      Q_OBJECT
    public:
      TimeLineDialog(Fragmentarium::GUI::MainWindow* parent);
      ~TimeLineDialog();
    public slots:
      
    protected slots:
        void customContextMenuRequested(QPoint);
        
    private slots:
      void itemChange(const QList<QRectF> & );
      void readTimeLineSettings();
      void saveTimeLineSettings();
      void restoreTimeLineSettings();
      void createKeyframeMap();
      void createEasingCurveMap();
      QPainterPath createCurve(QSize sz, int type);
    private:
      Ui::TimeLineDialog m_ui;
      MainWindow* mainWin;
      QBrush greenBrush;
      QBrush grayBrush;
      QBrush redBrush;
      QPen outlinePen;
      int frames;
      int keyframeCount;
      QStringList uNames;
      QStringList uSettings;
      int vCount;
      int yOff;
      double timeMax;
      int renderFPS;
      // TimeLineScene
      QGraphicsScene *scene;
      QGraphicsTextItem *text;
      QMap<int, QGraphicsTextItem *> textMap;
      QMap<int, QGraphicsRectItem*> rectMap;
      QMap<int, QGraphicsPathItem*> pathMap;
      QMap<int, QGraphicsItemGroup*> eGroupMap;
      QMap<int, KeyFrameInfo*> keyframeMap;
      QMap<int, EasingInfo*> easingMap;
      QRectF sceneMaxRect;
    };
    
  }
}

#endif
