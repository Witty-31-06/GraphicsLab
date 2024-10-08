#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <QPainter>
#include <QPixmap>
#include <QColor>
#include <QTimer>
#include <bits/stdc++.h>
#include <Eigen/Dense>
#include <QMouseEvent>
#define Delay delay(10)
#include <QElapsedTimer>
#define PI 3.14

// struct Edge {
//     int ymax;      // maximum y-value of the edge
//     float x;       // x-coordinate of the lower end
//     float slopeInv; // inverse slope of the edge
// };

MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    , ui(new Ui::MainWindow)
{
    ui->setupUi(this);
    ui->workArea->setMouseTracking(true);
    ui->workArea->installEventFilter(this);

    QPixmap canvas = ui->workArea->pixmap(Qt::ReturnByValue);
    if (canvas.isNull()) {
        canvas = QPixmap(ui->workArea->size().width(), ui->workArea->size().height());
        canvas.fill(Qt::white);
        ui->workArea->setPixmap(canvas);
    }
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::delay(int ms){
    QEventLoop loop;
    QTimer::singleShot(ms, &loop, &QEventLoop::quit);
    loop.exec();
}
QPoint MainWindow::convertCoordinates(int x, int y)
{
    int gridOffset = (ui->gridOffset->value()==0)?1:ui->gridOffset->value();
    int width = ui->workArea->width();
    int height = ui->workArea->height();
    int centerX=width/2;
    int centerY=height/2;
    int calcX = centerX+ x*gridOffset + gridOffset/2;
    int calcY = centerY -  y*gridOffset - gridOffset/2;
    return QPoint(calcX, calcY);

}
QPoint MainWindow::convertCoordinates(QPoint p)
{
    return convertCoordinates(p.x(), p.y());
}
void MainWindow::plotPointsTransform(QHash<QPoint, QColor> pointMap)
{
    int penwidth = ui->gridOffset->value();
    QPixmap canvas=justGridlines;
    QPainter painter(&canvas);


    for(auto it = pointMap.begin(); it != pointMap.end(); it++)
    {
        QColor col = it.value();
        QPoint pt = it.key();
        QPen pen=QPen(col,penwidth);
        painter.setPen(pen);
        painter.drawPoint(convertCoordinates(pt));
    }
    ui->workArea->setPixmap(canvas);
}
void MainWindow::plotPoints(QVector<QPoint> points, QColor col)
{
    int penwidth = ui->gridOffset->value();
    QPixmap canvas=ui->workArea->pixmap();
    QPainter painter(&canvas);
    QPen pen=QPen(col,penwidth);
    painter.setPen(pen);

    QVector<QPoint> transformedCoords;
    transformedCoords.reserve(points.size());

    for(const QPoint &pt: points)
    {
        pts[pt] = col;
        transformedCoords.push_back(convertCoordinates(pt));
    }
    painter.drawPoints(QPolygon(transformedCoords));

    ui->workArea->setPixmap(canvas);

}
void MainWindow::colorPoint(int x, int y, int r, int g, int b, int penwidth=1)
{

    QPixmap canvas=ui->workArea->pixmap();
    QPainter painter(&canvas);
    QPen pen=QPen(QColor(r,g,b),penwidth);
    painter.setPen(pen);
    painter.drawPoint(x, y);

    ui->workArea->setPixmap(canvas);
}

void MainWindow::plotPoint(int x, int y, int r, int g, int b)
{
    int gridOffset = (ui->gridOffset->value()==0)?1:ui->gridOffset->value();
    pts[QPoint(x, y)] =  QColor(r,g,b) ;
    QPoint coordinates = convertCoordinates(x, y);
    int calcX = coordinates.x();
    int calcY = coordinates.y();
    colorPoint(calcX, calcY, r,g,b, gridOffset);
}


void MainWindow::plotPoint(QPoint pt, QColor col)
{
    plotPoint(pt.x(), pt.y(), col.red(), col.green(), col.blue());
}
void MainWindow::on_showAxis_clicked() {
    int gridOffset = (ui->gridOffset->value() == 0) ? 1  : ui->gridOffset->value();
    int width = ui->workArea->width();
    int height = ui->workArea->height();
    int centerX = width / 2;
    int centerY = height / 2;
    int axisWidth=ui->gridOffset->value();
    // qDebug()<<width<<height<<centerX<<centerY<<axisWidth;
    // Draw horizontal axis
    for (int x = 0; x < width; ++x) {
        colorPoint(x, (centerY - gridOffset/2.0), 255, 0, 0, axisWidth); // Black color
    }

    // Draw vertical axis
    for (int y = 0; y < height; ++y) {
        colorPoint(centerX + gridOffset/2.0, y, 255, 0, 0, axisWidth); // Black color
    }

}
void MainWindow::on_gridlines_clicked() {

    QPixmap canvas = ui->workArea->pixmap(Qt::ReturnByValue);

    canvas = QPixmap(ui->workArea->size());
    canvas.fill(Qt::white);
    ui->workArea->setPixmap(canvas);

    int gridOffset = ui->gridOffset->value();
    int width = ui->workArea->width();
    int height = ui->workArea->height();
    if (gridOffset <= 0) return; // Prevent invalid grid offset

    int centerX = width / 2;
    int centerY = height / 2;

    int maxDistanceX = std::min(centerX, width - centerX);
    int maxDistanceY = std::min(centerY, height - centerY);
    qDebug()<<maxDistanceX<<maxDistanceY;
    QPainter painter(&canvas);
    if (gridOffset > 3)painter.setPen(QColor(128,128,128));
    else {painter.setPen(QColor(255,255,255));}

    for (int i = 0; i <= maxDistanceX; i += gridOffset)
    {

        QPoint qp1 = QPoint(centerX + i, 0);
        QPoint qp2 = QPoint(centerX + i, height);
        QPoint qp3 = QPoint(centerX - i, 0);
        QPoint qp4 = QPoint(centerX - i, height);



        painter.drawLine(qp1, qp2);
        painter.drawLine(qp3, qp4);

    }

    for(int i = 0; i<=maxDistanceY; i+=gridOffset)
    {
        QPoint qp5 = QPoint(0, centerY + i);
        QPoint qp6 = QPoint(width, centerY + i);
        QPoint qp7 = QPoint(0, centerY - i);
        QPoint qp8 = QPoint(width, centerY - i);
        painter.drawLine(qp5, qp6);
        painter.drawLine(qp7, qp8);
    }
    ui->workArea->setPixmap(canvas);

    if(pts.size() != 0)
    {
        QVector<QPoint> pastPoints = pts.keys();
        QColor col = pts.begin().value();
        plotPoints(pastPoints, col);
    }
    on_showAxis_clicked();
    justGridlines=canvas;
}



bool MainWindow::eventFilter(QObject *watched, QEvent *event) {
    if (watched == ui->workArea && event->type() == QEvent::MouseMove) {
        QMouseEvent *cursor = static_cast<QMouseEvent*>(event);
        int x = cursor->pos().x();
        int y = cursor->pos().y();
        int gridOffset = (ui->gridOffset->value()==0)?1:ui->gridOffset->value();
        int width = ui->workArea->width();
        int height = ui->workArea->height();
        int centerX=width/2;
        int centerY=height/2;
        ui->x_coordinate->setText(QString::number(floor((x-centerX)*1.0/gridOffset)));
        ui->y_coordinate->setText(QString::number(floor((centerY-y)*1.0/gridOffset)));
        return true; // Event handled
    }
    if(watched == ui->workArea && event->type() == QEvent::MouseButtonPress)
    {
        QMouseEvent *cursor = static_cast<QMouseEvent*>(event);
        int x = cursor->pos().x();
        int y = cursor->pos().y();
        int gridOffset = (ui->gridOffset->value()==0)?1:ui->gridOffset->value();
        int width = ui->workArea->width();
        int height = ui->workArea->height();
        int centerX=width/2;
        int centerY=height/2;
        int X = floor((x-centerX)*1.0/gridOffset);
        int Y = floor((centerY-y)*1.0/gridOffset);
        points.push_back({X, Y});
        polygonPoints.push_back(QPoint(X, Y));
        int calcX = centerX+ X*gridOffset + gridOffset/2;
        int calcY = centerY -  Y*gridOffset - gridOffset/2;
        colorPoint(calcX, calcY, 0,0, 255, gridOffset);

    }
    return QMainWindow::eventFilter(watched, event);

}

void MainWindow::draw_dda_line(float x1, float y1, float x2, float y2)
{
    float dx, dy, xinc, yinc, steps;
    // qDebug()<<x1<<y1<<x2<<y2;
    dx = x2 - x1;
    dy = y2 - y1;
    steps = std::max(abs(dx), abs(dy));  // Determine the number of steps based on the larger difference

    xinc = dx / steps;  // Increment in x
    yinc = dy / steps;  // Increment in y

    float x_float =  x1;
    float y_float =  y1;

    int xn = static_cast<int>(x_float);  // Initial x position in the grid
    int yn = static_cast<int>(y_float);  // Initial y position in the grid


    plotPoint(xn, yn, 165, 165, 0);

    for (int i = 1; i <= steps; i++)  // Loop to complete the straight line
    {
        x_float += xinc;
        y_float += yinc;

        int x_new = static_cast<int>(x_float);  // New x position in the grid
        int y_new = static_cast<int>(y_float);  // New y position in the grid

        if (x_new != xn || y_new != yn)  // If there is a change in the grid position
        {
            xn = x_new;
            yn = y_new;
            plotPoint(xn, yn, 165,165,0);// Color the new point
        }

        // qDebug() << x_new << y_new;  // Print the updated x, y values for debugging
    }
}


void MainWindow::draw_bresenham_line(int x1, int y1, int x2, int y2)
{
    int dx = abs(x2 - x1);
    int dy = abs(y2 - y1);
    int xsign = (x2 > x1) ? 1 : -1; // Determines if we move in positive or negative x direction
    int ysign = (y2 > y1) ? 1 : -1; // Determines if we move in positive or negative y direction

    bool isSteep = dy > dx;  // Check if the slope is steep

    // If the slope is steep, swap x and y coordinates to simplify drawing
    if (isSteep)
    {
        std::swap(x1, y1);
        std::swap(x2, y2);
        std::swap(dx, dy);
    }

    int p = 2 * dy - dx;  // Decision parameter
    int twoDy = 2 * dy;
    int twoDyDx = 2 * (dy - dx);

    int x = x1;
    int y = y1;

    // Plot the first point
    if (isSteep)
    {
        plotPoint(y, x, 0, 255, 255); // Swap x and y when slope is steep
    }
    else
    {
        plotPoint(x, y, 0, 255, 255); // Plot normally for shallow slopes
    }

    // Loop through the steps, plotting points
    for (int i = 0; i < dx; ++i)
    {
        x += xsign;  // Move in the x direction (either +1 or -1)
        if (p < 0)
        {
            p += twoDy;  // Adjust decision parameter
        }
        else
        {
            y += ysign;  // Move in the y direction (either +1 or -1 based on slope)
            p += twoDyDx;
        }

        // Plot the point
        if (isSteep)
        {
            plotPoint(y, x, 0, 255, 255);  // Swap x and y for steep slope
        }
        else
        {
            plotPoint(x, y, 0, 255, 255);  // Plot normally for shallow slope
        }
    }
}





void MainWindow::on_DDADraw_clicked()
{
    QElapsedTimer timer;
    timer.start();
    if(points.size()< 2 ) return;
    qint64 n  = points.size();
    auto coords1 = points[n-1];
    auto coords2 = points[n-2];
    draw_dda_line(coords1.x(), coords1.y(), coords2.x(), coords2.y());
    qint64 elapsedTime = timer.elapsed();
    ui->DDA_TIME->setText(QString::number(elapsedTime));
}


void MainWindow::on_store_clicked()
{
    temp = ui->workArea->pixmap(Qt::ReturnByValue);
}


void MainWindow::on_reset_clicked()
{
    ui->workArea->setPixmap(temp);
    points.clear();
}


void MainWindow::on_pushButton_clicked()
{
    QElapsedTimer timer;
    timer.start();
    if(points.size()< 2 ) return;
    int n = points.size();
    auto coords1 = points[n-1];
    auto coords2 = points[n-2];
    draw_bresenham_line(coords1.x(), coords1.y(), coords2.x(), coords2.y());
    qint64 elapsedTime = timer.elapsed();
    ui->Bresenham_TIME->setText(QString::number(elapsedTime));
}


void MainWindow::draw_polar_circle(int cx, int cy, int r)
{
    float delTheta = PI/(4*r);
    float angle = 0;
    while(angle <= PI/4)
    {
        float xn = round(cx + r*cos(angle));
        float yn = round(cy + r*sin(angle));
        plotPoint(xn, yn, 255,0,255);
        plotPoint(2*cx - xn, yn, 255,0,255);
        plotPoint(xn, 2*cy - yn, 255,0,255);
        plotPoint(2*cx - xn, 2*cy - yn, 255,0,255);
        float x_new = round(cx + r*sin(angle));
        float y_new = round(cy + r*cos(angle));
        plotPoint(x_new, y_new,255,0,255);
        plotPoint(2*cx - x_new, 2*cy-y_new,255,0,255);
        plotPoint(2*cx - x_new, y_new,255,0,255);
        plotPoint(x_new, 2*cy-y_new,255,0,255);
        angle += delTheta;
        // Delay;
    }
}


void MainWindow::on_polar_circle_clicked()
{
    if(points.size() == 0) return;
    auto coords = points[points.size() -1];
    int radius = ui->Radius->value();
    // qDebug()<<coords<<radius;
    QElapsedTimer timer;
    timer.start();
    draw_polar_circle(coords.x(), coords.y(), radius);
    qint64 time = timer.nsecsElapsed();
    ui->PolarTime->setText(QString::number(time) + QString("ns"));

}
void MainWindow::draw_bresenham_circle(int cx, int cy, int r)
{

    int p = (1 - r);
    int xn = cx;
    int yn = cy + r;
    plotPoint(cx, cy + r,0,255,255);
    plotPoint(cx, cy-r,0,255,255);
    plotPoint(cx-r, cy,0,255,255);
    plotPoint(cx+r, cy,0,255,255);
    xn++;
    while((xn-cx)< (yn - cy))
    {
        // qDebug() << (xn-cx) << (yn-cy);
        if(p< 0)
        {
            p = p + 2*(xn-cx) + 3;
        }

        else
        {
            yn--;
            p = p + 2*((xn-cx) - (yn-cy)) + 5;
        }
        plotPoint(xn, yn,0,255,255);
        plotPoint(-xn + 2*cx, -yn + 2*cy,0,255,255);
        plotPoint(xn, 2*cy - yn,0,255,255);
        plotPoint(2*cx-xn, yn,0,255,255);
        int x_new = yn - cy + cx;
        int y_new = xn - cx + cy;
        plotPoint(x_new, y_new,0,255,255);
        plotPoint(2*cx - x_new, 2*cy-y_new,0,255,255);
        plotPoint(2*cx - x_new, y_new,0,255,255);
        plotPoint(x_new, 2*cy-y_new,0,255,255);
        xn++;
        // Delay;
    }

}

void MainWindow::on_bresen_circle_clicked()
{
    if(points.size() == 0) return;
    auto coords = points[points.size() -1];
    int radius = ui->Radius->value();
    // qDebug()<<coords<<radius;
    QElapsedTimer timer;
    timer.start();
    draw_bresenham_circle(coords.x(), coords.y(), radius);
    qint64 time = timer.nsecsElapsed();
    ui->BresenHamTime->setText(QString::number(time)+QString("ns"));
}
void MainWindow::draw_polar_ellipse(float rx, float ry, int cx, int cy)
{
    QElapsedTimer t = QElapsedTimer();
    t.start();
    float theta = 1 /(std::max(rx, ry));
    float angle = 0;
    while(angle <= PI/2)
    {
        float x = round(rx*cos(angle));
        float y = round(ry*sin(angle));
        plotPoint(cx+x, cy+y, 124,231,0);
        plotPoint(cx-x, cy+y, 124,231,0);
        plotPoint(cx+x, cy-y, 124,231,0);
        plotPoint(cx-x, cy-y, 124,231,0);
        angle += theta;
    }
    qint64 time = t.nsecsElapsed();
    ui->PolEllipseTime->setText(QString::number(time));
}
void MainWindow::on_PolarEllipse_clicked()
{
    if(points.size() == 0) return;
    auto coords = points[points.size() -1];
    int rx = ui->Rx->value();
    int ry = ui->Ry->value();
    draw_polar_ellipse( rx,  ry, coords.x(), coords.y());
}

void MainWindow::draw_bresenham_ellipse(int rx, int ry, int cx, int cy)
{
    QElapsedTimer t = QElapsedTimer();
    t.start();
    int a = rx; // Major radius
    int b = ry; // Minor radius
    int x = 0;
    int y = b;
    int a2 = a * a;
    int b2 = b * b;
    int d1 = b2 - a2 * b + 0.25 * a2;
    int dx = 2 * b2 * x;
    int dy = 2 * a2 * y;

    // Plot the initial points
    plotPoint(cx + x, cy + y, 235, 124, 64);
    plotPoint(cx - x, cy + y, 235, 124, 64);
    plotPoint(cx + x, cy - y, 235, 124, 64);
    plotPoint(cx - x, cy - y, 235, 124, 64);

    // Region 1
    while (dx < dy) {
        x++;
        dx += 2 * b2;
        if (d1 < 0) {
            d1 += b2 + dx;
        } else {
            y--;
            dy -= 2 * a2;
            d1 += b2 + dx - dy;
        }
        // Plot the points
        plotPoint(cx + x, cy + y, 235, 124, 64);
        plotPoint(cx - x, cy + y, 235, 124, 64);
        plotPoint(cx + x, cy - y, 235, 124, 64);
        plotPoint(cx - x, cy - y, 235, 124, 64);
    }

    // Region 2
    int d2 = b2 * (x + 0.5) * (x + 0.5) + a2 * (y - 1) * (y - 1) - a2 * b2;
    while (y > 0) {
        y--;
        dy -= 2 * a2;
        if (d2 > 0) {
            d2 += a2 - dy;
        } else {
            x++;
            dx += 2 * b2;
            d2 += a2 - dy + dx;
        }
        // Plot the points
        plotPoint(cx + x, cy + y, 235, 124, 64);
        plotPoint(cx - x, cy + y, 235, 124, 64);
        plotPoint(cx + x, cy - y, 235, 124, 64);
        plotPoint(cx - x, cy - y, 235, 124, 64);
    }
    qint64 time = t.nsecsElapsed();
    ui->BresenhaEllipseTime->setText(QString::number(time));
}

void MainWindow::on_BresEllipse_clicked()
{
    if(points.size() == 0) return;
    auto coords = points[points.size() -1];
    int rx = ui->Rx->value();
    int ry = ui->Ry->value();
    draw_bresenham_ellipse(rx, ry, coords.x(), coords.y());
}



void MainWindow::on_ZoomIn_clicked()
{

}


void MainWindow::on_drawPoly_clicked()
{
    qDebug()<<"HAHA";
    int n = polygonPoints.size();
    if (n < 3) return;
    for(int i = 0; i<n-1; i++)
    {
        QPoint pt1 = polygonPoints[i];
        QPoint pt2 = polygonPoints[i+1];
        int x1 = pt1.x();
        int y1 = pt1.y();
        int x2 = pt2.x();
        int y2 = pt2.y();
        draw_dda_line(x1,y1,x2,y2);
    }
    QPoint pt1 = polygonPoints[0];
    QPoint pt2 = polygonPoints[n-1];
    int x1 = pt1.x();
    int y1 = pt1.y();
    int x2 = pt2.x();
    int y2 = pt2.y();
    draw_bresenham_line(x1,y1,x2,y2);
    // polygonPoints.clear();
}
void MainWindow::iterativeFloodFill(int x, int y, QVector<QPoint> &interior)
{
    std::queue<QPoint> q;
    q.push(QPoint(x, y));

    int xdir[] = {0, -1, 1, 0};
    int ydir[] = {-1, 0, 0, 1};

    while (!q.empty())
    {
        QPoint p = q.front();
        q.pop();

        if (visited.contains(p) || pts.contains(p))
            continue;

        interior.push_back(QPoint(p.x(), p.y()));
        visited.insert(p);

        for (int i = 0; i < 4; ++i)
        {
            QPoint neighbor = QPoint(p.x() + xdir[i], p.y() + ydir[i]);
            if (!visited.contains(neighbor) && !pts.contains(neighbor))
            {
                q.push(neighbor);
            }
        }
    }
}


void MainWindow::on_floodfill_clicked()
{
    if (polygonPoints.size() == 0) return;
    QPoint seed = polygonPoints[0];
    QVector<QPoint> interior_points;
    interior_points.reserve(1011);
    polygonPoints.clear();
    iterativeFloodFill(seed.x(), seed.y(), interior_points);
    plotPoints(interior_points, QColor(10,20,30));
}
void MainWindow::applyTranslation(Eigen::Matrix3d matrix) {
    QHash<QPoint, QColor> newPoints;
    Eigen::Vector3d vector;
    for(auto it = pts.begin(); it != pts.end(); it++)
    {
        QPoint p = it.key();
        QColor col = it.value();
        vector << p.x(), p.y(), 1;
        Eigen::Vector3d result = matrix * vector;
        newPoints[QPoint(result(0), result(1))] = col;
    }
    pts = newPoints;
    plotPointsTransform(newPoints);
}
void MainWindow::on_Translation_clicked()
{

    double tx = ui->transX->value();
    double ty = ui->transY->value();
    points.clear();
    Eigen::Matrix3d matrix;
    matrix <<   1, 0, tx,
                0, 1, ty,
                0,0,1;
    applyTranslation(matrix);

}
void MainWindow::createEdgeTable(QVector<QPoint>& polygonPoints,QVector<QVector<Edge>>& edgeTable, int& ymin, int& ymax)
{
    int n = polygonPoints.size();
    for (int i = 0; i < n; i++) {
        QPoint p1 = polygonPoints[i];
        QPoint p2 = polygonPoints[(i + 1) % n];  // Wrapping to connect the last point with the first one

        if (p1.y() == p2.y()) continue;  // Ignore horizontal edges

        // Ensure p1.y < p2.y (to sort edges correctly)
        if (p1.y() > p2.y()) std::swap(p1, p2);

        ymin = std::min(ymin, p1.y());
        ymax = std::max(ymax, p2.y());

        // Compute inverse slope (to increment x for each scanline)
        float slopeInv = (float)(p2.x() - p1.x()) / (p2.y() - p1.y());

        // Add edge to the edge table, using the shifted index
        Edge edge = {p2.y(), (float)p1.x(), slopeInv};
        edgeTable[p1.y()].push_back(edge);  // Offset y by ymin to ensure non-negative indexing
    }
}


void MainWindow::scanlineFill(QVector<QPoint>& polygonPoints, QRgb fillColor)
{
    // Determine ymin and ymax of the polygon
    int ymin = INT_MAX, ymax = INT_MIN;

    // Determine the bounds of the polygon
    for (const QPoint& point : polygonPoints) {
        ymin = std::min(ymin, point.y());
        ymax = std::max(ymax, point.y());
    }

    // Create an edge table with enough room for both negative and positive y-values
    QVector<QVector<Edge>> edgeTable(ymax - ymin + 1);  // Ensure space for all y-values (shifted by ymin)

    // Populate the edge table
    createEdgeTable(polygonPoints, edgeTable, ymin, ymax);

    // Active Edge List (AEL)
    QVector<Edge> activeEdgeList;

    // Loop over each scanline from ymin to ymax
    for (int y = ymin; y <= ymax; y++) {
        // Step 1: Add edges from edgeTable[y - ymin] to activeEdgeList where y == ymin of the edge
        activeEdgeList.append(edgeTable[y - ymin]);

        // Step 2: Remove edges from activeEdgeList where ymax == current y
        for (int i = 0; i < activeEdgeList.size(); i++) {
            if (activeEdgeList[i].ymax == y) {
                activeEdgeList.remove(i);
                i--;
            }
        }

        // Step 3: Sort the activeEdgeList by x value
        std::sort(activeEdgeList.begin(), activeEdgeList.end(), [](Edge e1, Edge e2) {
            return e1.x < e2.x;
        });

        // Step 4: Fill between pairs of intersections
        for (int i = 0; i < activeEdgeList.size(); i += 2) {
            if (i + 1 >= activeEdgeList.size()) break;

            int x1 = (int)activeEdgeList[i].x;
            int x2 = (int)activeEdgeList[i + 1].x;

            // Draw a horizontal line between x1 and x2 using plotPoint
            for (int x = x1; x < x2; x++) {
                if(pts.contains(QPoint(x,y))) continue;
                plotPoint(x, y, qRed(fillColor), qGreen(fillColor), qBlue(fillColor));
                Delay;
            }
        }

        // Step 5: Update x for each edge in the activeEdgeList (x += slopeInv)
        for (int i = 0; i < activeEdgeList.size(); i++) {
            activeEdgeList[i].x += activeEdgeList[i].slopeInv;
        }
    }
}
void MainWindow::on_scanlineFill_clicked()
{
    scanlineFill(polygonPoints,qRgb(20,200,10));
}


void MainWindow::on_Rotate_clicked()
{
    double theta = ui->rotAngle->value();
    theta = (theta*PI)/180;
    points.clear();
    Eigen::Matrix3d matrix;
    matrix << std::cos(theta), std::sin(theta), 0,
        -std::sin(theta),  std::cos(theta), 0,
        0, 0, 1;
    applyTranslation(matrix);
}


void MainWindow::on_Shear_clicked()
{
    double shx = ui->shx->value();
    double shy = ui->shy->value();
    // points.clear();
    Eigen::Matrix3d shearMatrix;
    shearMatrix << 1, shx, 0,
        shy, 1, 0,
        0, 0, 1;

    QVector<QPoint> newPoints;
    Eigen::Vector3d vector;
    for(QPoint p : points)
    {
        qDebug()<<p;
        vector << p.x(), p.y(), 1;
        Eigen::Vector3d result = shearMatrix * vector;
        newPoints.push_back(QPoint(result(0), result(1)));

    }
    ui->workArea->setPixmap(justGridlines);
    polygonPoints = newPoints;
    on_drawPoly_clicked();
    // pts = newPoints;
    // plotPointsTransform(newPoints);
}


void MainWindow::on_gridOffset_valueChanged(int arg1)
{
    return;
}


void MainWindow::on_ReflectX_clicked()
{
    Eigen::Matrix3d reflectXMatrix;
    reflectXMatrix << 1, 0, 0,
        0, -1, 0,
        0, 0, 1;
    applyTranslation(reflectXMatrix);
}


void MainWindow::on_ReflectY_clicked()
{
    Eigen::Matrix3d reflectYMatrix;
    reflectYMatrix << -1, 0, 0,
        0, 1, 0,
        0, 0, 1;
    applyTranslation(reflectYMatrix);
}


void MainWindow::on_scale_clicked()
{
    double sx = ui->sx->value();
    double sy = ui->sy->value();
    Eigen::Matrix3d scalingMatrix;
    scalingMatrix << sx, 0, 0,
        0, sy, 0,
        0, 0, 1;
    applyTranslation(scalingMatrix);
}


void MainWindow::on_boundaryFill_clicked()
{

}

