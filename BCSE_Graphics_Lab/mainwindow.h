#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <vector>
#include <Eigen/Dense>
QT_BEGIN_NAMESPACE
namespace Ui {
class MainWindow;
}
QT_END_NAMESPACE
struct Edge {
    int ymax;      // maximum y-value of the edge
    float x;       // x-coordinate of the lower end
    float slopeInv; // inverse slope of the edge
};
class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow();
protected:
    bool eventFilter(QObject *watched, QEvent *event) override;
private Q_SLOTS:
    void on_showAxis_clicked();

    void on_gridlines_clicked();


    void on_DDADraw_clicked();

    void on_store_clicked();

    void on_reset_clicked();

    void on_pushButton_clicked();


    void on_polar_circle_clicked();

    void on_bresen_circle_clicked();

    void on_PolarEllipse_clicked();

    void on_BresEllipse_clicked();
    void on_ZoomIn_clicked();

    void on_drawPoly_clicked();

    void on_floodfill_clicked();

    void on_Translation_clicked();

    void on_scanlineFill_clicked();

    void on_Rotate_clicked();

    void on_Shear_clicked();

    void on_gridOffset_valueChanged(int arg1);

    void on_ReflectX_clicked();

    void on_ReflectY_clicked();

    void on_scale_clicked();

    void on_boundaryFill_clicked();

private:
    int ymin;
    int ymax;
    int xmin;
    int xmax;
    Ui::MainWindow *ui;
    QPixmap temp;
    QVector<QPoint> polygonPoints;
    QSet<QPoint> visited;
    QHash<QPoint, QColor> pts;
    void draw_bresenham_ellipse(int, int, int, int);
    void draw_polar_ellipse(float, float, int, int);
    void plotPoint(int, int, int, int, int);
    void dfs(int a, int b, QVector<QPoint> & v);
    void createEdgeTable(QVector<QPoint>& polygonPoints, QVector<QVector<Edge>>& edgeTable, int& ymin, int& ymax);

    void draw_polar_circle(int ,int , int);
    void draw_bresenham_circle(int, int, int);
    void colorPoint(int x,int y,int r,int g, int b, int penwidth);
    void delay(int ms);
    void scanlineFill(QVector<QPoint>& polygonPoints, QRgb fillColor);
    void applyTranslation(Eigen::Matrix3d);
    void draw_dda_line(float x1,float y1,float x2,float y2);
    void draw_bresenham_line(int x1, int y1, int x2, int y2);
    void plotPoint(QPoint pt, QColor col);
    void iterativeFloodFill(int x, int y, QVector<QPoint>&);
    void plotPoints(QVector<QPoint>, QColor);
    void plotPoints(std::vector<QPoint>, QColor);
    void plotPointsTransform(QHash<QPoint, QColor>);
    // void applyTranslation(double tx, double ty);
    QPoint convertCoordinates(QPoint);
    QPoint convertCoordinates(int x, int y);
    QVector<QPoint> points;
    QPixmap justGridlines;
};
#endif // MAINWINDOW_H
