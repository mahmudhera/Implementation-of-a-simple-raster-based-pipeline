#include <iostream>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <cmath>
#include <stack>
#include <queue>
#include "bitmap_image.hpp"
using namespace std;

#define pi (2*acos(0.0))
#define epsilon (1.0e-6)

class homogeneous_point
{
public:
    double x, y, z, w;

    homogeneous_point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
        this->w = 1;
    }

    homogeneous_point() {
    }

    homogeneous_point(double x, double y, double z, double w)
    {
        assert (w != 0);
        this->x = x/w;
        this->y = y/w;
        this->z = z/w;
        this->w = 1;
    }

    homogeneous_point operator+ (const homogeneous_point& point)
    {
        double x = this->x + point.x;
        double y = this->y + point.y;
        double z = this->z + point.z;
        double w = this->w + point.w;
        homogeneous_point p(x, y, z, w);
        return p;
    }

    homogeneous_point operator- (const homogeneous_point& point)
    {
        double x = this->x - point.x;
        double y = this->y - point.y;
        double z = this->z - point.z;
        double w = this->w - point.w;
        homogeneous_point p(x, y, z, w);
    }

    void print()
    {
        cout << "Point: " << endl;
        cout << x << " " << y << " " << z << " " << w << endl;
    }

};

class Vector
{
public:
    double x, y, z;

    Vector(double x, double y, double z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    void normalize()
    {
        double r = sqrt(x*x + y*y + z*z);
        x = x / r;
        y = y / r;
        z = z / r;
    }

    Vector operator+(const Vector& v)
    {
        Vector v1(x+v.x, y+v.y, z+v.z);
        return v1;
    }

    Vector operator-(const Vector& v)
    {
        Vector v1(x-v.x, y-v.y, z-v.z);
        return v1;
    }

    Vector operator* (double m)
    {
        Vector v(x*m, y*m, z*m);
        return v;
    }

    static double dot(Vector a, Vector b)
    {
        return a.x*b.x + a.y*b.y + a.z*b.z;
    }

    static Vector cross(Vector a, Vector b)
    {
        Vector v(a.y*b.z - a.z*b.y, b.x*a.z - b.z*a.x, a.x*b.y - a.y*b.x);
        return v;
    }

    void print ()
    {
        cout << "Vector" << endl;
        cout << x << " " << y << " " << z << endl;
    }
};

class matrix
{
public:
    double values[4][4];
    int num_rows, num_cols;

    matrix(int rows, int cols)
    {
        assert (rows <= 4 && cols <= 4);
        num_rows = rows;
        num_cols = cols;
    }
    matrix(int n)
    {
        assert (n <= 4);
        num_rows = num_cols = n;
    }
    static matrix make_identity(int n)
    {
        assert (n <= 4);
        matrix m(n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                    m.values[i][j] = 1;
                else
                    m.values[i][j] = 0;
            }
        }
        return m;
    }

    void print()
    {
        cout << "Matrix:" << endl;
        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                cout << values[i][j] << "\t";
            }
            cout << endl;
        }
    }

    matrix operator+ (const matrix& m)
    {
        assert (this->num_rows == m.num_rows);
        assert (this->num_cols == m.num_cols);

        matrix m1(num_rows, num_cols);
        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                m1.values[i][j] = values[i][j] + m.values[i][j];
            }
        }
        return m1;
    }

    matrix operator- (const matrix& m)
    {
        assert (this->num_rows == m.num_rows);
        assert (this->num_cols == m.num_cols);

        matrix m1(num_rows, num_cols);
        for (int i = 0; i < num_rows; i++)
        {
            for (int j = 0; j < num_cols; j++)
            {
                m1.values[i][j] = values[i][j] - m.values[i][j];
            }
        }
        return m1;
    }

    matrix operator* (const matrix& m)
    {
        assert (this->num_cols == m.num_rows);
        matrix m1(this->num_rows, m.num_cols);

        for (int i = 0; i < m1.num_rows; i++) {
            for (int j = 0; j < m1.num_cols; j++) {
                double val = 0;
                for (int k = 0; k < this->num_cols; k++) {
                    val += this->values[i][k] * m.values[k][j];
                }
                m1.values[i][j] = val;
            }
        }
        return m1;
    }

    matrix operator* (double m)
    {
        matrix m1(this->num_rows, this->num_cols);
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                m1.values[i][j] = m * this->values[i][j];
            }
        }
        return m1;
    }

    homogeneous_point operator* (const homogeneous_point& p)
    {
        assert (this->num_rows == this->num_cols && this->num_rows == 4);

        matrix m(4, 1);
        m.values[0][0] = p.x;
        m.values[1][0] = p.y;
        m.values[2][0] = p.z;
        m.values[3][0] = p.w;

        matrix m1 = (*this)*m;
        homogeneous_point p1(m1.values[0][0], m1.values[1][0], m1.values[2][0], m1.values[3][0]);
        return p1;
    }

    matrix transpose()
    {
        matrix m(num_cols, num_rows);
        for (int i = 0; i < num_rows; i++) {
            for (int j = 0; j < num_cols; j++) {
                m.values[j][i] = values[i][j];
            }
        }
        return m;
    }

};

class color {
public:
    double r, g, b;
    color(double r, double g, double b) {
        this->r = r;
        this->g = g;
        this->b = b;
    }
    color() {
    }
};

matrix mat = matrix::make_identity(4);
stack <matrix> matrix_stack;
queue <color> color_list;
double eye_x, eye_y, eye_z;
double look_x, look_y, look_z;
double up_x, up_y, up_z;
double fov_x, fov_y, aspectRatio, near, far;
color backgroud;
int screen_x, screen_y;

matrix translation_matrix(double tx, double ty, double tz)
{
    matrix m = matrix::make_identity(4);
    m.values[0][3] = tx;
    m.values[1][3] = ty;
    m.values[2][3] = tz;
    return m;
}

matrix scale_matrix(double sx, double sy, double sz)
{
    matrix m = matrix::make_identity(4);
    m.values[0][0] = sx;
    m.values[1][1] = sy;
    m.values[2][2] = sz;
    return m;
}

matrix rotation_matrix(double angle, double ax, double ay, double az)
{
    Vector rotVector(ax, ay, az);
    rotVector.normalize();

    double cosine = cos(angle * pi / 180);
    double sine = sin(angle * pi / 180);

    Vector i(1, 0, 0), j(0, 1, 0), k(0, 0, 1);

    Vector new_i = i * cosine + rotVector * ( ( 1 - cosine ) * Vector::dot(i, rotVector) ) + Vector::cross(rotVector, i) * sine;
    Vector new_j = j * cosine + rotVector * ( ( 1 - cosine ) * Vector::dot(j, rotVector) ) + Vector::cross(rotVector, j) * sine;
    Vector new_k = k * cosine + rotVector * ( ( 1 - cosine ) * Vector::dot(k, rotVector) ) + Vector::cross(rotVector, k) * sine;

    matrix m = matrix::make_identity(4);

    m.values[0][0] = new_i.x;
    m.values[0][1] = new_i.y;
    m.values[0][2] = new_i.z;

    m.values[1][0] = new_j.x;
    m.values[1][1] = new_j.y;
    m.values[1][2] = new_j.z;

    m.values[2][0] = new_k.x;
    m.values[2][1] = new_k.y;
    m.values[2][2] = new_k.z;

    return m.transpose();

}

void final_stage_loop() {
    ifstream stage3;
    stage3.open("stage3.txt");

    color** pixels = new color*[screen_x];
    double** zs = new double*[screen_x];
    for (int i = 0; i < screen_x; i++) {
        pixels[i] = new color [screen_y];
        for (int j = 0; j < screen_y; j++) {
            pixels[i][j] = backgroud;
        }
        zs[i] = new double [screen_y];
        for (int j = 0; j < screen_y; j++) {
            zs[i][j] = +20;
        }
    }

    while (near < far)
    {
        color c = color_list.front();
        color_list.pop();
        double x, y, z;
        bool eof_found = false;
        homogeneous_point arr[3];
        for (int i = 0; i < 3; i++)
        {
            stage3 >> x;
            if (stage3.eof()) {
                eof_found = true;
                break;
            }
            stage3 >> y >> z;
            homogeneous_point p(x, y, z);
            arr[i] = p;
        }
        if (eof_found)
            break;
        for (int i = 0; i < 3; i++) {
            for (int j = i - 1; j >= 0; j--) {
                if (arr[j].y < arr[j+1].y) {
                    homogeneous_point temp = arr[j];
                    arr[j] = arr[j+1];
                    arr[j+1] = temp;
                }
            }
        }

        int index[3];
        for (int i = 0; i < 3; i++) {
            index[i] = screen_y / 2.0 * ( 1 - arr[i].y);
        }

        if (index[0] != index[1] && index[0] != index[2]) {
            for (int y = index[0] + 1; y <= index[1]; y++) {
                if (y < 0 || y >= screen_y) {
                    continue;
                }
                double y_scan = 1.0 - (2.0*y - 1) / screen_y;

                double x1 = (y_scan - arr[1].y) / (arr[0].y - arr[1].y) * (arr[0].x - arr[1].x) + arr[1].x;
                double z1 = (y_scan - arr[1].y) / (arr[0].y - arr[1].y) * (arr[0].z - arr[1].z) + arr[1].z;
                double x2 = (y_scan - arr[2].y) / (arr[0].y - arr[2].y) * (arr[0].x - arr[2].x) + arr[2].x;
                double z2 = (y_scan - arr[2].y) / (arr[0].y - arr[2].y) * (arr[0].z - arr[2].z) + arr[2].z;
                double x_l, x_r, z_l, z_r;
                if (x1 < x2) {
                    x_l = x1;
                    z_l = z1;
                    x_r = x2;
                    z_r = z2;
                } else {
                    x_r = x1;
                    z_r = z1;
                    x_l = x2;
                    z_l = z2;
                }

                int xindex0 = (x_l+1.0)*(screen_x/2.0);
                int xindex1 = (x_r+1.0)*(screen_x/2.0);

                for (int x = xindex0; x < xindex1; x++) {
                    if (x < 0 || x >= screen_x) {
                        continue;
                    }
                    double x_scan = (1.0 + 2 * x) / (double)screen_x - 1.0;
                    double z_scan = (x_scan - x_r) / (x_l - x_r) * (z_l - z_r) + z_r;
                    //if (z_scan - zs[x][y] < 1.0e-12 || -(z_scan - zs[x][y]) > 1.0e-12) {
                    if (z_scan - zs[x][y] <= epsilon && z_scan <= 1.0) {
                        pixels[x][y] = c;
                        zs[x][y] = z_scan;
                    }
                }
            }
        }
        if (index[2] != index[1] && index[0] != index[2]) {
            for (int y = index[1]+1; y <= index[2]; y++) {
                if (y < 0 || y >= screen_y) {
                    continue;
                }
                double y_scan = 1.0 - (2.0*y - 1) / screen_y;

                homogeneous_point* sorted = arr;
                double x1 = (y_scan - sorted[1].y) / (sorted[2].y - sorted[1].y) * (sorted[2].x - sorted[1].x) + sorted[1].x;
                double z1 = (y_scan - sorted[1].y) / (sorted[2].y - sorted[1].y) * (sorted[2].z - sorted[1].z) + sorted[1].z;
                double x2 = (y_scan - sorted[2].y) / (sorted[0].y - sorted[2].y) * (sorted[0].x - sorted[2].x) + sorted[2].x;
                double z2 = (y_scan - sorted[2].y) / (sorted[0].y - sorted[2].y) * (sorted[0].z - sorted[2].z) + sorted[2].z;
                double x_l, x_r, z_l, z_r;
                if (x1 < x2) {
                    x_l = x1;
                    z_l = z1;
                    x_r = x2;
                    z_r = z2;
                } else {
                    x_r = x1;
                    z_r = z1;
                    x_l = x2;
                    z_l = z2;
                }

                int xindex0 = (x_l+1)*(screen_x/2.0);
                int xindex1 = (x_r+1)*(screen_x/2.0);

                for (int x = xindex0; x < xindex1; x++) {
                    if (x < 0 || x >= screen_x) {
                        continue;
                    }
                    double x_scan = (1.0 + 2 * x) / (double)screen_x - 1.0;
                    double z_scan = (x_scan - x_r) / (x_l - x_r) * (z_l - z_r) + z_r;
                    //if (z_scan - zs[x][y] < 1.0e-12 || -(z_scan - zs[x][y]) > 1.0e-12) {
                    if (z_scan - zs[x][y] <= epsilon && z_scan <= 1.0) {
                        pixels[x][y] = c;
                        zs[x][y] = z_scan;
                    }
                }
            }
        }
    }


    bitmap_image image(screen_x, screen_y);
    for (int x = 0; x < screen_x; x++) {
        for (int y = 0; y < screen_y; y++) {
            image.set_pixel(x, y, pixels[x][y].r, pixels[x][y].g, pixels[x][y].b);
        }
    }
    image.save_image("out.bmp");

}

void stage3_main_loop()
{
    if (near == far) return;
    ifstream stage2;
    ofstream stage3;
    stage2.open ("stage2.txt");
    stage3.open ("temp.txt");
    stage3 << std::fixed;
    stage3 << std::setprecision(7);

    fov_x = fov_y * aspectRatio;
    double t = near * tan(pi / 180.0 * fov_y / 2.0);
    double r = near * tan(pi / 180.0 * fov_x / 2.0);

    matrix P = matrix::make_identity(4);

    P.values[0][0] = near / r;
    P.values[1][1] = near / t;
    P.values[2][2] = -(far + near)/(far - near);
    P.values[2][3] = -(2 * far * near)/(far - near);
    P.values[3][2] = -1;
    P.values[3][3] = 0;

    while (1)
    {
        double x, y, z;
        bool eof_found = false;
        homogeneous_point arr[3];
        for (int i = 0; i < 3; i++)
        {
            stage2 >> x;
            if (stage2.eof()) {
                eof_found = true;
                break;
            }
            stage2 >> y >> z;
            homogeneous_point p(x, y, z);
            arr[i] = p;
        }
        if (eof_found)
            break;

        color c = color_list.front();
        color_list.pop();

        bool flag[3];
        int counter = 0;
        for (int i = 0; i < 3; i++) {
            if (arr[i].z >= -near) {
                flag[i] = true;
                counter++;
            } else {
                flag[i] = false;
            }
        }

        if (counter == 1) {
            int k = -1;
            for (int i = 0; i < 3; i++) {
                if (flag[i] == true) {
                    k = i;
                    break;
                }
            }
            homogeneous_point points[3];
            for (int j = 0; j < 3; j++) {
                if (j == k) continue;
                double z = -near;
                double x = (z - arr[k].z) / (arr[j].z - arr[k].z) * (arr[j].x - arr[k].x) + arr[k].x;
                double y = (z - arr[k].z) / (arr[j].z - arr[k].z) * (arr[j].y - arr[k].y) + arr[k].y;
                homogeneous_point pt(x, y, z);
                points[j] = pt;
            }
            for (int i = 0; i < 3; i++) {
                if (i == k) continue;
                homogeneous_point p = arr[i];
                p = P*p;
                stage3 << p.x << " " << p.y << " " << p.z << endl;
            }
            for (int i = 0; i < 3; i++) {
                if (i == k) continue;
                homogeneous_point p = points[i];
                p = P*p;
                stage3 << p.x << " " << p.y << " " << p.z << endl;
                break;
            }
            stage3 << endl;
            for (int i = 0; i < 3; i++) {
                if (i == k) continue;
                homogeneous_point p = points[i];
                p = P*p;
                stage3 << p.x << " " << p.y << " " << p.z << endl;
            }
            for (int i = 2; i >= 0; i--) {
                if (i == k) continue;
                homogeneous_point p = arr[i];
                p = P*p;
                stage3 << p.x << " " << p.y << " " << p.z << endl;
                break;
            }
            stage3 << endl;
            color_list.push(c);
            color_list.push(c);
        } else if (counter == 2) {
            int k = -1;
            for (int i = 0; i < 3; i++) {
                if (flag[i] == false) {
                    k = i;
                    break;
                }
            }
            homogeneous_point points[2];
            int index = 0;
            for (int j = 0; j < 3; j++) {
                if (j == k) continue;
                double z = - near;
                double x = (z - arr[k].z) / (arr[j].z - arr[k].z) * (arr[j].x - arr[k].x) + arr[k].x;
                double y = (z - arr[k].z) / (arr[j].z - arr[k].z) * (arr[j].y - arr[k].y) + arr[k].y;
                homogeneous_point pt(x, y, z);
                points[index++] = pt;
            }
            index = 0;
            for (int i = 0; i < 3; i++) {
                homogeneous_point p = arr[i];
                if (i != k) p = points[index++];
                p = P*p;
                stage3 << p.x << " " << p.y << " " << p.z << endl;
            }
            stage3 << endl;
            color_list.push(c);
        } else if (counter == 0) {
            for (int i = 0; i < 3; i++) {
                homogeneous_point p = arr[i];
                p = P*p;
                stage3 << p.x << " " << p.y << " " << p.z << endl;
            }
            stage3 << endl;
            color_list.push(c);
        }

    }

    stage3.close();
    stage2.close();

}

void stage3_main_loop2()
{
    if (near == far) return;
    ifstream stage2;
    ofstream stage3;
    stage2.open ("temp.txt");
    stage3.open ("stage3.txt");
    stage3 << std::fixed;
    stage3 << std::setprecision(7);

    while (1)
    {
        double x, y, z;
        bool eof_found = false;
        homogeneous_point arr[3];
        for (int i = 0; i < 3; i++)
        {
            stage2 >> x;
            if (stage2.eof()) {
                eof_found = true;
                break;
            }
            stage2 >> y >> z;
            homogeneous_point p(x, y, z);
            arr[i] = p;
        }
        if (eof_found)
            break;

        color c = color_list.front();
        color_list.pop();

        bool flag[3];
        int counter = 0;
        for (int i = 0; i < 3; i++) {
            if (arr[i].z > 1.0) {
                flag[i] = true;
                counter++;
            } else {
                flag[i] = false;
            }
        }

        if (counter == 1) {
            int k = -1;
            for (int i = 0; i < 3; i++) {
                if (flag[i] == true) {
                    k = i;
                    break;
                }
            }
            homogeneous_point points[3];
            for (int j = 0; j < 3; j++) {
                if (j == k) continue;
                double z = 1.0;
                double x = (z - arr[k].z) / (arr[j].z - arr[k].z) * (arr[j].x - arr[k].x) + arr[k].x;
                double y = (z - arr[k].z) / (arr[j].z - arr[k].z) * (arr[j].y - arr[k].y) + arr[k].y;
                homogeneous_point pt(x, y, z);
                points[j] = pt;
            }
            for (int i = 0; i < 3; i++) {
                if (i == k) continue;
                homogeneous_point p = arr[i];
                stage3 << p.x << " " << p.y << " " << p.z << endl;
            }
            for (int i = 0; i < 3; i++) {
                if (i == k) continue;
                homogeneous_point p = points[i];
                stage3 << p.x << " " << p.y << " " << p.z << endl;
                break;
            }
            stage3 << endl;
            for (int i = 0; i < 3; i++) {
                if (i == k) continue;
                homogeneous_point p = points[i];
                stage3 << p.x << " " << p.y << " " << p.z << endl;
            }
            for (int i = 2; i >= 0; i--) {
                if (i == k) continue;
                homogeneous_point p = arr[i];
                stage3 << p.x << " " << p.y << " " << p.z << endl;
                break;
            }
            stage3 << endl;
            color_list.push(c);
            color_list.push(c);
        } else if (counter == 2) {
            int k = -1;
            for (int i = 0; i < 3; i++) {
                if (flag[i] == false) {
                    k = i;
                    break;
                }
            }
            homogeneous_point points[2];
            int index = 0;
            for (int j = 0; j < 3; j++) {
                if (j == k) continue;
                double z = 1.0;
                double x = (z - arr[k].z) / (arr[j].z - arr[k].z) * (arr[j].x - arr[k].x) + arr[k].x;
                double y = (z - arr[k].z) / (arr[j].z - arr[k].z) * (arr[j].y - arr[k].y) + arr[k].y;
                homogeneous_point pt(x, y, z);
                points[index++] = pt;
            }
            index = 0;
            for (int i = 0; i < 3; i++) {
                homogeneous_point p = arr[i];
                if (i != k) p = points[index++];
                stage3 << p.x << " " << p.y << " " << p.z << endl;
            }
            stage3 << endl;
            color_list.push(c);
        } else if (counter == 0) {
            for (int i = 0; i < 3; i++) {
                homogeneous_point p = arr[i];
                stage3 << p.x << " " << p.y << " " << p.z << endl;
            }
            stage3 << endl;
            color_list.push(c);
        }

    }

    stage3.close();
    stage2.close();

}

void stage2_main_loop()
{
    ifstream stage1;
    ofstream stage2;
    stage1.open ("stage1.txt");
    stage2.open ("stage2.txt");
    stage2 << std::fixed;
    stage2 << std::setprecision(7);

    Vector l(look_x - eye_x, look_y - eye_y, look_z - eye_z);
    l.normalize();
    Vector u(up_x, up_y, up_z);
    Vector r = Vector::cross(l, u);
    r.normalize();
    u = Vector::cross(r, l);

    matrix T = translation_matrix(-eye_x, -eye_y, -eye_z);
    matrix R = matrix::make_identity(4);

    R.values[0][0] = r.x;
    R.values[0][1] = r.y;
    R.values[0][2] = r.z;

    R.values[1][0] = u.x;
    R.values[1][1] = u.y;
    R.values[1][2] = u.z;

    R.values[2][0] = -l.x;
    R.values[2][1] = -l.y;
    R.values[2][2] = -l.z;

    matrix V = R*T;

    while (1)
    {
        double x, y, z;
        bool eof_found = false;
        for (int i = 0; i < 3; i++)
        {
            stage1 >> x;
            if (stage1.eof()) {
                eof_found = true;
                break;
            }
            stage1 >> y >> z;
            homogeneous_point p(x, y, z);
            p = V*p;
            stage2 << p.x << " " << p.y << " " << p.z << endl;
        }
        if (eof_found)
            break;
        stage2 << endl;
    }

    stage1.close();
    stage2.close();

}

void stage1_main_loop()
{
    ifstream scene;
    ofstream stage1;
    scene.open ("scene.txt");
    stage1.open ("stage1.txt");
    stage1 << std::fixed;
    stage1 << std::setprecision(7);

    string command;

    scene >> eye_x >> eye_y >> eye_z;
    scene >> look_x >> look_y >> look_z;
    scene >> up_x >> up_y >> up_z;
    scene >> fov_y >> aspectRatio >> near >> far;
    scene >> screen_x >> screen_y;
    scene >> backgroud.r >> backgroud.g >> backgroud.b;

    while (1)
    {
        scene >> command;
        if (command == "triangle")
        {
            color c;
            for (int i = 0; i < 3; i++)
            {
                double x, y, z;
                scene >> x >> y >> z;
                homogeneous_point p(x, y, z);
                p = mat*p;
                stage1 << p.x << " " << p.y << " " << p.z << endl;
            }
            stage1 << endl;
            scene >> c.r >> c.g >> c.b;
            color_list.push(c);
        }
        if (command == "translate")
        {
            double tx, ty, tz;
            scene >> tx >> ty >> tz;
            mat = mat * translation_matrix(tx, ty, tz);
        }
        if (command == "scale")
        {
            double sx, sy, sz;
            scene >> sx >> sy >> sz;
            mat = mat * scale_matrix(sx, sy, sz);
        }
        if (command == "rotate")
        {
            double ax, ay, az, angle;
            scene >> angle >> ax >> ay >> az;
            mat = mat * rotation_matrix(angle, ax, ay, az);
        }
        if (command == "push")
        {
            matrix_stack.push(mat);
        }
        if (command == "pop")
        {
            if (!matrix_stack.empty())
            {
                mat = matrix_stack.top();
                matrix_stack.pop();
            }
        }
        if (command == "end")
        {
            break;
        }
    }

    scene.close();
    stage1.close();

}

int main()
{
    cout << std::fixed;
    cout << std::setprecision(4);

    stage1_main_loop();
    stage2_main_loop();
    stage3_main_loop();
    stage3_main_loop2();
    final_stage_loop();

    return 0;
}
