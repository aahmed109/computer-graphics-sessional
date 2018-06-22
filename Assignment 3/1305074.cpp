#include<stdio.h>
#include<fstream>
#include<string>
#include<iostream>
#include<sstream>
#include<vector>
#include<algorithm>
#include<cmath>
#include<iomanip>
#include<random>
#include "bitmap_image.hpp"

#define PI 3.14159265

using namespace std;

bool comp (double i, double j) { return (i > j); }
bool revcomp (double i, double j) { return (i < j); }
/*struct comp {
  bool operator() (double i, double j) { return (i > j);}
} mycomp;
*/
ofstream fout1("test.txt");

typedef vector<double> Point;
typedef vector<unsigned char> Color;

int Screen_Width;
int Screen_Height;
double x_left_limit;
double x_right_limit;
double y_top_limit;
double y_bottom_limit;
double z_front_limit;
double z_rear_limit;

double z_max;

double dx;
double dy;

double Top_y;
double Left_x;

double** Z_buffer;


class Equations{
public:
    vector<double> xy;
    vector<double> yz;
    vector<double> zx;
};

class Triangle{
private:
    double find_gcd(double a, double b)
    {
        if (a < b){

            return find_gcd(b, a);
        }

        if (fabs(b) < 0.001){
            return a;
        }

        else{
            return (find_gcd(b, a - floor(a / b) * b));
        }
    }

    Point GCD(Point A)
    {
        double s = A[0];
        Point B = A;
        for(int i = 1; i < 3; i++){
            //cout<<"mew"<<endl;
            //printf("%f %f\n", A[i], s);
            if(s != 0.0000000 && A[i] != 0.0000000) s = find_gcd(fabs(A[i]), fabs(s));
        }
        //printf("%f\n", s);
        if(s != 0.0000000){
            for(int i = 0; i < 3; i++){
                B[i] /= s;
            }
        }
        bool neg = true;
        for(int i = 0; i < 3; i++){
            if(B[i] > 0.0000000)    neg = false;
        }
        if(neg){
            for(int i = 0; i < 3; i++){
                B[i] = -B[i];
            }
        }
        return B;
    }

    Equations line_to_eqn(Point A, Point B)
    {
        Equations eqn;
        Point C;    ///Actually Vector, bit of confusion in here
        C.push_back(B[0]-A[0]);
        C.push_back(B[1]-A[1]);
        C.push_back(B[2]-A[2]);
        Point F = GCD(C);
        //cout<<"printing C"<<endl;
        //for(int j = 0; j < 3; j++) printf("%f ", C[j]);
        //cout<<endl;
        //cout<<"printing F"<<endl;
        //for(int j = 0; j < 3; j++) printf("%f ", F[j]);
        //cout<<endl;
        eqn.xy.push_back(F[1]);
        eqn.xy.push_back(-F[0]);
        eqn.xy.push_back(-F[0] * A[1] + F[1] * A[0]);

        eqn.yz.push_back(F[2]);
        eqn.yz.push_back(-F[1]);
        eqn.yz.push_back(-F[1] * A[2] + F[2] * A[1]);

        eqn.zx.push_back(F[0]);
        eqn.zx.push_back(-F[2]);
        eqn.zx.push_back(-F[2] * A[0] + F[0] * A[2]);

        return eqn;
        /*///the following loop will be more simple, what a dumb you are!
        for(int j = 0; j < 3; j++){
            if(j == 0){
                eqn.xy.push_back(F[1]);
                eqn.xy.push_back(-F[0]);
                eqn.xy.push_back(F[0] * A[1] - F[1] * A[0]);

                eqn.yz.push_back(-9999.0);

                eqn.zx.push_back(F[0]);
                eqn.zx.push_back(-F[2]);
                eqn.zx.push_back(F[2] * A[0] - F[0] * A[2]);
            }

            else if(j == 1){
                eqn.xy.push_back(F[1]);
                eqn.xy.push_back(-F[0]);
                eqn.xy.push_back(F[0] * A[1] - F[1] * A[0]);

                eqn.yz.push_back(F[2]);
                eqn.yz.push_back(-F[1]);
                eqn.yz.push_back(F[1] * A[2] - F[2] * A[1]);

                eqn.zx.push_back(-9999.0);
            }

            else if(j == 2){
                eqn.xy.push_back(-9999.0);

                eqn.yz.push_back(F[2]);
                eqn.yz.push_back(-F[1]);
                eqn.yz.push_back(F[1] * A[2] - F[2] * A[1]);

                eqn.zx.push_back(F[0]);
                eqn.zx.push_back(-F[2]);
                eqn.zx.push_back(F[2] * A[0] - F[0] * A[2]);
            }
        }*/
    }

public:
    vector<Point> theTriangle;
    Color color;
    vector<Equations> eqtn;
    Triangle(){
    }
    void generate_color(){
        for(int i = 0; i < 3; i++)  color.push_back((unsigned char)(rand() % 255));
    }

    void generate_equation()
    {

        for(int i = 0; i < 2; i++){
            Point A = theTriangle[i];
            //cout<<A[0]<<endl;
            Point B = theTriangle[i+1];
            //cout<<B[0]<<endl;
            eqtn.push_back(line_to_eqn(A, B));
        }

        eqtn.push_back(line_to_eqn(theTriangle[2], theTriangle[0]));
        /*for(int i = 0; i < 2; i++){
            cout<<"wut " << i+1<<endl;
            Point A = triangle.theTriangle[i];
            Point B = triangle.theTriangle[i+1];
            Point C;    //Actually Vector
            C.push_back(B[0]-A[0]);
            C.push_back(B[1]-A[1]);
            C.push_back(B[2]-A[2]);
            Point F = GCD(C);
            cout<<"printing C"<<endl;
            for(int j = 0; j < 3; j++) printf("%f ", C[j]);
            cout<<endl;
            cout<<"printing F"<<endl;
            for(int j = 0; j < 3; j++) printf("%f ", F[j]);
            cout<<endl;
            ///nicher part ta k function banaite hobe
            for(int j = 0; j < 3; j++){
                if(F[j] == 0.0000000){
                    if(j == 0)  triangle.eqn.x = A[0];
                    else if(j == 1) triangle.eqn.y = A[1];
                    else if(j == 2) triangle.eqn.z = A[2];
                }
                else if(F[j] == 1.0000000){
                    if(j == 0){
                        if(F[1] > 0.0000000){
                            triangle.eqn.xy.push_back(F[1]);
                            triangle.eqn.xy.push_back(-1.0000000);
                            triangle.eqn.xy.push_back(A[1] - F[1] * A[0]);
                        }
                        if(F[2] > 0.0000000){
                            triangle.eqn.zx.push_back(F[2]);
                            triangle.eqn.zx.push_back(-1.0000000);
                            triangle.eqn.zx.push_back(A[2] - F[2] * A[0]);
                        }
                    }
                    else if(j == 1){
                        if(F[0] > 0.0000000){
                            triangle.eqn.xy.push_back(1.0000000);
                            triangle.eqn.xy.push_back(-F[0]);
                            triangle.eqn.xy.push_back(F[0] * A[1] - A[0]);
                        }
                        if(F[2] > 0.0000000){
                            triangle.eqn.yz.push_back(F[2]);
                            triangle.eqn.yz.push_back(-1.0000000);
                            triangle.eqn.yz.push_back(A[2] - F[2] * A[1]);
                        }
                    }
                    else if(j == 2){
                        if(F[1] > 0.0000000){
                            triangle.eqn.yz.push_back(1.0000000);
                            triangle.eqn.yz.push_back(-F[1]);
                            triangle.eqn.yz.push_back(F[1] * A[2] - A[1]);
                        }
                        if(F[0] > 0.0000000){
                            triangle.eqn.zx.push_back(F[0]);
                            triangle.eqn.zx.push_back(-1.0000000);
                            triangle.eqn.zx.push_back(A[0] - F[0] * A[2]);
                        }
                    }
                }
            }
        }*/
    }

    ~Triangle(){
        //delete theTriangle;
    }
};

void read_config_file(string configs[])
{
    ifstream fin("config.txt");
    string str;

    int lc = 0;
    while(getline(fin, str)) configs[lc++] = str;
    fin.close();

    for(int i = 0; i < 4; i++){
        stringstream ss(configs[i]);

        if(i == 0) {
            ss>>Screen_Width;
            ss>>Screen_Height;
        }
        else if(i == 1) {
            ss>>x_left_limit;
            x_right_limit = -x_left_limit;
        }
        else if(i == 2) {
            ss>>y_bottom_limit;
            y_top_limit = -y_bottom_limit;
        }
        else{
            ss>>z_front_limit;
            ss>>z_rear_limit;
        }
    }
}

void set_Z_buffer()
{
    Z_buffer = new double*[Screen_Height];
    for(int i = 0; i < Screen_Height; i++){
        Z_buffer[i] = new double[Screen_Width];
    }

    for(int i = 0; i < Screen_Height; i++){
        for(int j = 0; j < Screen_Width; j++){
            Z_buffer[i][j] = z_max;
        }
    }
}

void read_stage_3_file(string triangles[])
{
    ifstream fin("stage3.txt");
    string str;

    int lc = 0;
    while(getline(fin, str)){
        if(str.compare("")) triangles[lc++] = str;
    }

    fin.close();
}

void find_triangles(Triangle tri[], Equations eqn[], string triangles[], int n)
{
    Point point;
    for(int i = 0; i < n; i++){
        for(int j = 0; j < 3; j++){
            stringstream ss(triangles[i * 3 + j]);
            double s;
            for(int k = 0; k < 3; k++){
                ss>>s;
                point.push_back(s);
            }
            tri[i].theTriangle.push_back(point);
            point.clear();
        }
        tri[i].generate_color();

        tri[i].generate_equation();

        //cout<<"printing equations"<<endl;
        //for(int k = 0; k < 3; k++){
            //Equations eqnt = tri[i].eqtn[k];
            //for(int j = 0; j < 3; j++) printf("%f ", eqnt.xy[j]);
            //cout<<endl;
            //for(int j = 0; j < 3; j++) printf("%f ", eqnt.yz[j]);
            //cout<<endl;
            //for(int j = 0; j < 3; j++) printf("%f ", eqnt.zx[j]);
            //cout<<endl;
            //cout<<endl;
        //}
        /*cout<<endl;
        for(int k = 0; k < 3; k++){
            for(int j = 0; j < 3; j++) printf("%f ", tri[i].eqtn[k].yz[j]);
            cout<<endl;
        }
        cout<<endl;
        for(int k = 0; k < 3; k++){
            for(int j = 0; j < 3; j++) printf("%f ", tri[i].eqtn[k].zx[j]);
            cout<<endl;
        }
        cout<<endl;*/
    }
}

void print_triangle(Triangle tri[], int n)
{
        for(int i = 0; i < n; i++){
        for(int j = 0; j < 3; j++){
            for(int k = 0; k < 3; k++){
                printf("%f ", tri[i].theTriangle[j][k]);
            }
            cout<<endl;
        }
        cout<<endl;

        for(int j = 0; j < 3; j++){
            printf("%u ", tri[i].color[j]);
        }
        cout<<endl;
        cout<<endl;
    }
}

Point intersection_check(Triangle one, int n, double scan_line)
{
    Point intersect;
    Equations eqt = one.eqtn[n];
    //Point intersect;
    if((eqt.xy[0] == 0.000000 || eqt.xy[0] == -0.000000) && (eqt.xy[1] != 0.000000 || eqt.xy[1] != -0.000000)){
        double y = eqt.xy[2] / eqt.xy[1];
        if(y == scan_line){
            intersect.push_back(-0.000000);
            intersect.push_back(scan_line);
            if(eqt.zx[0] != 0.000000 || eqt.zx[0] != -0.000000) intersect.push_back(eqt.zx[2] / eqt.zx[0]);
            else intersect.push_back(-0.000000);
        }
    }

    else if((eqt.xy[1] == 0.000000 || eqt.xy[1] == -0.000000) && (eqt.xy[0] != 0.000000 || eqt.xy[0] != -0.000000)){
        if(scan_line == 0.000000){
            double x = eqt.xy[2] / eqt.xy[0];
            intersect.push_back(x);
            intersect.push_back(scan_line);
            if(eqt.yz[1] != 0.000000 || eqt.yz[1] != -0.000000) intersect.push_back(eqt.yz[2] / eqt.yz[1]);
            else intersect.push_back(0.000000);
        }
    }

    else if((eqt.xy[0] == 0.000000 || eqt.xy[0] == -0.000000) && (eqt.xy[1] == 0.000000 || eqt.xy[1] == -0.000000)){
        if(scan_line == 0.000000){
            intersect.push_back(-0.000000);
            intersect.push_back(-0.000000);
            if(eqt.yz[1] != 0.000000 || eqt.yz[1] != -0.000000) intersect.push_back(eqt.yz[2] / eqt.yz[1]);
            else intersect.push_back(-0.000000);
        }
    }

    else{
        //if()
        intersect.push_back((eqt.xy[2] - eqt.xy[1] * scan_line) / eqt.xy[0]);
        intersect.push_back(scan_line);
        if(eqt.yz[1] != 0.000000 || eqt.yz[1] != -0.000000) intersect.push_back((eqt.yz[2] - eqt.yz[0] * scan_line) / eqt.yz[1]);
        else intersect.push_back(-0.000000);
    }

    return intersect;
}

void apply_Z_buffer_algo(bitmap_image *image, Triangle tri[], int n)
{
    for(int i = 0; i < n; i++){
        //cout<<"for the "<<i+1<<"-th time"<<endl;
        Triangle one = tri[i];
        //triangleLines(one);
        vector<double> Ys, Xs;
        for(int j = 0; j < 3; j++){
            Xs.push_back(one.theTriangle[j][0]);
            Ys.push_back(one.theTriangle[j][1]);
        }
        sort(Ys.begin(), Ys.end(), comp);
        sort(Xs.begin(), Xs.end(), revcomp);
        double top_scan_line = -9999.0;
        double bottom_scan_line = 9999.0;
        for(int j = 0; j < 3; j++){
            if(one.theTriangle[j][1] > top_scan_line) top_scan_line = one.theTriangle[j][1];
            if(one.theTriangle[j][1] < bottom_scan_line) bottom_scan_line = one.theTriangle[j][1];
        }
        printf("before top and bottom %f %f\n", top_scan_line, bottom_scan_line);
        if(top_scan_line > y_top_limit) top_scan_line = y_top_limit;
        if(bottom_scan_line < y_bottom_limit) bottom_scan_line = y_bottom_limit;

        printf("after top and bottom %f %f\n", top_scan_line, bottom_scan_line);

        ///triangle equation generation done, hope it's correct -_- ///
        double j = top_scan_line;
        cout<<"first ever ";
        printf("%f\n", top_scan_line);
        double important_left, important_right;
        double left_X = 9999.0, right_X = -9999.0;
        while(j >= bottom_scan_line){
            //cout<<"inside1"<<endl;
            int intersect_count = 0;
            Point intersect_points[3];
            for(int k = 0; k < 3; k++){
                //cout<<"intersections"<<endl;
                //intersection_check(intersect_points[k], one, k, j);
                intersect_points[k] = intersection_check(one, k, j);
                if(j == top_scan_line && intersect_points[k].size()>0){
                    printf("%f %f %f\n", intersect_points[k][0], intersect_points[k][1], j);
                }
                //cout<<"intersection done"<<endl;
                //cout<<k<<endl;
                //cout<<endl;
            }
            for(int k = 0; k < 3; k++){
                //fout1<<k<<endl;
                //cout<<isnan(intersect_points[k][0])<<endl;
                //cout<<"mew"<<endl;
                if(intersect_points[k].size() > 0){
                    //cout<<"inside"<<endl;
                    for(int l = 0; l < intersect_points[k].size(); l++){
                        //printf("%f ", intersect_points[k][l]);
                        //fout1<<intersect_points[k][l]<<" ";
                    }
                //cout<<endl;
                //fout1<<endl;
                }
                /*else {
                    fout1<<"here"<<endl;
                    cout<<"here"<<endl;
                }*/
            }
            //double left_X = 9999.0, right_X = -9999.0;
            double left_Z, right_Z;
            vector<double> xs;
            for(int k = 0; k < 3; k++){
                if(intersect_points[k].size() > 0){
                    //cout<<"the x "<<intersect_points[k][0]<<" the z "<<intersect_points[k][2]<<endl;
                    /*if(k == 0 || k == 1){
                        left_X = intersect_points[k][0];
                        left_Z = intersect_points[k][2];
                    }
                    if(intersect_points[k][0] > right_X){
                        right_X = intersect_points[k][0];
                        right_Z = intersect_points[k][2];
                    }*/
                    xs.push_back(intersect_points[k][0]);
                    xs.push_back(intersect_points[k][2]);
                }
            }
            for(int k = 0; k < xs.size(); k++){
                if(xs.size() == 6){
                    if(j <= Ys[1]){
                        bool found = true;

                        for(int u = 0; u < 3; u++){
                            if(xs[2*u] < Xs[0]){
                                found = false;
                                xs.erase(xs.begin()+2*u, xs.begin()+2*u+2);
                            }
                            else if(xs[2*u] > Xs[2]){
                                found = false;
                                xs.erase(xs.begin()+2*u, xs.begin()+2*u+2);
                            }
                        }

                        if(found) xs.erase(xs.begin()+2, xs.begin()+4);
                    }

                    else if(j <= Ys[0]){
                        //cout<<"okokok"<<endl;
                        xs.erase(xs.begin(), xs.begin()+2);
                    }
                    //cout<<"SIZE "<<xs.size()<<endl;
                }
                /*
                if(xs.size() == 6){
                    //cout<<Left_x<<" "<<j<<" "<<xs[0]<<" "<<xs[2]<<" "<<xs[4]<<endl;
                    //cout<<j<<endl;
                    //xs.erase(xs.begin(), xs.begin()+2);
                    if(xs[4] >= Xs[2]){
                        xs.erase(xs.begin()+4, xs.begin()+6);
                    }
                    else if(xs[2] >= Xs[2]){
                        xs.erase(xs.begin()+2, xs.begin()+4);
                    }
                    else if(xs[0] >= Xs[2]){
                        xs.erase(xs.begin(), xs.begin()+2);
                    }

                    else if(xs[4] <= Xs[0]){
                        cout<<11<<endl;
                        xs.erase(xs.begin()+4, xs.begin()+6);
                    }

                    else if(xs[2] <= Xs[0]){
                        cout<<22<<endl;
                        xs.erase(xs.begin()+2, xs.begin()+4);
                    }

                    else if(xs[0] <= Xs[0]){
                        cout<<33<<endl;
                        xs.erase(xs.begin(), xs.begin()+2);
                    }

                    else if(j <= Ys[1]){
                        xs.erase(xs.begin()+2, xs.begin()+4);
                    }

                    else if(j <= Ys[0]){
                        //cout<<"okokok"<<endl;
                        xs.erase(xs.begin(), xs.begin()+2);
                    }
                    //cout<<"SIZE "<<xs.size()<<endl;
                }*/

                if(xs.size() == 4){
                    if(xs[0] < xs[2]){
                        left_X = xs[0];
                        left_Z = xs[1];
                        right_X = xs[2];
                        right_Z = xs[3];
                    }
                    else{
                        left_X = xs[2];
                        left_Z = xs[3];
                        right_X = xs[0];
                        right_Z = xs[1];
                    }
                }
                important_left = left_X;
                important_right = right_X;
                /*else if(xs.size() == 6){
                    if(xs[0] <= xs[2] && xs[0] <= xs[4]){
                        if(xs[0] < Left_x){
                            if(xs[2]<=xs[4]){
                                left_X = xs[2];
                                left_Z = xs[3];
                                right_X = xs[4];
                                right_Z = xs[5];
                            }
                            else{
                                left_X = xs[4];
                                left_Z = xs[5];
                                right_X = xs[2];
                                right_Z = xs[3];
                            }
                        }
                        else{
                            left_X = xs[0];
                            left_Z = xs[1];


                            if(xs[2] >= xs[4]){
                                if(xs[2]>x_right_limit){
                                    right_X = xs[4];
                                    right_Z = xs[5];
                                }
                                else{
                                    right_X = xs[2];
                                    right_Z = xs[3];
                                }
                            }
                            else{
                                if(xs[4]>x_right_limit){
                                    right_X = xs[2];
                                    right_Z = xs[3];
                                }
                                else{
                                    right_X = xs[4];
                                    right_Z = xs[5];
                                }
                            }
                        }
                    }
                    else if(xs[2] <= xs[0] && xs[2] <= xs[4]){
                        if(xs[2] < Left_x){
                            if(xs[0]<=xs[4]){
                                left_X = xs[0];
                                left_Z = xs[1];
                                right_X = xs[4];
                                right_Z = xs[5];
                            }
                            else{
                                left_X = xs[4];
                                left_Z = xs[5];
                                right_X = xs[0];
                                right_Z = xs[1];
                            }
                        }
                        else{
                            left_X = xs[2];
                            left_Z = xs[3];


                            if(xs[0] >= xs[4]){
                                if(xs[0] > x_right_limit){
                                    right_X = xs[4];
                                    right_Z = xs[5];
                                }
                                else{
                                    right_X = xs[0];
                                    right_Z = xs[1];
                                }
                            }
                            else{
                                if(xs[4] > x_right_limit){
                                    right_X = xs[0];
                                    right_Z = xs[1];
                                }
                                else{
                                    right_X = xs[4];
                                    right_Z = xs[5];
                                }
                            }
                        }
                    }
                    else if(xs[4] <= xs[0] && xs[4] <= xs[2]){
                        if(xs[4] <= Left_x){
                            if(xs[2]<xs[0]){
                                left_X = xs[2];
                                left_Z = xs[3];
                                right_X = xs[0];
                                right_Z = xs[1];
                            }
                            else{
                                left_X = xs[0];
                                left_Z = xs[1];
                                right_X = xs[2];
                                right_Z = xs[3];
                            }
                        }
                        else{
                            left_X = xs[4];
                            left_Z = xs[5];
                        }

                        if(xs[2] >= xs[0]){
                            if(xs[2] > x_right_limit){
                                right_X = xs[0];
                                right_Z = xs[1];
                            }
                            else{
                                right_X = xs[2];
                                right_Z = xs[3];
                            }
                        }
                        else{
                            if(xs[0] > x_right_limit){
                                right_X = xs[2];
                                right_Z = xs[3];
                            }
                            else{
                                right_X = xs[0];
                                right_Z = xs[1];
                            }
                        }
                    }
                }*/
            }
            xs.clear();
            //printf("left %f right %f\n", left_X, right_X);
            if(left_X < x_left_limit){
                left_X = x_left_limit;
            }
            if(right_X > x_right_limit){
                right_X = x_right_limit;
            }

            printf("theLeft %f left %f right %f\n", Left_x, left_X, right_X);
            /*if(right_X >= 0.048740) {
                cout<<"ok"<<endl;
                break;
            }*/
            double k = left_X;
            //cout<<k<<endl;
            while(k <= right_X){
                //cout<<"inside"<<endl;
                //cout<<"ll "<<left_X<<" "<<Left_x<<endl;
                //double l = Left_x>x_left_limit?Left_x:x_left_limit;
                //double lf = l>Left_x?l:left_X;
                //double lf = l;
                //double t = Top_y<y_top_limit?Top_y:y_top_limit;
                //double r = right_X<x_right_limit?right_X:x_right_limit;

                double z = left_Z;
                if(important_left != important_right) z += ((k - important_left)*(left_Z - right_Z)) / (important_left - important_right);

                //cout<<"the z is "<<z<<" "<<endl;
                int buffer_column = ((k - Left_x) / dx);
                //if(buffer_column < 0) buffer_column = -buffer_column;
                int buffer_row = ((Top_y - j) / dy);
                fout1<<buffer_column<<" "<<buffer_row<<" "<<k<<" "<<z<<" "<<dx<<" "<<dy<<endl;
                //if(buffer_row < 0) buffer_row = -buffer_row;
                //fout1<<buffer_column<<" "<<buffer_row<<" "<<k<<" "<<z<<" "<<dx<<" "<<dy<<endl;
                /*if(buffer_column > Screen_Width - 1 || buffer_row > Screen_Height -1) {
                    //cout<<"wtf?"<<endl;
                    break;
                }*/
                //if(buffer_column > 499) buffer_column -= 1;
                //if(buffer_row > 499) buffer_row -= 1;
                //cout<<buffer_column<<" "<<buffer_row<<" "<<z<<endl;

                //if(z_max > z && Z_buffer[buffer_row][buffer_column] < z){
                if(Z_buffer[buffer_row][buffer_column] > z){
                    fout1<<buffer_column<<" "<<buffer_row<<" "<<k<<" "<<z<<" "<<dx<<" "<<dy<<endl;
                    //cout<<"ok"<<endl;

                    Z_buffer[buffer_row][buffer_column] = z;
                    image->set_pixel(buffer_column, buffer_row, (int) one.color[0], (int) one.color[1], (int) one.color[2]);
                    //cout<<(int) one.color[0]<<" "<< (int) one.color[1]<<" "<< (int) one.color[2]<<endl;

                }

                k += dx;
            }
            j -= dy; ///not sure if it's the correct value
        }
        //cout<<"out"<<endl;
        printf("%f %f %f %f\n", top_scan_line, bottom_scan_line, Left_x, right_X);
    }
}

void write_Z_buffer()
{
    ofstream fout("z_buffer.txt");
    for(int i = 0; i < Screen_Height; i++){
        for(int j = 0; j < Screen_Width; j++){
            if(Z_buffer[i][j] < z_max){
                fout <<std::fixed <<std::setprecision(6)<< Z_buffer[i][j]<<"\t";
            }
        }
        fout<<"\n";
    }
    fout.close();
}

int main()
{

    string configs[4];
    read_config_file(configs);

    dx = (x_right_limit - x_left_limit) / Screen_Width;
    dy = (y_top_limit - y_bottom_limit) / Screen_Height;
    //printf("%f %f\n", dx, dy);

    Top_y = y_top_limit - (dy / 2);
    Left_x = x_left_limit + (dx / 2);
    //printf("%f %f\n", Top_y, Left_x);
    fout1<<Top_y<<" "<<Left_x<<endl;
    z_max = z_front_limit>z_rear_limit?z_front_limit:z_rear_limit;
    //cout<<z_max<<endl;
    set_Z_buffer();

    bitmap_image image(Screen_Width, Screen_Height);
    ///not sure if the following nested loop is necessary, i.e., if we
    ///really need to set the background color black, since the image
    ///looks black just after being initialized. will see later
    for(int i = 0; i < Screen_Width; i++){
        for(int j = 0; j < Screen_Height; j++){
            image.set_pixel(i, j, 0, 0, 0);
        }
    }

    ifstream fin("stage3.txt");
    string str;

    int lines = 0;

    while(getline(fin, str)) {
        if(str.compare("")) lines++;
    }
    fin.close();

    string triangles[lines];
    read_stage_3_file(triangles);

    Triangle tri[lines / 3];
    Equations eqn[lines / 3];

    find_triangles(tri, eqn, triangles, lines / 3);

    //print_triangle(tri, lines / 3);
    //cout<<"hello"<<endl;
    apply_Z_buffer_algo(&image, tri, lines/3);

    image.save_image("1.bmp");

    write_Z_buffer();
}
