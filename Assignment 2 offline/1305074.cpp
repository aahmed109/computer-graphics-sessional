#include<stdio.h>
#include<fstream>
#include<string>
#include<iostream>
#include<sstream>
#include<stack>
#include<vector>
#include<cmath>
#include<iomanip>
#define PI 3.14159265

using namespace std;

ofstream fout1("stage1.txt");
ofstream fout2("stage2.txt");
ofstream fout3("stage3.txt");

typedef vector<vector<double>> mat;
typedef vector<double> row;
stack<mat> theStack;
row stackEntry;
mat multiply_view_projection(mat s, row rown)
{
    row thePoint;
    mat mult;
    /*cout<<"input1"<<endl;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++) {
            printf("%f ", s[i][j]);
            //fout1 <<std::fixed <<std::setprecision(7)<< mult[i][j]<<" ";
        }
        cout<<endl;
    }

    cout<<"input2"<<endl;
    for(int i = 0; i < 4; i++) printf("%f ", rown[i]);
    cout<<endl;*/
    for(int i = 0; i < 4; i++){
        double sss = 0.0;
        for(int k = 0; k < 4; k++){
            sss += s[i][k] * rown[k];
        }
        thePoint.push_back(sss);
        mult.push_back(thePoint);
        thePoint.clear();
    }
    /*cout<<"multiplied"<<endl;
    for(int i = 0; i < 4; i++){
        printf("%f ", mult[i][0]);
    }
    cout<<endl;*/
    return mult;
}
void multiply(mat matrix)
{
    mat s = theStack.top();
    mat mult;
    row nn;
    /*for(int i = 0; i < 4; i++){
        cout<<"here"<<endl;
        for(int j = 0; j < 4; j++) mult[i][j] = 0.0;
    }*/
    /*cout<<"stack top"<<endl;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            printf("%f ", s[i][j]);
        }
        cout<<endl;
    }*/

    for(int i = 0; i < 4; i++){
        //cout<<"here"<<endl;
        for(int j = 0; j < 4; j++){
            //cout<<"here"<<endl;
            double sss = 0.0;
            for(int k = 0; k < 4; k++){
                //cout<<"here"<<endl;
                sss += s[i][k] * matrix[k][j];
            }
            nn.push_back(sss);
        }
        mult.push_back(nn);
        nn.clear();
    }
    /*cout<<"mutipled"<<endl;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            printf("%f ", mult[i][j]);
        }
        cout<<endl;
    }*/
    theStack.push(mult);
}

mat only_mult_is_real(mat R, mat T)
{
    mat V;
    row nn;

    for(int i = 0; i < 4; i++){
        //cout<<"here"<<endl;
        for(int j = 0; j < 4; j++){
            //cout<<"here"<<endl;
            double sss = 0.0;
            for(int k = 0; k < 4; k++){
                //cout<<"here"<<endl;
                sss += R[i][k] * T[k][j];
            }
            nn.push_back(sss);
        }
        V.push_back(nn);
        nn.clear();
    }
    /*cout<<"mutipled"<<endl;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            printf("%f ", V[i][j]);
        }
        cout<<endl;
    }*/

    cout<<"V"<<endl;
    for(int i1 = 0; i1 < 4; i1++){
        for(int j = 0; j < 4; j++){
            printf("%f ", V[i1][j]);
        }
        cout<<endl;
    }

    return V;
}

void multiply_for_triangle(mat V, mat P, row rown)
{
    mat s = theStack.top();
    /*cout<<"stack top"<<endl;
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++) printf("%f ", s[i][j]);
        cout<<endl;
    }

    for(int i = 0; i < 4; i++){
        printf("%f ", rown[i]);
    }
    cout<<endl;*/
    row thePoint;
    mat mult;
    row rowns, proj;
    for(int i = 0; i < 4; i++){
        double sss = 0.0;
        for(int k = 0; k < 4; k++){
            sss += s[i][k] * rown[k];
        }
        thePoint.push_back(sss);
        mult.push_back(thePoint);
        thePoint.clear();
    }
    //cout<<"Testing"<<endl;
    for(int i = 0; i < 4; i++){
        //printf("%f\n", mult[i][3]);
        for(int j = 0; j <1; j++){
            //if(mult[i][3] > 1)  mult[i][j] /= mult[i][3];
            mult[i][j] /= mult[3][0];
        }
    }
    //cout<<"tested"<<endl;
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < 1; j++) {
            //printf("%f ", mult[i][j]);
            fout1 <<std::fixed <<std::setprecision(7)<< mult[i][j]<<" ";
        }
    }
    fout1<<"\n";
    cout<<endl;
    for(int i = 0; i < 4; i++){
        rowns.push_back(mult[i][0]);
    }
    /*cout<<"V"<<endl;
    for(int i1 = 0; i1 < 4; i1++){
        for(int j = 0; j < 4; j++){
            printf("%f ", V[i1][j]);
        }
        cout<<endl;
    }
    cout<<"rowns"<<endl;
    for(int i = 0; i < 3; i++){
        printf("%f ", rowns[i]);
    }
    cout<<endl;*/
    mat sq = multiply_view_projection(V, rowns);

    for(int i = 0; i < 4; i++){
        for(int j = 0; j <1; j++){
            sq[i][j] /= sq[3][0];
        }
    }
    //cout<<"\nView"<<endl;
    for(int ut = 0; ut < 3; ut++){
        //printf("%f ", sq[ut][0]);
        fout2 <<std::fixed <<std::setprecision(7)<< sq[ut][0]<<" ";
    }
    fout2<<"\n";

    for(int i = 0; i < 4; i++){
        proj.push_back(sq[i][0]);
    }

    /*cout<<"\nPROJ"<<endl;
    for(int i = 0; i < 4; i++){
        printf("%f ", proj[i]);
    }
    cout<<endl;*/
    mat pp = multiply_view_projection(P, proj);

    for(int i = 0; i < 4; i++){
        for(int j = 0; j <1; j++){
            pp[i][j] /= pp[3][0];
        }
    }
    for(int ut = 0; ut < 3; ut++){
        //printf("%f ", pp[ut][0]);
        fout3 <<std::fixed <<std::setprecision(7)<< pp[ut][0]<<" ";
    }
    fout3<<"\n";
}

row normalize(row a)
{
    double mag = 0.0;
    for(int i = 0; i < 3; i++){
        mag += a[i] * a[i];
    }
    double magr = sqrt(mag);
    for(int i = 0; i < 3; i++){
        a[i] /= magr;
    }
    return a;
}
row Rodrigues(row x, row a, double angle)
{
    row s;
    //printf("%f\n",angle);
    double ang = angle * PI / 180.0;

    double xx = (cos(ang) * x[0]) + ((1 - cos(ang)) * (a[0] * x[0] + a[1] * x[1] + a[2] * x[2]) * a[0]) + (sin(ang) * (a[1] * x[2] - a[2] * x[1]));
    double yy = (cos(ang) * x[1]) + ((1 - cos(ang)) * (a[0] * x[0] + a[1] * x[1] + a[2] * x[2]) * a[1]) + (sin(ang) * (a[2] * x[0] - a[0] * x[2]));
    double zz = (cos(ang) * x[2]) + ((1 - cos(ang)) * (a[0] * x[0] + a[1] * x[1] + a[2] * x[2]) * a[2]) + (sin(ang) * (a[0] * x[1] - a[1] * x[0]));
    s.push_back(xx);
    s.push_back(yy);
    s.push_back(zz);

    return s;
}

int main()
{
    mat identity;
    row id;

    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            if(i == j) id.push_back(1);
            else id.push_back(0);
        }
        identity.push_back(id);
        id.clear();
    }
    theStack.push(identity);

    int lc = 0;
    int lines = 0;
    ifstream fin("scene.txt");
    //ofstream fout("stage1.txt");
    string str;
    while(getline(fin, str)){
        lines++;
    }
    fin.close();

    string commands[lines];
    fin.open("scene.txt");
    while(getline(fin, str)) commands[lc++] = str;
    fin.close();
    mat T, R, V, P;
    for(int i = 0; i < lines; i++){

        //printf("here?\n");
        if(i == 0){
            row eye, look, up;

            stringstream ss(commands[i]);
            double s;
            for(int il = 0; il < 3; il++){
                ss>>s;
                eye.push_back(s);
            }
            stringstream ss1(commands[i + 1]);

            for(int il = 0; il < 3; il++){
                ss1>>s;
                look.push_back(s);
            }
            stringstream ss2(commands[i + 2]);

            for(int il = 0; il < 3; il++){
                ss2>>s;
                up.push_back(s);
            }
            row l, r, u;
            row aa;
            l.push_back(look[0] - eye[0]);
            l.push_back(look[1] - eye[1]);
            l.push_back(look[2] - eye[2]);
            l = normalize(l);
            r.push_back(l[1] * up[2] - l[2] * up[1]);
            r.push_back(l[2] * up[0] - l[0] * up[2]);
            r.push_back(l[0] * up[1] - l[1] * up[0]);
            r = normalize(r);
            u.push_back(r[1] * l[2] - r[2] * l[1]);
            u.push_back(r[2] * l[0] - r[0] * l[2]);
            u.push_back(r[0] * l[1] - r[1] * l[0]);

            //view transform for lc = 0-2, projection transform for lc = 3
            //lc++;
            /*cout<<"l"<<endl;
            for(int i1 = 0; i1 < 3; i1++) printf("%f ", l[i1]);
            cout<<endl;
            cout<<"r"<<endl;
            for(int i1 = 0; i1 < 3; i1++) printf("%f ", r[i1]);
            cout<<endl;
            cout<<"u"<<endl;
            for(int i1 = 0; i1 < 3; i1++) printf("%f ", u[i1]);
            cout<<endl;*/
            for(int i1 = 0; i1 < 4; i1++){
                for(int j = 0; j < 4; j++){
                    if(i1 == j) aa.push_back(1);
                    else if(i1 != 3 && j == 3) aa.push_back(-eye[i1]);
                    else aa.push_back(0);
                }
                T.push_back(aa);
                aa.clear();
            }

            /*cout<<"T"<<endl;
            for(int i1 = 0; i1 < 4; i1++){
                for(int j = 0; j < 4; j++){
                    printf("%f ", T[i1][j]);
                }
                cout<<endl;
            }*/
            for(int i1 = 0; i1 < 4; i1++){
                for(int j = 0; j < 4; j++){
                    if(i1 == 0 && j == 0) aa.push_back(r[0]);
                    else if(i1 == 0 && j == 1) aa.push_back(r[1]);
                    else if(i1 == 0 && j == 2) aa.push_back(r[2]);
                    else if(i1 == 0 && j == 3) aa.push_back(0);
                    else if(i1 == 1 && j == 0) aa.push_back(u[0]);
                    else if(i1 == 1 && j == 1) aa.push_back(u[1]);
                    else if(i1 == 1 && j == 2) aa.push_back(u[2]);
                    else if(i1 == 1 && j == 3) aa.push_back(0);
                    else if(i1 == 2 && j == 0) aa.push_back(-l[0]);
                    else if(i1 == 2 && j == 1) aa.push_back(-l[1]);
                    else if(i1 == 2 && j == 2) aa.push_back(-l[2]);
                    else if(i1 == 2 && j == 3) aa.push_back(0);
                    else if(i1 == 3 && j == 0) aa.push_back(0);
                    else if(i1 == 3 && j == 1) aa.push_back(0);
                    else if(i1 == 3 && j == 2) aa.push_back(0);
                    else if(i1 == 3 && j == 3) aa.push_back(1);
                }
                R.push_back(aa);
                aa.clear();
            }
            /*cout<<"R"<<endl;
            for(int i1 = 0; i1 < 4; i1++){
                for(int j = 0; j < 4; j++){
                    printf("%f ", R[i1][j]);
                }
                cout<<endl;
            }*/
            V = only_mult_is_real(R, T);
            i += 2;
        }

        else if(i == 3){

            row rown, aa;
            stringstream ss(commands[i]);
            for(int l = 0; l < 4; l++){
                double s;
                ss>>s;
                rown.push_back(s);
            }

            double fovX = rown[0] * rown[1];
            double t = rown[2] * tan(((rown[0]/2) * PI)/180.0);
            double r = rown[2] * tan(((fovX/2) * PI)/180.0);
            for(int i1 = 0; i1 < 4; i1++){
                for(int j = 0; j < 4; j++){
                    if(i1 == 0 && j == 0)    aa.push_back(rown[2]/r);
                    if(i1 == 0 && j == 1)    aa.push_back(0);
                    if(i1 == 0 && j == 2)    aa.push_back(0);
                    if(i1 == 0 && j == 3)    aa.push_back(0);
                    if(i1 == 1 && j == 0)    aa.push_back(0);
                    if(i1 == 1 && j == 1)    aa.push_back(rown[2]/t);
                    if(i1 == 1 && j == 2)    aa.push_back(0);
                    if(i1 == 1 && j == 3)    aa.push_back(0);
                    if(i1 == 2 && j == 0)    aa.push_back(0);
                    if(i1 == 2 && j == 1)    aa.push_back(0);
                    if(i1 == 2 && j == 2)    aa.push_back(-((rown[3] + rown[2]) / (rown[3] - rown[2])));
                    if(i1 == 2 && j == 3)    aa.push_back(-((2 * rown[3] * rown[2]) / (rown[3] - rown[2])));
                    if(i1 == 3 && j == 0)    aa.push_back(0);
                    if(i1 == 3 && j == 1)    aa.push_back(0);
                    if(i1 == 3 && j == 2)    aa.push_back(-1);
                    if(i1 == 3 && j == 3)    aa.push_back(0);
                }
                P.push_back(aa);
                aa.clear();
            }
            /*cout<<"P"<<endl;
            for(int i1 = 0; i1 < 4; i1++){
                for(int j = 0; j < 4; j++){
                    printf("%f ", P[i1][j]);
                }
                cout<<endl;
            }*/
        }

        else if(commands[i].compare("triangle") == 0){
            cout<<"triangle"<<endl;
            /*
            cout<<"P"<<endl;
            for(int i1 = 0; i1 < 4; i1++){
                for(int j = 0; j < 4; j++){
                    printf("%f ", P[i1][j]);
                }
                cout<<endl;
            }
            */
            mat points;
            row rown;
            cout<<"multiplication result"<<endl;
            for(int l = 0; l < 3; l++){
                stringstream ss(commands[i+l+1]);
                double s;
                for(int j = 0; j < 3; j++){
                    ss>>s;
                    rown.push_back(s);
                }
                rown.push_back(1);
                multiply_for_triangle(V, P, rown);
                points.push_back(rown);
                rown.clear();
            }
            fout1<<"\n";
            fout2<<"\n";
            fout3<<"\n";
            cout<<endl;
            /*cout<<"TRIANGLE"<<endl;
            for(int iii = 0; iii < 3; iii++){
                for(int j = 0; j < 3; j++){
                    printf("%f ", points[iii][j]);
                }
                cout<<endl;
            }*/
        }

        else if(commands[i].compare("translate") == 0){
            cout<<"translate"<<endl;
            row rown;
            mat trans;
            row id;
            //double trans[3];
            stringstream ss(commands[i+1]);
            for(int j = 0; j < 3; j++){
                double sss;
                ss>>sss;
                rown.push_back(sss);
            }
            //for(int j = 0; j < 3; j++)  printf("%f ", rown[j]);
            cout<<endl;
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    if(i == j) id.push_back(1);
                    else if(j == 3) id.push_back(rown[i]);
                    else id.push_back(0);
                }
                trans.push_back(id);
                id.clear();
            }

            multiply(trans);
            /*for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    printf("%f ", trans[i][j]);
                }
                cout<<endl;
            }*/
        }

        else if(commands[i].compare("scale") == 0){
            cout<<"scale"<<endl;
            row rown;
            mat scal;
            row id;
            //double trans[3];
            stringstream ss(commands[i+1]);
            for(int j = 0; j < 3; j++){
                double sss;
                ss>>sss;
                rown.push_back(sss);
            }
            //for(int j = 0; j < 3; j++)  printf("%f ", rown[j]);
            cout<<endl;
            for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    if(i == j && i < 3) id.push_back(rown[i]);
                    else if(i == 3 && j == 3) id.push_back(1);
                    else id.push_back(0);
                }
                scal.push_back(id);
                id.clear();
            }
            multiply(scal);

            /*for(int i = 0; i < 4; i++){
                for(int j = 0; j < 4; j++){
                    printf("%f ", trans[i][j]);
                }
                cout<<endl;
            }*/        }

        else if(commands[i].compare("rotate") == 0){
            cout<<"rotate"<<endl;
            row rown;
            row a,b;
            mat rot;
            row id;
            row ii, jj, kk;
            row c1, c2, c3;
            ii.push_back(1);
            ii.push_back(0);
            ii.push_back(0);
            jj.push_back(0);
            jj.push_back(1);
            jj.push_back(0);
            kk.push_back(0);
            kk.push_back(0);
            kk.push_back(1);
            //double trans[3];
            stringstream ss(commands[i+1]);
            for(int j = 0; j < 4; j++){
                double sss;
                ss>>sss;
                rown.push_back(sss);
                if(j > 0) b.push_back(sss);
            }
            //for(int j = 0; j < 4; j++)  printf("%f ", rown[j]);
            cout<<endl;
            a = normalize(b);
            c1 = Rodrigues(ii, a, rown[0]);
            c2 = Rodrigues(jj, a, rown[0]);
            c3 = Rodrigues(kk, a, rown[0]);
            /*for(int iiiii = 0; iiiii < 3; iiiii++) printf("%f ", c1[iiiii]);
            cout<<endl;
            for(int iiiii = 0; iiiii < 3; iiiii++) printf("%f ", c2[iiiii]);
            cout<<endl;
            for(int iiiii = 0; iiiii < 3; iiiii++) printf("%f ", c3[iiiii]);
            cout<<endl;*/
            for(int j = 0; j < 4; j++){
                for(int k = 0; k < 4; k++){
                    if(j == 3 && k == 3) id.push_back(1);
                    else if(j == 3 || k == 3) id.push_back(0);
                    else if(j == 0 && k == 0) id.push_back(c1[0]);
                    else if(j == 0 && k == 1) id.push_back(c2[0]);
                    else if(j == 0 && k == 2) id.push_back(c3[0]);
                    else if(j == 1 && k == 0) id.push_back(c1[1]);
                    else if(j == 1 && k == 1) id.push_back(c2[1]);
                    else if(j == 1 && k == 2) id.push_back(c3[1]);
                    else if(j == 2 && k == 0) id.push_back(c1[2]);
                    else if(j == 2 && k == 1) id.push_back(c2[2]);
                    else if(j == 2 && k == 2) id.push_back(c3[2]);
                }
                rot.push_back(id);
                id.clear();
            }
            /*cout<<"Rotating Matrix"<<endl;
            for(int m = 0; m < 4; m++){
                for(int nn = 0; nn < 4; nn++){
                    printf("%f ", rot[m][nn]);
                }
                cout<<endl;
            }*/
            multiply(rot);
        }

        else if(commands[i].compare("push") == 0){
            cout<<"push"<<endl;
            stackEntry.push_back(theStack.size());
            cout<<theStack.size()<<endl;
        }

        else if(commands[i].compare("pop") == 0){
            cout<<"pop"<<endl;
            int q = theStack.size();
            //cout<<q<<" "<<stackEntry.back()<<endl;
            if(q == stackEntry.back()){
                for(int ii = 0; ii < q - 1; ii++){
                    theStack.pop();
                }
            }
            else{
                for(int ii = 0; ii < q - stackEntry.back(); ii++){
                    theStack.pop();
                }
            }
            //cout<<theStack.size()<<endl;
            //if(theStack.size() == 0) theStack.push(identity);
        }

        else if(commands[i].compare("end") == 0) break;
    }
    //for(int i = 0; i < lines; i++) cout<<commands[i]<<endl;
}
