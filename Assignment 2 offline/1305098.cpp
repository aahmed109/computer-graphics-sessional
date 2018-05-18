#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>
//#include <stack>
#define N 4
using namespace std;
struct Vec
{
	double x,y,z;
};
double aDotxai(double a[], double x[])
{
    double aDotx=a[0]*x[0]+ a[1]*x[1]+ a[2]*x[2];
    return aDotx*a[0];
}

double aDotxaj(double a[], double x[])
{
    double aDotx=a[0]*x[0]+ a[1]*x[1]+ a[2]*x[2];
    return aDotx*a[1];
}

double aDotxak(double a[], double x[])
{
    double aDotx=a[0]*x[0]+ a[1]*x[1]+ a[2]*x[2];
    return aDotx*a[2];
}
double aCrossxi(double a[], double x[])
{
    return a[1]*x[2]-a[2]*x[1];
}

double aCrossxj(double a[], double x[])
{
    return a[2]*x[0]-a[0]*x[2];
}
double aCrossxk(double a[], double x[])
{
    return a[0]*x[1]-a[1]*x[0];
}
void Rod(double x[], double a[],double angle, double ret[])
{
    ret[0]=cos(angle)*x[0]+ (1-cos(angle))*aDotxai(a,x) + sin(angle)*aCrossxi(a,x);
    ret[1]=cos(angle)*x[1]+ (1-cos(angle))*aDotxaj(a,x) + sin(angle)*aCrossxj(a,x);
    ret[2]=cos(angle)*x[2]+ (1-cos(angle))*aDotxak(a,x) + sin(angle)*aCrossxk(a,x);
    if((ret[0]>=0 && ret[0]<0.00001)|| (ret[0]>-0.00001 && ret[0]<=0))
        ret[0]=0;
    if((ret[1]>=0 && ret[1]<0.00001)|| (ret[1]>-0.00001 && ret[1]<=0))
        ret[1]=0;
    if((ret[2]>=0 && ret[2]<0.00001)|| (ret[2]>-0.00001 && ret[2]<=0))
        ret[2]=0;

}
void multiply(double mat1[][4], double mat2[][4], double res[][4])
{
    int i, j, k;
    for (i = 0; i < 4; i++)
    {
        for (j = 0; j < 4; j++)
        {
            res[i][j] = 0;
            for (k = 0; k < 4; k++)
                res[i][j] += mat1[i][k]*mat2[k][j];
        }
    }
}
int main()
{
    int height=-1;
    double stk[20][4][4];
    for(int k=0;k<20;k++)
    {
        for(int i=0;i<4;i++)
        {
            for(int j=0;j<4;j++)
            {
                stk[k][i][j]=0;
            }
        }
    }
    std::ifstream inFile;
    std::ofstream outFile1; std::ofstream outFile2; std::ofstream outFile3;
    inFile.open("scene.txt");
    outFile1.open("stage1.txt");
    outFile2.open("stage2.txt");
    outFile3.open("stage3.txt");
    if(!inFile){
        cerr<<"unable";
        //exit(1);
    }
    int cameraParameters=0;
    double eyeX,eyeY,eyeZ,upX,upY,upZ,lookX,lookY,lookZ,fovY,aspectRatio,near,far;
    inFile>> eyeX>> eyeY>> eyeZ;
    inFile>> lookX>> lookY>> lookZ;
    inFile>> upX>> upY>> upZ;
    inFile>>fovY>>aspectRatio>>near>>far;
    //cout<< fovY<<" "<<aspectRatio<<" "<<near<<" "<< far<<endl;
    double matrix[4][4]={{1,0,0,0},{0,1,0,0},{0,0,1,0},{0,0,0,1}};
            double fovX=0;double t=0; double r=0;
            fovY=fovY*3.14159265358979/180;
            fovX=fovY*aspectRatio;
            t=near*tan(fovY/2);
            r=near*tan(fovX/2);
            double P[4][4]={{near/r,0,0,0},{0,near/t,0,0},{0,0,(-(far+near)/(far-near)),(-(2*far*near)/(far-near))},{0,0,-1,0}};
            double lx=0;double ly=0; double lz=0;double rx=0; double ry=0; double rz=0;double ux=0; double uy=0;double uz=0;
            lx=lookX-eyeX; ly=lookY-eyeY; lz=lookZ-eyeZ;
            double maan=sqrt((lx*lx)+(ly*ly)+(lz*lz)); //normalize korbo
            lx=lx/maan;
            ly=ly/maan;
            lz=lz/maan;

            double l[3]={lx,ly,lz};
            double up[3]={upX,upY,upZ};
            rx=aCrossxi(l,up);
            ry=aCrossxj(l,up);
            rz=aCrossxk(l,up);

            double man=sqrt((rx*rx)+(ry*ry)+(rz*rz)); //normalize korbo
            rx=rx/man;
            ry=ry/man;
            rz=rz/man;
            double rr[3]={rx,ry,rz};
            ux=aCrossxi(rr,l);
            uy=aCrossxj(rr,l);
            uz=aCrossxk(rr,l);
            double V[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
            double T[4][4]={{1,0,0,-eyeX},{0,1,0,-eyeY},{0,0,1,-eyeZ},{0,0,0,1}};
            double R[4][4]={{rx,ry,rz,0},{ux,uy,uz,0},{-lx,-ly,-lz,0},{0,0,0,1}};
            multiply(R,T,V);
    for(std:: string line; getline(inFile,line);)
    {
        if(line.compare("triangle")==0){
            double Point1[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{1,0,0,0}};
            double Point2[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{1,0,0,0}};
            double Point3[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{1,0,0,0}};


            cout<<"printing P"<<endl;
            for(int i=0;i<4;i++)
            {
                for(int j=0;j<4;j++)
                {
                    cout<<P[i][j]<<" ";
                }
                cout<< endl;
            }

            cout<<"triangle"<<endl;
            inFile>>Point1[0][0]>>Point1[1][0]>>Point1[2][0];
            inFile>>Point2[0][0]>>Point2[1][0]>>Point2[2][0];
            inFile>>Point3[0][0]>>Point3[1][0]>>Point3[2][0];
            double result[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
            multiply(matrix,Point1,result);
            for(int i=0;i<3;i++)
            {
                   outFile1<<std::fixed <<std::setprecision(7)<<result[i][0]<<" ";
            }
            outFile1<<endl;
            //stage 2 te pathai 1st point ke


            double result3[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
            double result2[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};

            multiply(V,result,result2);

            if(result2[3][0]!=1)
            {
                for(int i=0;i<4;i++)
                {

                   result2[i][0]=result2[i][0]/result2[3][0];
                }
            }

            for(int i=0;i<3;i++)
            {
                   outFile2<<std::fixed <<std::setprecision(7)<<result2[i][0]<<" ";
            }
            outFile2<<endl;
            multiply(P,result2,result3);
            cout<<"printing result3"<<endl;
            for(int i=0;i<4;i++)
            {
                for(int j=0;j<4;j++)
                {
                    cout<<P[i][j]<<" ";
                }
                cout<< endl;
            }
            if(result3[3][0]!=1)
            {
                for(int i=0;i<4;i++)
                {
                   result3[i][0]=result3[i][0]/result3[3][0];
                }
            }
            for(int i=0;i<3;i++)
            {
                   outFile3<<std::fixed <<std::setprecision(7)<<result3[i][0]<<" ";
            }
            outFile3<<endl;
            //second point ber kori
            multiply(matrix,Point2,result);
            for(int i=0;i<3;i++)
            {
                   outFile1<<result[i][0]<<" ";
            }
            outFile1<<endl;
            multiply(V,result,result2);
            if(result2[3][0]!=1)
            {
                for(int i=0;i<4;i++)
                {
                   result2[i][0]=result2[i][0]/result2[3][0];
                }
            }
            for(int i=0;i<3;i++)
            {
                   outFile2<<std::fixed <<std::setprecision(7)<<result2[i][0]<<" ";
            }
            outFile2<<endl;
            multiply(P,result2,result3);
            if(result3[3][0]!=1)
            {
                for(int i=0;i<4;i++)
                {
                   result3[i][0]=result3[i][0]/result3[3][0];
                }
            }
            for(int i=0;i<3;i++)
            {
                   outFile3<<std::fixed <<std::setprecision(7)<<result3[i][0]<<" ";
            }
            outFile3<<endl;
            //3rd point er ta ber kori
            multiply(matrix,Point3,result);
            for(int i=0;i<3;i++)
            {
                   outFile1<<result[i][0]<<" ";
            }
            outFile1<<endl;
            outFile1<<endl;
            multiply(V,result,result2);
            if(result2[3][0]!=1)
            {
                for(int i=0;i<4;i++)
                {
                   result2[i][0]=result2[i][0]/result2[3][0];
                }
            }
            for(int i=0;i<3;i++)
            {
                   outFile2<<std::fixed <<std::setprecision(7)<<result2[i][0]<<" ";
            }
            outFile2<<endl;
            outFile2<<endl;
            multiply(P,result2,result3);
            if(result3[3][0]!=1)
            {
                for(int i=0;i<4;i++)
                {
                   result3[i][0]=result3[i][0]/result3[3][0];
                }
            }
            for(int i=0;i<3;i++)
            {
                   outFile3<<std::fixed <<std::setprecision(7)<<result3[i][0]<<" ";
            }
            outFile3<<endl;
            outFile3<<endl;
            //outFile1<< result[0][0]<<" "<<result[0][1]<< " "<<result[0][2]<<endl;

        }
        else if(line.compare("translate")==0){
            cout<<"translate"<<endl;
            double tx=0; double ty=0; double tz=0;
            inFile>>tx>>ty>>tz;
            double T[4][4]={{1,0,0,tx},{0,1,0,ty},{0,0,1,tz},{0,0,0,1}};
            double result[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
            multiply(matrix,T,result);
            for(int i=0; i<4;i++)
            {
              for(int j=0;j<4;j++)
              {
                  matrix[i][j]=result[i][j];
              }
            }

        }
        else if(line.compare("scale")==0){
            cout<<"scale"<<endl;
            double sx=0; double sy=0; double sz=0;
            inFile>>sx>>sy>>sz;
            double S[4][4]={{sx,0,0,0},{0,sy,0,0},{0,0,sz,0},{0,0,0,1}};
            double result[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
            multiply(matrix,S,result);
            for(int i=0; i<4;i++)
            {
              for(int j=0;j<4;j++)
              {
                  matrix[i][j]=result[i][j];
              }
            }
        }
        else if(line.compare("rotate")==0){
            cout<<"rotate"<<endl;
            double ax=0; double ay=0; double az=0; double angle=0;
            inFile>>angle>>ax>> ay>> az;
            angle=angle*3.14159265358979/180;
            double maan=sqrt((ax*ax)+(ay*ay)+(az*az)); //normalize korbo
            ax=ax/maan;
            ay=ay/maan;
            az=az/maan;
            double R[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,1}};
            double r[3]={0,0,0};
            double a[3]={ax,ay,az};
            double i[3]={1,0,0}; double j[3]={0,1,0}; double k[3]={0,0,1};
            Rod(i,a,angle,r);
            for(int i=0;i<3;i++)
            {
                R[i][0]=r[i];
            }
            Rod(j,a,angle,r);
            for(int i=0;i<3;i++)
            {
                R[i][1]=r[i];
            }
            Rod(k,a,angle,r);
            for(int i=0;i<3;i++)
            {
                R[i][2]=r[i];
            }
            for(int i=0;i<4;i++)
            {
                for(int j=0;j<4;j++)
                {
                    cout<<R[i][j]<<" ";
                }
                cout<<endl;
            }

            double result[4][4]={{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};
            multiply(matrix,R,result);
            for(int i=0; i<4;i++)
            {
              for(int j=0;j<4;j++)
              {
                  matrix[i][j]=result[i][j];
              }
            }


        }
        else if(line.compare("push")==0){
            cout<<"push"<<endl;
            height++;
            for(int i=0;i<4;i++)
            {
                for(int j=0;j<4;j++)
                {
                    stk[height][i][j]=matrix[i][j];
                }
            }
        }
        else if(line.compare("pop")==0){
            cout<<"pop"<<endl;
            for(int i=0;i<4;i++)
            {
                for(int j=0;j<4;j++)
                {
                    matrix[i][j]=stk[height][i][j];
                }
            }
            height--;
        }
        else if(line.compare("end")==0){
            cout<<"end"<<endl;
            break;
        }
    }
    return 0;
}
