#include <fstream>//file input and output
#include <string> //using string
#include <vector> //using vector
#include <iostream> //basic io in cpp, cout and cin
#include <sstream> //use string stream
#include <math.h> //using sqrt
//#include <Eigen> //for linear algebra
#include <eigen3/Eigen/Dense>
#include <time.h> //measure run time of this program
#include <chrono> 

using namespace std;
using namespace Eigen;

int readDimension(const string& s, vector <double>& v){
    istringstream is(s);
    double n;
    while (is >> n) {
        v.push_back(n);
    }
    return v.size();
}

void load_file(const char* file_name, vector <double> &v, int& rows, int&cols){
    ifstream file;
    string line;

    file.open(file_name);
    if (file.is_open()){
        int i = 0;
        getline(file, line);

        cols = readDimension(line, v);

        for (i = 1; i < 100; i++){
            if (getline(file, line) == 0) break;
            readDimension(line, v);
        }
        rows = i;
        if (rows > 100) cout << "Input data is not valid";
        file.close();
        file.clear();
    }
    else{
        cout << "Can not open file";
    }


    for (int i = 0; i < rows; i++){
        for (int j = 0; j < cols; j++){
            cout << v[i*cols + j];
        }
        // cout << v << endl;
        cout << endl;
    }
}


double clc_miu(double p, double miu_sc = 0.6){
    return miu_sc + 0.00134*log(p);
}
double clc_B(double p, double c = 3e-5, double p_sc = 14.7){
    return pow((1.0+ c*(p-p_sc)), -1);
}

double clc_rho(double p, double c = 3e-5, double rho_sc = 42.0){
    return rho_sc*(1.0+c*(p-14.7));
}

double clc_trans_x(double kx, double kx1, double h, double h1, double delta_x, double delta_x1, double delta_y, double delta_y1, double miu, double B, double beta_c = 1.127e-3){
    double Ax = delta_y*h;
    double Ax1 = delta_y*h1;
    return 2.0*beta_c/(miu*B)/(delta_x/(Ax*kx) + (delta_x1/(Ax1*kx1)));
}

double clc_trans_y(double ky, double ky1, double h, double h1, double delta_x, double delta_x1, double delta_y, double delta_y1, double miu, double B, double beta_c = 1.127e-3){
    double Ay = delta_x*h;
    double Ay1 = delta_x1*h1;
    return 2.0*beta_c/(miu*B)/(delta_y/(Ay*ky) + (delta_y1/(Ay1*ky1)));
}

double clc_trans_x_(double kx, double kx1, double h, double h1, double delta_x, double delta_x1, double delta_y, double delta_y1, double rho, double miu, double B, double beta_c = 1.127e-3){
    double Ax = delta_y*h;
    double Ax1 = delta_y*h1;
    return 2.0*beta_c*rho/(144.0*miu*B)/(delta_x/(Ax*kx) + (delta_x1/(Ax1*kx1)));
}

double clc_trans_y_(double ky, double ky1, double h, double h1, double delta_x, double delta_x1, double delta_y, double delta_y1, double rho, double miu, double B, double beta_c = 1.127e-3){
    double Ay = delta_x*h;
    double Ay1 = delta_x1*h1;
    return 2.0*beta_c*rho/(144.0*miu*B)/(delta_y/(Ay*ky) + (delta_y1/(Ay1*ky1)));
}


double clc_gamma(double delta_x, double delta_y, double h, double phi, double delta_t = 30.0, double ct = 3.1e-5, double alpha_c = 5.615){
    return delta_x*delta_y*h*phi*ct/(alpha_c*delta_t);
}

double sum(int a, double b){

    return a+2*b;
}







int main() {
    clock_t start, end;
    
    // vector <double> v;
    // int rows = 0;
    // int cols = 0;
    // load_file("cpp_kx.txt", v, rows, cols);

    ifstream inFile;
    double kx[15][29];
    double ky[15][29];
    double dx[15][29];
    double dy[15][29];
    double h[15][29];
    double top[15][29];
    double poro[15][29];
    int row;
    int col;
    inFile.open("cpp_kx.txt");
    {
        if (!inFile){
            cout << "File is not found" << endl;
        }
        else{
            for (row =0; row < 15; row++){
                for (col =0; col < 29; col++){
                    inFile >> kx[row][col];
                }
            }

        }
    }
    inFile.close();


    inFile.open("cpp_ky.txt");
    {
        if (!inFile){
            cout << "File is not found" << endl;
        }
        else{
            for (row =0; row < 15; row++){
                for (col =0; col < 29; col++){
                    inFile >> ky[row][col];
                }
            }

        }
    }
    inFile.close();


    inFile.open("cpp_dx.txt");
    {
        if (!inFile){
            cout << "File is not found" << endl;
        }
        else{
            for (row =0; row < 15; row++){
                for (col =0; col < 29; col++){
                    inFile >> dx[row][col];
                }
            }

        }
    }
    inFile.close();


    inFile.open("cpp_dy.txt");
    {
        if (!inFile){
            cout << "File is not found" << endl;
        }
        else{
            for (row =0; row < 15; row++){
                for (col =0; col < 29; col++){
                    inFile >> dy[row][col];
                }
            }

        }
    }
    inFile.close();


    inFile.open("cpp_h.txt");
    {
        if (!inFile){
            cout << "File is not found" << endl;
        }
        else{
            for (row =0; row < 15; row++){
                for (col =0; col < 29; col++){
                    inFile >> h[row][col];
                }
            }

        }
    }
    inFile.close();


    inFile.open("cpp_top.txt");
    {
        if (!inFile){
            cout << "File is not found" << endl;
        }
        else{
            for (row =0; row < 15; row++){
                for (col =0; col < 29; col++){
                    inFile >> top[row][col];
                }
            }

        }
    }
    inFile.close();

    inFile.open("cpp_poro.txt");
    {
        if (!inFile){
            cout << "File is not found" << endl;
        }
        else{
            for (row =0; row < 15; row++){
                for (col =0; col < 29; col++){
                    inFile >> poro[row][col];
                }
            }

        }
    }
    inFile.close();


    int cnt = 0;
    int order[15][29];
    double Pre[15][29];
    double N[15][29];
    double S[15][29];
    double W[15][29];
    double E[15][29];
    double N_[15][29];
    double S_[15][29];
    double W_[15][29];
    double E_[15][29];
    double C[15][29];
    double Q[15][29];
    double gamma[15][29];
    for (int i =0; i< 15; i++){
        for (int j = 0; j < 29; j++){
            if (kx[i][j] > 0){
                cnt += 1;
                order[i][j] = cnt;
                Pre[i][j] = 5000.0;
            }
        }
    }
    start = clock();
    for (int time =0; time < 360*2+30; time += 30){
        for (int i =0; i < 15; i++){
            for (int j =0; j < 29; j++){
                if (order[i][j] > 0){
                    if (order[i-1][j] > 0){
                        N[i][j] = clc_trans_x(kx[i-1][j], kx[i][j], h[i-1][j], h[i][j], dx[i-1][j], dx[i][j], dy[i-1][j], dy[i][j], clc_miu((Pre[i-1][j] + Pre[i][j])/2), clc_B((Pre[i-1][j] + Pre[i][j])/2));
                        N_[i][j] = clc_trans_x_(kx[i-1][j], kx[i][j], h[i-1][j], h[i][j], dx[i-1][j], dx[i][j], dy[i-1][j], dy[i][j],clc_rho((Pre[i-1][j] + Pre[i][j])/2), clc_miu((Pre[i-1][j] + Pre[i][j])/2), clc_B((Pre[i-1][j] + Pre[i][j])/2));
                    }
                    if (order[i+1][j] > 0){
                        S[i][j] = clc_trans_x(kx[i+1][j], kx[i][j], h[i+1][j], h[i][j], dx[i+1][j], dx[i][j], dy[i+1][j], dy[i][j], clc_miu((Pre[i+1][j] + Pre[i][j])/2), clc_B((Pre[i+1][j] + Pre[i][j])/2));
                        S_[i][j] = clc_trans_x_(kx[i+1][j], kx[i][j], h[i+1][j], h[i][j], dx[i+1][j], dx[i][j], dy[i+1][j], dy[i][j],clc_rho((Pre[i+1][j] + Pre[i][j])/2), clc_miu((Pre[i+1][j] + Pre[i][j])/2), clc_B((Pre[i+1][j] + Pre[i][j])/2));
                    }
                    if (order[i][j-1] > 0){
                        W[i][j] = clc_trans_y(ky[i][j-1], ky[i][j], h[i][j-1], h[i][j], dx[i][j-1], dx[i][j], dy[i][j-1], dy[i][j], clc_miu((Pre[i][j-1] + Pre[i][j])/2), clc_B((Pre[i][j-1] + Pre[i][j])/2));
                        W_[i][j] = clc_trans_y_(ky[i][j-1], ky[i][j], h[i][j-1], h[i][j], dx[i][j-1], dx[i][j], dy[i][j-1], dy[i][j],clc_rho((Pre[i][j-1] + Pre[i][j])/2), clc_miu((Pre[i][j-1] + Pre[i][j])/2), clc_B((Pre[i][j-1] + Pre[i][j])/2));
                    }
                    if (order[i][j+1] > 0){
                        E[i][j] = clc_trans_y(ky[i][j+1], ky[i][j], h[i][j+1], h[i][j], dx[i][j+1], dx[i][j], dy[i][j+1], dy[i][j], clc_miu((Pre[i][j+1] + Pre[i][j])/2), clc_B((Pre[i][j+1] + Pre[i][j])/2));
                        E_[i][j] = clc_trans_y_(ky[i][j+1], ky[i][j], h[i][j+1], h[i][j], dx[i][j+1], dx[i][j], dy[i][j+1], dy[i][j],clc_rho((Pre[i][j+1] + Pre[i][j])/2), clc_miu((Pre[i][j+1] + Pre[i][j])/2), clc_B((Pre[i][j+1] + Pre[i][j])/2));
                    }
                    gamma[i][j] = clc_gamma(dx[i][j], dy[i][j], h[i][j], poro[i][j]);
                    C[i][j] = -(N[i][j] + S[i][j] + W[i][j] + E[i][j] +gamma[i][j]);
                    Q[i][j] = N_[i][j]*(top[i-1][j] - top[i][j]) + S_[i][j]*(top[i+1][j] - top[i][j])+W_[i][j]*(top[i][j-1] - top[i][j]) + E_[i][j]*(top[i][j+1] - top[i][j]) - gamma[i][j]*Pre[i][j];
                }
            }
        }
        double re1 = 0.14*sqrt(pow((ky[4][6]/kx[4][6]),0.5)*pow(dx[4][6],2) + pow((kx[4][6]/ky[4][6]),0.5)*pow(dy[4][6],2))/(pow((ky[4][6]/kx[4][6]),0.25)+pow((kx[4][6]/ky[4][6]),0.25));
        double omega1 = 2.0*M_PI*1.127e-3*sqrt(kx[4][6]*ky[4][6])*h[4][6]/(clc_miu(Pre[4][6])*clc_B(Pre[4][6])*log(re1/0.25));
        // cout<<omega1<<endl;
        double re2 = 0.14*sqrt(pow((ky[9][6]/kx[9][6]),0.5)*pow(dx[9][6],2) + pow((kx[9][6]/ky[9][6]),0.5)*pow(dy[9][6],2))/(pow((ky[9][6]/kx[9][6]),0.25)+pow((kx[9][6]/ky[9][6]),0.25));
        double omega2 = 2.0*M_PI*1.127e-3*sqrt(kx[9][6]*ky[9][6])*h[9][6]/(clc_miu(Pre[9][6])*clc_B(Pre[9][6])*log(re1/0.25));

        double re3 = 0.14*sqrt(pow((ky[4][20]/kx[4][20]),0.5)*pow(dx[4][20],2) + pow((kx[4][20]/ky[4][20]),0.5)*pow(dy[4][20],2))/(pow((ky[4][20]/kx[4][20]),0.25)+pow((kx[4][20]/ky[4][20]),0.25));
        double omega3 = 2.0*M_PI*1.127e-3*sqrt(kx[4][20]*ky[4][20])*h[4][20]/(clc_miu(Pre[4][20])*clc_B(Pre[4][20])*log(re1/0.25));

        double re4 = 0.14*sqrt(pow((ky[8][10]/kx[8][10]),0.5)*pow(dx[8][10],2) + pow((kx[8][10]/ky[8][10]),0.5)*pow(dy[8][10],2))/(pow((ky[8][10]/kx[8][10]),0.25)+pow((kx[8][10]/ky[8][10]),0.25));
        double omega4 = 2.0*M_PI*1.127e-3*sqrt(kx[8][10]*ky[8][10])*h[8][10]/(clc_miu(Pre[8][10])*clc_B(Pre[8][10])*log(re1/0.25));

        double re5 = 0.14*sqrt(pow((ky[12][13]/kx[12][13]),0.5)*pow(dx[12][13],2) + pow((kx[12][13]/ky[12][13]),0.5)*pow(dy[12][13],2))/(pow((ky[12][13]/kx[12][13]),0.25)+pow((kx[12][13]/ky[12][13]),0.25));
        double omega5 = 2.0*M_PI*1.127e-3*sqrt(kx[12][13]*ky[12][13])*h[12][13]/(clc_miu(Pre[12][13])*clc_B(Pre[12][13])*log(re1/0.25));

        double re6 = 0.14*sqrt(pow((ky[6][14]/kx[6][14]),0.5)*pow(dx[6][14],2) + pow((kx[6][14]/ky[6][14]),0.5)*pow(dy[6][14],2))/(pow((ky[6][14]/kx[6][14]),0.25)+pow((kx[6][14]/ky[6][14]),0.25));
        double omega6 = 2.0*M_PI*1.127e-3*sqrt(kx[6][14]*ky[6][14])*h[6][14]/(clc_miu(Pre[6][14])*clc_B(Pre[6][14])*log(re1/0.25));
        C[4][6] -= omega1;
        C[9][6] -= omega2;
        C[4][20] -= omega3;
        C[8][10] -= omega4;
        C[12][13] -= omega5;
        C[6][14] -= omega6;
        Q[4][6] -= omega1*1000.0;
        Q[9][6] -= omega2*1000.0;
        Q[4][20] -= omega3*1000.0;
        Q[8][10] -= omega4*1000.0;
        Q[12][13] -= omega5*1000.0;
        Q[6][14] -= omega6*1000.0;


        // cout << clc_trans_x(100.0, 100.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, clc_miu(5000.0), clc_B(5000.0)) << endl;
        // cout << Q[4][6] <<endl;

        MatrixXf A(cnt,cnt);
        VectorXf b(cnt);
        for (int i=0; i < 15; i++){
            for (int j =0; j < 29; j++){
                if (order[i][j] > 0){
                    A(order[i][j] - 1, order[i][j] - 1) = C[i][j];
                    if (order[i-1][j] > 0){
                        A(order[i][j] - 1, order[i-1][j] - 1) = N[i][j];
                    }
                    if (order[i+1][j] > 0){
                        A(order[i][j] - 1,order[i+1][j] - 1) = S[i][j];
                    }
                    if (order[i][j-1] > 0){
                        A(order[i][j] - 1, order[i][j-1] - 1) = W[i][j];
                    }
                    if (order[i][j+1] > 0){
                        A(order[i][j] - 1,order[i][j+1] - 1) = E[i][j];
                    }
                    b(order[i][j] - 1) = Q[i][j];
                }
            }
        }
        VectorXf x = A.lu().solve(b);
        // for (int i = 0; i < cnt; i++){
        //     cout << x(i) <<endl;
        //
        // }
        cout << Pre[4][6] << endl;
        // cout << x(order[4][6]-1) << endl;
        // cout <<M_PI << endl;

        for (int i = 0; i < 15; i++){
            for (int j =0; j < 29; j++ ){
                if (order[i][j] > 0){
                    Pre[i][j] = x[order[i][j] - 1];
                }
            }
        }

    

    }
    end = clock();
    double cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    cout << "Time take by this program is :" << cpu_time_used <<setprecision(9);
    cout << "sec" <<endl;
}
