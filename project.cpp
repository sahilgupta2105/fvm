#include<iostream>
#include<fstream>
#include<cmath>
#include<algorithm>
#include <iomanip>

using namespace std;

double residuals1[96][48], residuals2[96][48], residuals3[96][48], u[98][50], v[98][50], p[98][50], volume[96][48], area_x[97][48], area_y[96][49]={0.0};
double dt=0.0;
const double beta=20.0;
const double rho=1000;
const double tolerance= 1e-6;
const double c =0.5;

// defining a new variable type for coordinates and normals

struct coordinates{

double x_coord;
double y_coord;

};

coordinates coord[97][49];

coordinates normals_x[97][48];
coordinates normals_y[96][49];

//

void normals();
void volumes();
void areas();
void residuals();
double mass_e_pos(int i, int j);
double mass_e_neg(int i, int j);
double mass_w_pos(int i, int j);
double mass_w_neg(int i, int j);
double mass_n_pos(int i, int j);
double mass_n_neg(int i, int j);
double mass_s_pos(int i, int j);
double mass_s_neg(int i, int j);

main(){



// reading the grid file

 ifstream inFile;

    // Open the file.
  inFile.open("bumpgrid.dat");

  // Read the coordinates from the file.

  for(int j=0;j<49;j++){
    for(int i=0;i<97;i++){
        inFile>>coord[i][j].x_coord;
        inFile>>coord[i][j].y_coord;
    }
  }

   // Close the file.
  inFile.close();

  cout<<coord[0][0].x_coord<<endl;

// calculating volumes, areas, normals

normals();
areas();
volumes();

// iterative loop

for(int it=0;it>0;it++){

        //calculating residuals

        residuals();

        //solving the discretized equations

        for(int j=0;j<50;j++){
            for(int i=0;i<98;i++){
                p[i][j]= p[i][j] - (dt/volume[i][j])*(beta*beta)*(residuals1[i][j]);
                u[i][j] = u[i][j] - (dt/volume[i][j])*(residuals2[i][j]/rho - (residuals1[i][j]*u[i][j])/rho );
                v[i][j] = v[i][j] - (dt/volume[i][j])*(residuals3[i][j]/rho - (residuals1[i][j]*v[i][j])/rho );
            }
        }

        //updating boundary conditions cell values

        for(int k=0;k<50;k++){
            v[0][k]= v[97][k] =0.0;
            p[0][k]= p[97][k] = 101325.0;
            u[0][k]= u[97][k] = 20;
        }
        for(int k=0;k<98;k++){
            v[0][k]= v[49][k] =0.0;
            p[0][k]= p[49][k] = 101325.0;
            u[0][k]= u[49][k] = 20;
        }


        //checking residual

        double maximum=0.0;

        for(int j=0;j<48;j++){
            for(int i=0;i<96;i++){

                if(maximum<(abs(residuals1[i][j])+abs(residuals1[i][j])+abs(residuals1[i][j])))
                    maximum= abs(residuals1[i][j])+abs(residuals1[i][j])+abs(residuals1[i][j]);
            }
        }

        if(maximum<tolerance)
            break;
}


}

// function for calculating normals

void normals(){

// for vertical faces

for(int j=0;j<48;j++){
    for(int i=0;i<97;i++){

        normals_x[i][j].x_coord= coord[i][j+1].y_coord - coord[i][j].y_coord;
        normals_x[i][j].y_coord= coord[i][j].x_coord - coord[i][j+1].x_coord;
    }

}

// for horizontal faces

for(int j=0;j<49;j++){
    for(int i=0;i<96;i++){

        normals_y[i][j].x_coord= coord[i][j].y_coord - coord[i][j+1].y_coord;
        normals_y[i][j].y_coord= coord[i+1][j].x_coord - coord[i][j].x_coord;
    }

}

}

//function for calculating areas

void areas(){

//assuming unit depth

// x-direction areas

for(int j=0;j<48;j++){
    for(int i=0;i<97;i++){

        area_x[i][j]= sqrt( pow(coord[i][j+1].x_coord - coord[i][j].x_coord,2) + pow(coord[i][j+1].y_coord - coord[i][j].y_coord,2)); // distance formula
    }
}

// y-direction areas

for(int j=0;j<49;j++){
    for(int i=0;i<96;i++){

        area_y[i][j]= sqrt( pow( coord[i+1][j].x_coord - coord[i][j].x_coord,2) + pow( coord[i+1][j].y_coord - coord[i][j].y_coord,2)); // distance formula
    }
}

}

//function for calculating volumes

void volumes(){

// assuming unit depth

int dummy_x=0;
int dummy_y=0;

for(int j=0;j<48;j++){
    for(int i=0;i<96;i++){
        volume[i][j]= 0.5*(sqrt( pow((coord[dummy_x][dummy_y].x_coord - coord[dummy_x+1][dummy_y+1].x_coord)*(coord[dummy_x][dummy_y].y_coord - coord[dummy_x+1][dummy_y+1].y_coord),2) + pow((coord[dummy_x+1][dummy_y].x_coord - coord[dummy_x][dummy_y+1].x_coord)*(coord[dummy_x+1][dummy_y].y_coord - coord[dummy_x][dummy_y+1].y_coord),2) ));
        dummy_x= dummy_x+1;
    }
    dummy_y= dummy_y+1;
}

}

// function for calculating the components of residual matrix

void residuals(){

for(int j=0;j<48;j++){
    for(int i=0;i<96;i++){
        residuals1[i][j]= mass_e_pos(i,j) + mass_w_pos(i,j) + mass_n_pos(i,j) + mass_s_pos(i,j)+ mass_e_neg(i,j) + mass_w_neg(i,j) + mass_n_neg(i,j) + mass_s_neg(i,j);
        residuals2[i][j]= mass_e_pos(i,j)*u[i][j] + mass_e_neg(i,j)*u[i+1][j] + mass_w_neg(i,j)*u[i][j] + mass_w_pos(i,j)*u[i-1][j] + mass_n_pos(i,j)*u[i][j] + mass_n_neg(i,j)*u[i][j+1] + mass_s_pos(i,j)*u[i-1][j] + mass_s_neg(i,j)*u[i][j] + 0.5*((p[i+1][j]+p[i][j])*area_x[i+1][j]+ (p[i-1][j]+p[i][j])*area_x[i][j]);
        residuals3[i][j]= mass_e_pos(i,j)*v[i][j] + mass_e_neg(i,j)*v[i+1][j] + mass_w_neg(i,j)*v[i][j] + mass_w_pos(i,j)*v[i-1][j] + mass_n_pos(i,j)*v[i][j] + mass_n_neg(i,j)*v[i][j+1] + mass_s_pos(i,j)*v[i-1][j] + mass_s_neg(i,j)*v[i][j] + 0.5*((p[i][j+1]+p[i][j])*area_y[i+1][j]+ (p[i][j-1]+p[i][j])*area_y[i][j]);
    }
}


}
/*
// function for calculating mass flux due to positive eigen values

double mass_x(int i, int j){

double uk_e, uk_w, uk_n, uk_s, lamda_e, lamda_w, lamda_n, lamda_s, answer=0.0;

uk_e= sqrt( pow(u[i][j]*normals_x[i+1][j].x_coord,2) + pow(v[i][j]*normals_x[i+1][j].y_coord,2) );

uk_w= sqrt( pow(u[i][j]*normals_x[i][j].x_coord,2) + pow(v[i][j]*normals_x[i][j].y_coord,2) );

uk_n= sqrt( pow(u[i][j]*normals_y[i][j+1].x_coord,2) + pow(v[i][j]*normals_y[i][j+1].y_coord,2) );

uk_s= sqrt( pow(u[i][j]*normals_y[i][j].x_coord,2) + pow(v[i][j]*normals_y[i][j].y_coord,2) );

lamda_e= (abs(uk_e)+ sqrt( pow(uk_e,2) + 4.0*pow(beta,2) ) )/2.0;

lamda_w= (abs(uk_w)+ sqrt( pow(uk_w,2) + 4.0*pow(beta,2) ) )/2.0;

lamda_n= (abs(uk_n)+ sqrt( pow(uk_n,2) + 4.0*pow(beta,2) ) )/2.0;

lamda_s= (abs(uk_s)+ sqrt( pow(uk_s,2) + 4.0*pow(beta,2) ) )/2.0;

answer= (rho*max(0,u[i][j]) + (c/2.0/lamda_e)*(p[i][j]-p[i+1][j]))*area_x[i+1][j] + (rho*min(0,u[i+1][j]) + (c/2.0/lamda_e)*(p[i][j]-p[i+1][j]))*area_x[i+1][j] ;

return answer;

}
*/
// function for calculating mass flux at east face

double mass_e_pos(int i, int j){

double uk, lamda, answer=0.0;

uk= sqrt( pow(u[i][j]*normals_x[i+1][j].x_coord,2) + pow(v[i][j]*normals_x[i+1][j].y_coord,2) );

lamda= (abs(uk)+ sqrt( pow(uk,2) + 4.0*pow(beta,2) ) )/2.0;

answer= ((rho*(0>u[i][j]? 0:u[i][j]) + (c/2.0/lamda)*(p[i][j]-p[i+1][j]))*area_x[i+1][j] )*normals_x[i+1][j].x_coord + ((rho*(0>v[i][j]? 0:v[i][j]) + (c/2.0/lamda)*(p[i][j]-p[i+1][j]))*area_x[i+1][j] )*normals_x[i+1][j].y_coord;

return answer;

}

double mass_e_neg(int i, int j){

double uk, lamda, answer=0.0;

uk= sqrt( pow(u[i][j]*normals_x[i+1][j].x_coord,2) + pow(v[i][j]*normals_x[i+1][j].y_coord,2) );

lamda= (abs(uk)+ sqrt( pow(uk,2) + 4.0*pow(beta,2) ) )/2.0;

answer= ( (rho*(0<u[i+1][j]? 0:u[i][j]) + (c/2.0/lamda)*(p[i][j]-p[i+1][j]))*area_x[i+1][j])*normals_x[i+1][j].x_coord + ( (rho*(0<v[i][j]? 0:v[i][j]) + (c/2.0/lamda)*(p[i][j]-p[i+1][j]))*area_x[i+1][j])*normals_x[i+1][j].y_coord;

return answer;

}

// function for calculating mass flux at west face

double mass_w_pos(int i, int j){

double uk, lamda, answer=0.0;

uk= sqrt( pow(u[i][j]*normals_x[i][j].x_coord,2) + pow(v[i][j]*normals_x[i][j].y_coord,2) );

lamda= (abs(uk)+ sqrt( pow(uk,2) + 4.0*pow(beta,2) ) )/2.0;

answer= -((rho*(0>u[i-1][j]? 0:u[i-1][j]) + (c/2.0/lamda)*(p[i-1][j]-p[i][j]))*area_x[i][j] )*normals_x[i][j].x_coord - ((rho*(0>v[i-1][j]? 0:v[i-1][j]) + (c/2.0/lamda)*(p[i-1][j]-p[i][j]))*area_x[i][j] )*normals_x[i][j].y_coord;

return answer;

}

double mass_w_neg(int i, int j){

double uk, lamda, answer=0.0;

uk= sqrt( pow(u[i][j]*normals_x[i][j].x_coord,2) + pow(v[i][j]*normals_x[i][j].y_coord,2) );

lamda= (abs(uk)+ sqrt( pow(uk,2) + 4.0*pow(beta,2) ) )/2.0;

answer= -((rho*(0<u[i][j]? 0:u[i][j]) + (c/2.0/lamda)*(p[i-1][j]-p[i][j]))*area_x[i][j])*normals_x[i][j].x_coord - ((rho*(0<v[i][j]? 0:v[i][j]) + (c/2.0/lamda)*(p[i-1][j]-p[i][j]))*area_x[i][j])*normals_x[i][j].y_coord;

return answer;

}

// function for calculating mass flux at north face

double mass_n_pos(int i, int j){

double uk, lamda, answer=0.0;

uk= sqrt( pow(u[i][j]*normals_y[i][j+1].x_coord,2) + pow(v[i][j]*normals_y[i][j+1].y_coord,2) );

lamda= (abs(uk)+ sqrt( pow(uk,2) + 4.0*pow(beta,2) ) )/2.0;

answer= ((rho*(0>u[i][j]? 0:u[i][j])+ (c/2.0/lamda)*(p[i][j]-p[i][j+1]))*area_y[i][j+1] )*normals_y[i][j+1].x_coord + ((rho*(0>v[i][j]? 0:v[i][j]) + (c/2.0/lamda)*(p[i][j]-p[i][j+1]))*area_y[i][j+1] )*normals_y[i][j].y_coord;

return answer;

}

double mass_n_neg(int i, int j){

double uk, lamda, answer=0.0;

uk= sqrt( pow(u[i][j]*normals_y[i][j+1].x_coord,2) + pow(v[i][j]*normals_y[i][j+1].y_coord,2) );

lamda= (abs(uk)+ sqrt( pow(uk,2) + 4.0*pow(beta,2) ) )/2.0;

answer= ( (rho*(0<u[i][j+1]? 0:u[i][j+1]) + (c/2.0/lamda)*(p[i][j]-p[i][j+1]))*area_y[i][j+1])*normals_y[i][j+1].x_coord + ( (rho*(0<v[i][j+1]? 0:v[i][j+1]) + (c/2.0/lamda)*(p[i][j]-p[i][j+1]))*area_y[i][j+1])*normals_y[i][j].y_coord;

return answer;

}

// function for calculating mass flux at south face

double mass_s_pos(int i, int j){

double uk, lamda, answer=0.0;

uk= sqrt( pow(u[i][j]*normals_y[i][j].x_coord,2) + pow(v[i][j]*normals_y[i][j].y_coord,2) );

lamda= (abs(uk)+ sqrt( pow(uk,2) + 4.0*pow(beta,2) ) )/2.0;

answer= -((rho*(0>u[i][j-1]? 0:u[i][j-1]) + (c/2.0/lamda)*(p[i][j-1]-p[i][j]))*area_y[i][j] )*normals_y[i][j].x_coord + ((rho*(0>v[i][j-1]? 0:v[i][j-1]) + (c/2.0/lamda)*(p[i][j-1]-p[i][j]))*area_y[i][j] )*normals_y[i][j].y_coord;

return answer;

}

double mass_s_neg(int i, int j){

double uk, lamda, answer=0.0;

uk= sqrt( pow(u[i][j]*normals_y[i][j].x_coord,2) + pow(v[i][j]*normals_y[i][j].y_coord,2) );

lamda= (abs(uk)+ sqrt( pow(uk,2) + 4.0*pow(beta,2) ) )/2.0;

answer= -( (rho*(0<u[i][j]? 0:u[i][j]) + (c/2.0/lamda)*(p[i][j-1]-p[i][j]))*area_y[i][j])*normals_y[i][j].x_coord + ((rho*(0<v[i][j]? 0:v[i][j]) + (c/2.0/lamda)*(p[i][j-1]-p[i][j]))*area_y[i][j])*normals_y[i][j].y_coord;

return answer;

}
