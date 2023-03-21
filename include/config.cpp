#include "config.h"
//
// Created by lc06 on 3/18/2023.
//
std::string path = "D:\\__Code__\\Code_Repository\\sph\\";

const double pi {3.14159265358979323846};

double radius = 0.10;
double padding = 0.10;

int boundaryParticleNum = 504;
int fluidParticleNum = 2730;
int obstacleParticleNum = 61;

float x_bias = 0.3;
float y_bias = 0.3;

float gravity = -9.81;
float dt = 0.005;

double h = 0.13;
double grid_size = 2 * h;

double c_0 = 200.0;
double gamma_ = 7.0;
double rho_0 = 1000.0;
double B = c_0 * c_0 * rho_0 / gamma_;

double alpha = 1.0;

double r_0 = 0.1;
double D = 400;
double n1 = 12;
double n2 = 4;

double epsilon = 0.3;

double CFL_V = 0.25;
double CFL_A = 0.05;

double mass = 4. / 3. * radius * radius * radius * rho_0;