//
// Created by lc06 on 3/18/2023.
//

#ifndef SPH_CONFIG_H
#define SPH_CONFIG_H

#include <string>

extern std::string path;

extern const double pi;

extern double padding;
extern double radius;

extern int boundaryParticleNum;
extern int fluidParticleNum;
extern int obstacleParticleNum;

extern float x_bias;
extern float y_bias;

extern float gravity;
extern float dt;

extern double h;
extern double grid_size;

extern double c_0;
extern double gamma_;
extern double rho_0;
extern double B;

extern double alpha;

extern double r_0;
extern double D;
extern double n1;
extern double n2;

extern double epsilon;
extern double CFL_V;
extern double CFL_A;

extern double mass;

#endif //SPH_CONFIG_H
