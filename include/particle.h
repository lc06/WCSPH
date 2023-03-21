//
// Created by 25466 on 3/18/2023.
//

#ifndef SPH_PARTICLE_H
#define SPH_PARTICLE_H

#include "config.h"
#include <Eigen/Dense>
#include <iostream>

double alphaD(double h);
double kernelFunc(double h, double r);
double kernelFuncDerivative(double h, double r);

class Particle {
private:
    int type {};
    int id {-1};
    double p_x {0.0}, p_y {0.0};
    double v_x {0.0}, v_y {0.0};
    double a_x {0.0}, a_y {0.0};
//    double f_x {0.0}, f_y {0.0};
    double rho {1000.0};
    double pressure {0.0};

public:
    Particle();
    Particle(double p_x, double p_y);

    void update_rho(Particle* p);
    void update_pressure();
    void calc_gravity_force_acce();
    void calc_pressure_force_acce(Particle* p);
    void calc_viscosity_force_acce(Particle* p);
    void calc_repulsive_force_acce(Particle* p);

//    void calc_acceleration();
    void calc_velocity();
    void calc_position();
    void calc_xsph_correction(Particle* p);

    void set_type(int type);
    int get_type();

    void set_id(int id);
    int get_id();

    void set_p_x(double p_x);
    double get_p_x();

    void set_p_y(double p_y);
    double get_p_y();

    void clear_a_x();
    void clear_a_y();

    double get_v();
    double get_a();
    double get_rho();
};

#endif //SPH_PARTICLE_H
