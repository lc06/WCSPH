//
// Created by 25466 on 3/18/2023.
//
#include "particle.h"

Particle::Particle() = default;
Particle::Particle(double p_x, double p_y) : p_x {p_x}, p_y {p_y} {}

void Particle::update_rho(Particle *p) {
    Eigen::Vector2d u_ij, r_ij;
    u_ij << this->v_x - p->v_x, this->v_y - p->v_y;
    r_ij << this->p_x - p->p_x, this->p_y - p->p_y;

    rho += mass * u_ij.dot(r_ij) / r_ij.norm() * kernelFuncDerivative(h, r_ij.norm()) * dt;
}

void Particle::update_pressure() {
    pressure = B * (pow((this->rho / rho_0), gamma_) - 1.);
}

void Particle::calc_gravity_force_acce() {
    a_y += gravity;
}

void Particle::calc_pressure_force_acce(Particle* p) {
    Eigen::Vector2d r_ij;
    r_ij << this->p_x - p->p_x, this->p_y - p->p_y;

    double param = -mass * (this->pressure / (this->rho * this->rho) + p->pressure / (p->rho * p->rho)) * kernelFuncDerivative(h, r_ij.norm());
    a_x += param * (r_ij(0) / r_ij.norm());
    a_y += param * (r_ij(1) / r_ij.norm());
}

void Particle::calc_viscosity_force_acce(Particle* p) {
    Eigen::Vector2d u_ij, r_ij;
    u_ij << this->v_x - p->v_x, this->v_y - p->v_y;
    r_ij << this->p_x - p->p_x, this->p_y - p->p_y;

    double dotVal = u_ij.dot(r_ij);

    if(dotVal < 0) {
        double miu_ij = h * dotVal / (r_ij.norm() * r_ij.norm() + 0.01 * h * h);
        double param = 2. * mass * alpha * miu_ij * (c_0 - miu_ij) / (this->rho + p->rho) * kernelFuncDerivative(h, r_ij.norm());

        a_x += param * (r_ij(0) / r_ij.norm());
        a_y += param * (r_ij(1) / r_ij.norm());
    }
}

void Particle::calc_repulsive_force_acce(Particle* p){
    Eigen::Vector2d r_ij;
    r_ij << this->p_x - p->p_x, this->p_y - p->p_y;

    if(r_ij.norm() < r_0) {
        double param = D * (pow(r_0 / r_ij.norm(), n1) - pow(r_0 / r_ij.norm(), n2));
        a_x += param * r_ij(0) / (r_ij.norm() * r_ij.norm()) / mass;
        a_y += param * r_ij(1) / (r_ij.norm() * r_ij.norm()) / mass;
    }
}

void Particle::calc_velocity() {
    v_x += a_x * dt;
    v_y += a_y * dt;
}

void Particle::calc_position() {
    p_x += v_x * dt;
    p_y += v_y * dt;
}

void Particle::calc_xsph_correction(Particle* p) {
    Eigen::Vector2d u_ji, r_ij;
    u_ji << p->v_x - this->v_x, p->v_y - this->v_y;
    r_ij << this->p_x - p->p_x, this->p_y - p->p_y;

    double param = 2. * epsilon * mass / (this->rho + p->rho) * kernelFunc(h, r_ij.norm());
    p_x += param * u_ji(0) * dt;
    p_y += param * u_ji(1) * dt;

}

void Particle::set_type(int type) {
    this->type = type;
}
int Particle::get_type() {
    return this->type;
}

void Particle::set_id(int id) {
    this->id = id;
}
int Particle::get_id() {
    return this->id;
}

void Particle::set_p_x(double p_x) {
    this->p_x = p_x;
}
double Particle::get_p_x() {
    return this->p_x;
}

void Particle::set_p_y(double p_y) {
    this->p_y = p_y;
}
double Particle::get_p_y() {
    return this->p_y;
}

double Particle::get_v() {
    return sqrt(this->v_x * this->v_x + this->v_y * this->v_y);
}

double Particle::get_a() {
    return sqrt(this->a_x * this->a_x + this->a_y * this->a_y);
}

double Particle::get_rho() {
    return this->rho;
}

void Particle::clear_a_x() {
    this->a_x = 0.0;
}
void Particle::clear_a_y() {
    this->a_y = 0.0;
}

double alphaD(double h) {
    return 15. / (7. * pi * h * h);
}

double kernelFunc(double h, double r){
    double R = r / h;
    double kernel = 0.0;

    if(R < 0) std::cerr << "R cannot be a negative value!" << std::endl;
    if(R < 1) {
        kernel = 2. / 3. - R * R + 1 / 2. * R * R * R;
    } else if(R < 2.) {
        kernel = 1. / 6. * (2. - R) * (2. - R) * (2. - R);
    } else {}

    return alphaD(h) * kernel;
}

double kernelFuncDerivative(double h, double r){
    double R = r / h;
    double derivation = 0.0;

    if(R < 0.) std::cerr << "R cannot be a negative value!" << std::endl;
    if(R < 1.) {
        derivation = -2. * R + 3. / 2. * R * R;
    } else if (R < 2.) {
        derivation = 1. / 2. * (2. - R) * (2. - R);
    } else {}

    return alphaD(h) * derivation;
}