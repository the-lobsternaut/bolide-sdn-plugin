/**
 * Bolide Plugin — 6DOF Integration Tests
 */

#include "bolide/sixdof_core.h"
#include <iostream>
#include <cassert>
#include <cmath>

using namespace sixdof;

void testMeteorTumble() {
    State s;
    s.pos = {0, 0, 80e3};
    s.vel = {-20000, 0, -10000};
    s.quat = qidentity();
    s.omega = {5.0, 3.0, -2.0};
    s.mass = 10000;

    InertiaTensor I = inertiaDiag(500, 1200, 900); // Diagonal for stability
    double dt = 0.0001, t = 0; // Smaller step for 22 km/s entry

    auto aeroFn = [](const State& st, double) -> ForcesTorques {
        ForcesTorques ft;
        double alt = st.pos[2];
        double rho = 1.225 * std::exp(-alt / 8500.0);
        double V = v3norm(st.vel);
        double qbar = 0.5 * rho * V * V;
        double coeff = qbar * 0.1;
        ft.torque_body = {-coeff*0.3, coeff*0.5, -coeff*0.2};
        return ft;
    };

    for (int i = 0; i < 10000; i++) { s = rk4Step(s, I, dt, t, aeroFn); t += dt; }

    assert(std::abs(qnorm(s.quat) - 1.0) < 1e-6);

    std::cout << "  Meteor tumble ✓ (omega=[" << s.omega[0] << "," << s.omega[1] << "," << s.omega[2] << "])\n";
}

void testFragmentSpin() {
    State s;
    s.quat = qidentity();
    s.omega = {0, 0, 10.0};
    s.mass = 5;
    InertiaTensor I = inertiaDiag(0.5, 0.8, 0.3);
    auto coastFn = [](const State&, double) -> ForcesTorques { return {}; };
    double dt = 0.01, t = 0;

    double T0 = 0.5*(I[0]*s.omega[0]*s.omega[0]+I[1]*s.omega[1]*s.omega[1]+I[2]*s.omega[2]*s.omega[2]);
    for (int i = 0; i < 10000; i++) { s = rk4Step(s, I, dt, t, coastFn); t += dt; }
    double T1 = 0.5*(I[0]*s.omega[0]*s.omega[0]+I[1]*s.omega[1]*s.omega[1]+I[2]*s.omega[2]*s.omega[2]);

    assert(std::abs(T1 - T0) / T0 < 1e-4);
    std::cout << "  Fragment spin ✓ (KE drift=" << std::abs(T1-T0)/T0*100 << "%)\n";
}

int main() {
    std::cout << "=== bolide 6DOF tests ===\n";
    testMeteorTumble();
    testFragmentSpin();
    std::cout << "All bolide 6DOF tests passed.\n";
    return 0;
}
