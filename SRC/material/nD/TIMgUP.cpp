/*
Developed by: Giuseppe Lombardi (giuseppe.lombardi@uniroma1.it), Davide Noè Gorini (davidenoe.gorini@unitn.it)
Date: 2026
Description: this file contains the implementation of the TIM-UP simulating the three-axial force-deformation response at the foundation-superstructure contact of a shallow foundation and the soil interacting with it
Reference: Gorini, D.N., Lombardi, G., Callisto, L. (2025): "A thermodynamic-based macroelement framework for soil-foundation systems with coupled hydro-mechanical behaviour", Computers and Geotechnics, https://doi.org/10.1016/j.compgeo.2025.107681
*/

#include "TIMgUP.h"
#include <MaterialResponse.h>
#include <OPS_Globals.h>
#include <NDMaterial.h>
#include <Parameter.h>
#include <Information.h>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <elementAPI.h>

Matrix TIMgUP::Tangent(6, 6);
Vector TIMgUP::Stress(6);
Vector TIMgUP::Strain(6);

static inline double dotProduct(const Vector& a, const Vector& b) {
    double sum = 0.0;
    int sz = a.Size();
    for (int i = 0; i < sz; i++) sum += a(i) * b(i);
    return sum;
}

extern "C" void* OPS_TIMgUP(void) {
    NDMaterial* theMaterial = nullptr;
    int numArgs = OPS_GetNumRemainingInputArgs();

    if (numArgs < 24) {
        opserr << "ndMaterial TIMgUP incorrect num args: want tag, N, B, L, tf, Q3act "
            << "HinD1 HinD2 HinD3 HinD4 HinD5 HinD6 HinU01 HinU01 HinU02 HinU03 HinU04 HinU05 HinU06 "
            << "aD1 aD2 QlimD3 aD4 aD5 aD6 br QlimU03 QlimUC3 k_ev dt Tf_min Tf_max a_U Mq_min Mq_max dC0 "
            << "<-r_D rD1 rD2 rD3 rD4 rD5 rD6> <-r_U0 rU01 rU02 rU03 rU04 rU05 rU06> <-gamma gamma>\n";
        return nullptr;
    }

    int tag = 0, N = 10;
    int numData = 1;
    if (OPS_GetIntInput(&numData, &tag) != 0) return nullptr;
    if (OPS_GetIntInput(&numData, &N) != 0) return nullptr;

    double geom[3];
    numData = 3;
    if (OPS_GetDoubleInput(&numData, geom) != 0) return nullptr;
    double B = geom[0], L = geom[1], tf = geom[2];

    double Q3act = 0.0;
    numData = 1;
    if (OPS_GetDoubleInput(&numData, &Q3act) != 0) return nullptr;

    Vector HinD(6), HinU0(6);
    double hdataD[6], hdataU[6];

    numData = 6;
    if (OPS_GetDoubleInput(&numData, hdataD) != 0) return nullptr;
    for (int i = 0; i < 6; i++) HinD(i) = hdataD[i];

    numData = 6;
    if (OPS_GetDoubleInput(&numData, hdataU) != 0) return nullptr;
    for (int i = 0; i < 6; i++) HinU0(i) = hdataU[i];

    Vector QlimD(6);
    double qdata[6];
    numData = 6;
    if (OPS_GetDoubleInput(&numData, qdata) != 0) return nullptr;
    for (int i = 0; i < 6; i++) QlimD(i) = qdata[i];

    double br = 0.0;
    numData = 1;
    if (OPS_GetDoubleInput(&numData, &br) != 0) return nullptr;

    double Q3lim_U0 = 0.0, Q3lim_UC = 0.0;
    numData = 1;
    if (OPS_GetDoubleInput(&numData, &Q3lim_U0) != 0) return nullptr;
    if (OPS_GetDoubleInput(&numData, &Q3lim_UC) != 0) return nullptr;

    double k_ev_val = 1.0e-3, dt_val = 0.01;
    numData = 1;
    if (OPS_GetDoubleInput(&numData, &k_ev_val) != 0) return nullptr;
    if (OPS_GetDoubleInput(&numData, &dt_val) != 0) return nullptr;

    double Tf_min = 0.86, Tf_max = 6.0e6, a_U = 2.997e-5;
    double Mq_min = 0.1, Mq_max = 1.0, dC0 = 0.005;

    double u_params[3], c_params[3];
    numData = 3;
    if (OPS_GetDoubleInput(&numData, u_params) == 0) {
        Tf_min = u_params[0]; Tf_max = u_params[1]; a_U = u_params[2];
    }
    numData = 3;
    if (OPS_GetDoubleInput(&numData, c_params) == 0) {
        Mq_min = c_params[0]; Mq_max = c_params[1]; dC0 = c_params[2];
    }

    Vector r_D(6);
    r_D(0) = 0.01; r_D(1) = 0.01; r_D(2) = 0.1;
    r_D(3) = 0.025; r_D(4) = 0.025; r_D(5) = 0.1;
    Vector r_U0(6);
    r_U0(0) = 0.01; r_U0(1) = 0.01; r_U0(2) = 0.1;
    r_U0(3) = 0.025; r_U0(4) = 0.025; r_U0(5) = 0.1;
    double gamma = 0.10;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetString();
        if (strcmp(opt, "-rD") == 0 || strcmp(opt, "-r_D") == 0) {
            double rd[6]; numData = 6;
            if (OPS_GetDoubleInput(&numData, rd) == 0) for (int i = 0; i < 6; i++) r_D(i) = rd[i];
        }
        else if (strcmp(opt, "-rU0") == 0 || strcmp(opt, "-r_U0") == 0) {
            double ru[6]; numData = 6;
            if (OPS_GetDoubleInput(&numData, ru) == 0) for (int i = 0; i < 6; i++) r_U0(i) = ru[i];
        }
        else if (strcmp(opt, "-gamma") == 0) {
            numData = 1;
            OPS_GetDoubleInput(&numData, &gamma);
        }
    }

    theMaterial = new TIMgUP(tag, N, B, L, tf, Q3act, HinD, HinU0, QlimD, br, Q3lim_U0, Q3lim_UC,
        k_ev_val, dt_val, Tf_min, Tf_max, a_U, Mq_min, Mq_max, dC0,
        r_D, r_U0, gamma);
    return theMaterial;
}

TIMgUP::TIMgUP(int tag, int N, double B, double L, double tf, double Q3act,
    const Vector& Hin_0D_in, const Vector& Hin_0U0_in, const Vector& Qlim_D_in, double br_in, double Q3lim_U0_in, double Q3lim_UC_in,
    double k_ev_in, double dt_in, double tfmin, double tfmax, double au,
    double mqmin, double mqmax, double dc0,
    const Vector& r_D_in, const Vector& r_U0_in, double gamma)
    : NDMaterial(tag, ND_TAG_TIMgUP),
    N_surf(N), B_foot(B), L_foot(L), tf_foot(tf), Q3_act(Q3act), gamma_ratio(gamma),
    k_ev(k_ev_in), dt_default(dt_in), Tf_min(tfmin), Tf_max(tfmax), a_U(au),
    Mq_min(mqmin), Mq_max(mqmax), dC0(dc0),
    Qlim_D(Qlim_D_in), br(br_in), Qlim_U0(6), Qlim_UC(6),
    Hin0_D(6, 6), Hin0_U0(6, 6), Hin0_UC(6, 6),
    r_disp_D(r_D_in), r_disp_U0(r_U0_in), qlim_D(6), qlim_U0(6), scale(6),
    a_D(N, 6), a_U0(N, 6), a_UC(N, 6), H_0D(6, 6), H_0U0(6, 6), H_0UC(6, 6),
    H_eD(6 * N, 6), H_eU0(6 * N, 6), H_eUC(6 * N, 6),
    H_pD(6 * N, 6), H_pU0(6 * N, 6), H_pUC(6 * N, 6),
    trialq(6), trialQ(6), trialTangent(6, 6), trialqp(6 * N), trialcs(6 * N), trialO(6 * N), trialQint(6 * N),
    commq(6), commQ(6), commTangent(6, 6), commqp(6 * N), commcs(6 * N), commO(6 * N), commQint(6 * N),
    commpf_flag(N), trialpf_flag(N), comm_qpU(6), trial_qpU(6),
    commU(1.0), trialU(1.0), commC(0.0), trialC(0.0), commdC(0.0), trialdC(0.0), commMq(0.0), trialMq(0.0),
    trialdQU(6), trialdQC(6), commdQU(6), commdQC(6), commdq(6), trialdq(6)
{
    scale(0) = 1.0; scale(1) = 1.0; scale(2) = 1.0;
    scale(3) = B_foot; scale(4) = L_foot; scale(5) = tf_foot;

    if (r_disp_D.Size() < 6) { r_disp_D.resize(6); for (int i = 0; i < 6; i++) r_disp_D(i) = 0.001; }
    if (r_disp_U0.Size() < 6) { r_disp_U0.resize(6); for (int i = 0; i < 6; i++) r_disp_U0(i) = 0.001; }

    Hin0_D.Zero();
    Hin0_U0.Zero();
    for (int i = 0; i < 6; i++) {
        Hin0_D(i, i) = Hin_0D_in(i);
        Hin0_U0(i, i) = Hin_0U0_in(i);
        if (i >= 3) {
            Qlim_D(i) /= scale(i);
            Hin0_D(i, i) /= (scale(i) * scale(i));
            Hin0_U0(i, i) /= (scale(i) * scale(i));
        }
    }

    Qlim_U0 = Qlim_D;
    Qlim_UC = Qlim_D;
    Qlim_U0(2) = Q3lim_U0_in;
    Qlim_UC(2) = Q3lim_UC_in;

    for (int i = 0; i < 6; i++) {
        if (i != 2) {
            Qlim_U0(i) = Qlim_D(i) * Qlim_U0(2) / Qlim_D(2);
            Qlim_UC(i) = Qlim_D(i) * Qlim_UC(2) / Qlim_D(2);
        }
    }

    Hin0_UC = Hin0_U0;

    qlim_D(0) = r_disp_D(0) * B_foot; qlim_D(1) = r_disp_D(1) * L_foot;
    qlim_D(2) = r_disp_D(2) * B_foot; qlim_D(3) = r_disp_D(3) * B_foot;
    qlim_D(4) = r_disp_D(4) * L_foot; qlim_D(5) = r_disp_D(5) * tf_foot;

    qlim_U0(0) = r_disp_U0(0) * B_foot; qlim_U0(1) = r_disp_U0(1) * L_foot;
    qlim_U0(2) = r_disp_U0(2) * B_foot; qlim_U0(3) = r_disp_U0(3) * B_foot;
    qlim_U0(4) = r_disp_U0(4) * L_foot; qlim_U0(5) = r_disp_U0(5) * tf_foot;

    calibrateModel();

    trialTangent = H_0D; commTangent = H_0D;
    trialq.Zero(); commq.Zero();
    trialQ.Zero(); commQ.Zero();
    trialqp.Zero(); commqp.Zero();
    trialcs.Zero(); commcs.Zero();
    trialO.Zero(); commO.Zero();
    trialQint.Zero(); commQint.Zero();
    trialpf_flag.Zero(); commpf_flag.Zero();
    trial_qpU.Zero(); comm_qpU.Zero();
    trialdQU.Zero(); commdQU.Zero();
    trialdQC.Zero(); commdQC.Zero();
    commdq.Zero();
    trialdq.Zero();

    Matrix a_ev0 = var_UC(a_D, a_U0, a_UC, 1.0, 0.0);
    for (int n = 0; n < N_surf; n++) {
        commO(n * 6 + 2) = a_ev0(n, 2);
        trialO(n * 6 + 2) = a_ev0(n, 2);
    }

    this->commitState();
}

TIMgUP::TIMgUP() : NDMaterial(0, ND_TAG_TIMgUP),
N_surf(0), B_foot(0), L_foot(0), tf_foot(0), Q3_act(0), gamma_ratio(0),
k_ev(1e-3), dt_default(0.01), Tf_min(0.86), Tf_max(6e6), a_U(2.997e-5),
Mq_min(0.1), Mq_max(1.0), dC0(0.005),
Qlim_D(6), br(0), Qlim_U0(6), Qlim_UC(6), Hin0_D(6, 6), Hin0_U0(6, 6), Hin0_UC(6, 6),
r_disp_D(6), r_disp_U0(6), qlim_D(6), qlim_U0(6), scale(6),
a_D(1, 6), a_U0(1, 6), a_UC(1, 6), H_0D(6, 6), H_0U0(6, 6), H_0UC(6, 6),
H_eD(6, 6), H_eU0(6, 6), H_eUC(6, 6), H_pD(6, 6), H_pU0(6, 6), H_pUC(6, 6),
trialq(6), trialQ(6), trialTangent(6, 6), trialqp(6), trialcs(6), trialO(6), trialQint(6),
commq(6), commQ(6), commTangent(6, 6), commqp(6), commcs(6), commO(6), commQint(6),
commpf_flag(1), trialpf_flag(1), comm_qpU(6), trial_qpU(6),
commU(1.0), trialU(1.0), commC(0.0), trialC(0.0), commdC(0.0), trialdC(0.0), commMq(0.0), trialMq(0.0),
trialdQU(6), trialdQC(6), commdQU(6), commdQC(6), commdq(6), trialdq(6) {
    for (int i = 0; i < 6; i++) scale(i) = 1.0;
}

TIMgUP::~TIMgUP() {}

void TIMgUP::calibrateHelper(const Vector& Qlim_in, const Matrix& Hin_in, const Vector& qlim_in,
    Matrix& a_out, Matrix& H0_out, Matrix& Hp_out, Matrix& He_out,
    bool is_undrained, const Matrix& a_D_ref) {

    a_out.Zero();
    for (int j = 0; j < 6; j++) a_out(N_surf - 1, j) = Qlim_in(j);
    a_out(N_surf - 1, 2) = Qlim_in(2) / 2.0;

    for (int j = 0; j < 6; j++) {
        a_out(0, j) = gamma_ratio * a_out(N_surf - 1, j);
    }

    for (int i = 1; i < N_surf - 1; i++) {
        double factor = gamma_ratio + ((double)i / (N_surf - 1.0)) * (1.0 - gamma_ratio);
        for (int j = 0; j < 6; j++) {
            a_out(i, j) = factor * a_out(N_surf - 1, j);
        }
    }

    Vector Q_max(6), alpha(6);
    Matrix q_cal(N_surf, 6), Q_cal(N_surf, 6);

    Q_max(2) = Qlim_in(2);
    alpha(2) = 1.0 - Q_max(2) / (Hin_in(2, 2) * qlim_in(2));

    for (int i = 0; i < N_surf; i++) {
        Q_cal(i, 2) = 2.0 * a_out(i, 2);
        q_cal(i, 2) = (Q_cal(i, 2) / Hin_in(2, 2)) * (1.0 / (1.0 - alpha(2) * Q_cal(i, 2) / Q_max(2)));
    }

    Vector chi3_cal(N_surf);

    if (!is_undrained) {
        double Q3_maxD_cal = Q3_act + a_out(0, 2);
        double Q3_Du_cal = a_out(0, 2);

        for (int n = 0; n < N_surf; n++) {
            if (Q3_maxD_cal > 2.0 * a_out(n, 2)) {
                chi3_cal(n) = 2.0 * a_out(n, 2) - Q3_Du_cal;
            }
            else {
                chi3_cal(n) = Q3_maxD_cal - Q3_Du_cal;
            }
        }
    }
    else {
        double Q3_maxU0_cal = Q3_act + a_D_ref(0, 2);
        double Q3_U0u_cal = a_D_ref(0, 2);
        double Q3_maxD_cal = Q3_act + a_D_ref(0, 2);

        Vector O_3U0(N_surf);
        O_3U0(0) = a_out(0, 2);

        if ((Q3_maxU0_cal - Q3_U0u_cal) > 2.0 * a_out(0, 2)) {
            O_3U0(0) = Q3_maxD_cal - a_D_ref(0, 2);
            chi3_cal(0) = (Q3_maxD_cal - Q3_U0u_cal) - (O_3U0(0) - a_out(0, 2));
        }
        else {
            chi3_cal(0) = Q3_maxU0_cal - Q3_U0u_cal;
        }

        for (int n = 1; n < N_surf; n++) {
            O_3U0(n) = a_out(n, 2);
            if ((Q3_maxU0_cal - Q3_U0u_cal) > 2.0 * a_out(n, 2)) {
                O_3U0(n) = Q3_maxD_cal - a_D_ref(n, 2);
            }
            else {
                chi3_cal(n) = Q3_maxU0_cal - Q3_U0u_cal;
            }

            if ((O_3U0(n) + a_out(n, 2)) < (O_3U0(n - 1) + a_out(n - 1, 2))) {
                O_3U0(n) = (O_3U0(n - 1) + a_out(n - 1, 2)) - a_out(n, 2);
                chi3_cal(n) = (Q3_maxD_cal - Q3_U0u_cal) - (O_3U0(n) - a_out(n, 2));
            }
        }
    }

    for (int j = 0; j < 6; j++) {
        if (j != 2) {
            double term_last = (1.0 - std::pow(chi3_cal(N_surf - 1)/a_out(N_surf - 1, 2) - 1, 2)) / std::pow(br + 1 - chi3_cal(N_surf - 1) / a_out(N_surf - 1, 2), 2);
            Q_max(j) = a_out(N_surf - 1, j) * std::sqrt(std::max(0.0, term_last));
            alpha(j) = 1.0 - Q_max(j) / (Hin_in(j, j) * qlim_in(j));

            for (int i = 0; i < N_surf; i++) {
                double term = (1.0 - std::pow(chi3_cal(i)/a_out(i, 2) - 1, 2)) / std::pow(br + 1 - chi3_cal(i) / a_out(i, 2), 2);
                Q_cal(i, j) = a_out(i, j) * std::sqrt(std::max(0.0, term));
                q_cal(i, j) = (Q_cal(i, j) / Hin_in(j, j)) * (1.0 / (1.0 - alpha(j) * Q_cal(i, j) / Q_max(j)));

            }
        }
    }

    H0_out.Zero();
    for (int j = 0; j < 6; j++) H0_out(j, j) = Q_cal(0, j) / q_cal(0, j);

    Hp_out.Zero();
    He_out.Zero();

    for (int j = 0; j < 6; j++) {
        Hp_out((N_surf - 1) * 6 + j, j) = 0.03 * H0_out(j, j);
    }

    Matrix I(6, 6); I.Zero();
    for (int k = 0; k < 6; k++) I(k, k) = 1.0;

    Matrix H0_inv(6, 6);
    H0_out.Solve(I, H0_inv);

    for (int n = 0; n < N_surf - 1; n++) {
        Matrix Hp_sub(6, 6); Hp_sub.Zero();
        for (int j = 0; j < 6; j++) {
            Hp_out(n * 6 + j, j) = (Q_cal(n + 1, j) - Q_cal(n, j)) / (q_cal(n + 1, j) - q_cal(n, j));
            Hp_sub(j, j) = Hp_out(n * 6 + j, j);
        }

        Matrix Hp_inv(6, 6);
        Hp_sub.Solve(I, Hp_inv);

        Matrix C_eD = Hp_inv - H0_inv;
        if (n > 0) {
            for (int j = 0; j < n; j++) {
                Matrix He_sub(6, 6); He_sub.Zero();
                for (int k = 0; k < 6; k++) He_sub(k, k) = He_out(j * 6 + k, k);
                Matrix He_inv(6, 6); He_sub.Solve(I, He_inv);
                C_eD -= He_inv;
            }
        }
        Matrix He_new(6, 6);
        C_eD.Solve(I, He_new);
        for (int k = 0; k < 6; k++) He_out(n * 6 + k, k) = He_new(k, k);
    }

    Matrix Hp_last(6, 6); Hp_last.Zero();
    for (int k = 0; k < 6; k++) Hp_last(k, k) = Hp_out((N_surf - 1) * 6 + k, k);

    Matrix Hp_last_inv(6, 6);
    Hp_last.Solve(I, Hp_last_inv);

    Matrix C_eD_last = Hp_last_inv - H0_inv;
    
    for (int j = 1; j < N_surf - 1; j++) {
        Matrix He_sub(6, 6); He_sub.Zero();
        for (int k = 0; k < 6; k++) He_sub(k, k) = He_out((j - 1) * 6 + k, k);
        Matrix He_inv(6, 6); He_sub.Solve(I, He_inv);
        C_eD_last -= He_inv;
    }
    Matrix He_last(6, 6);
    C_eD_last.Solve(I, He_last);

    for (int k = 0; k < 6; k++) He_out((N_surf - 1) * 6 + k, k) = He_last(k, k);
    
}

void TIMgUP::calibrateModel(void) {
    calibrateHelper(Qlim_D, Hin0_D, qlim_D, a_D, H_0D, H_pD, H_eD, false, a_D);

    calibrateHelper(Qlim_U0, Hin0_U0, qlim_U0, a_U0, H_0U0, H_pU0, H_eU0, true, a_D);

    Matrix dummy_H0(6, 6), dummy_Hp(N_surf * 6, 6), dummy_He(N_surf * 6, 6);
    calibrateHelper(Qlim_UC, Hin0_UC, qlim_U0, a_UC, dummy_H0, dummy_Hp, dummy_He, true, a_D);
    H_0UC = H_0U0;
    H_pUC = H_pU0;
    H_eUC = H_eU0;
}

Matrix TIMgUP::var_UC(const Matrix& x1, const Matrix& x2, const Matrix& x3, double U, double C) const {
    return x1 * U + (x2 * (1.0 - C) + x3 * C) * (1.0 - U);
}

Vector TIMgUP::var_UC_vec(const Vector& x1, const Vector& x2, const Vector& x3, double U, double C) const {
    return x1 * U + (x2 * (1.0 - C) + x3 * C) * (1.0 - U);
}

double TIMgUP::var_UC_scalar(double x1, double x2, double x3, double U, double C) const {
    return x1 * U + (x2 * (1.0 - C) + x3 * C) * (1.0 - U);
}

double TIMgUP::DU_func(const Vector& Q_i, const Vector& cs_a, const Matrix& a_ev, const Matrix& a_D_mat, const Matrix& a_U0_mat, const Matrix& a_UC_mat, double C_val, double br, int n) const {
    double sum = 0.0;
    Vector da_dU(6);
    for (int i = 0; i < 6; ++i) {
        da_dU(i) = a_D_mat(n, i) - (a_U0_mat(n, i) * (1.0 - C_val) + a_UC_mat(n, i) * C_val);
    }
    double dbr_dU = br - (br * (1.0 - C_val) + br * C_val);
    double fn_1 = std::pow((Q_i(0) - cs_a(n * 6)), 2) / std::pow(a_ev(n, 0), 2);
    double fn_2 = std::pow((Q_i(1) - cs_a(n * 6 + 1)), 2) / std::pow(a_ev(n, 1), 2);
    double fn_31 = std::pow(br + 1.0 - (Q_i(2) - cs_a(n * 6 + 2)) / a_ev(n, 2), 2);
    double fn_4 = std::pow((Q_i(3) - cs_a(n * 6 + 3)), 2) / std::pow(a_ev(n, 3), 2);
    double fn_5 = std::pow((Q_i(4) - cs_a(n * 6 + 4)), 2) / std::pow(a_ev(n, 4), 2);
    double fn_6 = std::pow((Q_i(5) - cs_a(n * 6 + 5)), 2) / std::pow(a_ev(n, 5), 2);
    double DU1 = -2.0 * da_dU(0) * std::pow((Q_i(0) - cs_a(n * 6)), 2) / std::pow(a_ev(n, 0), 3);
    double DU2 = -2.0 * da_dU(1) * std::pow((Q_i(1) - cs_a(n * 6 + 1)), 2) / std::pow(a_ev(n, 1), 3);
    double DU31 = 2.0 * (br + 1.0 - (Q_i(2) - cs_a(n * 6 + 2)) / a_ev(n, 2)) * (dbr_dU + da_dU(2) * (Q_i(2) - cs_a(n * 6 + 2)) / std::pow(a_ev(n, 2), 2));
    double DU32 = 2.0 * ((Q_i(2) - cs_a(n * 6 + 2)) / a_ev(n, 2) - 1.0) * da_dU(2) * ((Q_i(2) - cs_a(n * 6 + 2)) / std::pow(a_ev(n, 2), 2));
    double DU4 = -2.0 * da_dU(3) * std::pow((Q_i(3) - cs_a(n * 6 + 3)), 2) / std::pow(a_ev(n, 3), 3);
    double DU5 = -2.0 * da_dU(4) * std::pow((Q_i(4) - cs_a(n * 6 + 4)), 2) / std::pow(a_ev(n, 4), 3);
    double DU6 = -2.0 * da_dU(5) * std::pow((Q_i(5) - cs_a(n * 6 + 5)), 2) / std::pow(a_ev(n, 5), 3);
    sum = (DU1 + DU2 + DU4 + DU5 + DU6) * fn_31 + (fn_1 + fn_2 + fn_4 + fn_5 + fn_6) * DU31 + DU32;
    return sum;
}

double TIMgUP::DC_func(const Vector& Q_i, const Vector& cs_a, const Matrix& a_ev, const Matrix& a_U0_mat, const Matrix& a_UC_mat, double U_val, double br, int n) const {
    double sum = 0.0;
    Vector da_dC(6);
    for (int i = 0; i < 6; ++i) {
        da_dC(i) = (a_UC_mat(n, i) - a_U0_mat(n, i)) * (1.0 - U_val);
    }
    double dbr_dC = (- br + br) * (1.0 - U_val);
    double fn_1 = std::pow((Q_i(0) - cs_a(n * 6)), 2) / std::pow(a_ev(n,0), 2);
    double fn_2 = std::pow((Q_i(1) - cs_a(n * 6 + 1)), 2) / std::pow(a_ev(n, 1), 2);
    double fn_31 = std::pow(br + 1.0 - (Q_i(2) - cs_a(n * 6 + 2)) / a_ev(n, 2), 2);
    double fn_4 = std::pow((Q_i(3) - cs_a(n * 6 + 3)), 2) / std::pow(a_ev(n, 3), 2);
    double fn_5 = std::pow((Q_i(4) - cs_a(n * 6 + 4)), 2) / std::pow(a_ev(n, 4), 2);
    double fn_6 = std::pow((Q_i(5) - cs_a(n * 6 + 5)), 2) / std::pow(a_ev(n, 5), 2);
    double DC1 = -2.0 * da_dC(0) * std::pow((Q_i(0) - cs_a(n * 6)), 2) / std::pow(a_ev(n, 0), 3);
    double DC2 = -2.0 * da_dC(1) * std::pow((Q_i(1) - cs_a(n * 6 + 1)), 2) / std::pow(a_ev(n, 1), 3);
    double DC31 = 2.0 * (br + 1.0 - (Q_i(2) - cs_a(n * 6 + 2)) / a_ev(n, 2)) * (dbr_dC + da_dC(2) * (Q_i(2) - cs_a(n * 6 + 2)) / std::pow(a_ev(n, 2), 2));
    double DC32 = 2.0 * ((Q_i(2) - cs_a(n * 6 + 2)) / a_ev(n, 2) - 1.0) * da_dC(2) * ((Q_i(2) - cs_a(n * 6 + 2)) / std::pow(a_ev(n, 2), 2));
    double DC4 = -2.0 * da_dC(3) * std::pow((Q_i(3) - cs_a(n * 6 + 3)), 2) / std::pow(a_ev(n, 3), 3);
    double DC5 = -2.0 * da_dC(4) * std::pow((Q_i(4) - cs_a(n * 6 + 4)), 2) / std::pow(a_ev(n, 4), 3);
    double DC6 = -2.0 * da_dC(5) * std::pow((Q_i(5) - cs_a(n * 6 + 5)), 2) / std::pow(a_ev(n, 5), 3);
    sum = (DC1 + DC2 + DC4 + DC5 + DC6) * fn_31 + (fn_1 + fn_2 + fn_4 + fn_5 + fn_6) * DC31 + DC32;
    return sum;
}

double TIMgUP::yieldFunction(const Vector& Q_i, const Vector& cs_surf, const Matrix& a_ev, double br, int n) const {
    double fn = 0.0;
    double fn_1 = std::pow((Q_i(0) - cs_surf(0)), 2) / std::pow(a_ev(n, 0), 2);
    double fn_2 = std::pow((Q_i(1) - cs_surf(1)), 2) / std::pow(a_ev(n, 1), 2);
    double fn_31 = std::pow(br + 1.0 - (Q_i(2) - cs_surf(2)) / a_ev(n, 2), 2);
    double fn_32 = std::pow(((Q_i(2) - cs_surf(2)) / a_ev(n, 2)) - 1.0, 2);
    double fn_4 = std::pow((Q_i(3) - cs_surf(3)), 2) / std::pow(a_ev(n, 3), 2);
    double fn_5 = std::pow((Q_i(4) - cs_surf(4)), 2) / std::pow(a_ev(n, 4), 2);
    double fn_6 = std::pow((Q_i(5) - cs_surf(5)), 2) / std::pow(a_ev(n, 5), 2);
    fn = (fn_1 + fn_2 + fn_4 + fn_5 + fn_6) * fn_31 + fn_32 - 1.0;
    return fn;
}

Vector TIMgUP::getDQ(const Vector& Q_i, const Vector& cs_surf, const Matrix& a_ev, double br, int n) const {
    double fn_1 = std::pow((Q_i(0) - cs_surf(0)), 2) / std::pow(a_ev(n, 0), 2);
    double fn_2 = std::pow((Q_i(1) - cs_surf(1)), 2) / std::pow(a_ev(n, 1), 2);
    double fn_31 = std::pow(br + 1.0 - (Q_i(2) - cs_surf(2)) / a_ev(n, 2), 2);
    double fn_4 = std::pow((Q_i(3) - cs_surf(3)), 2) / std::pow(a_ev(n, 3), 2);
    double fn_5 = std::pow((Q_i(4) - cs_surf(4)), 2) / std::pow(a_ev(n, 4), 2);
    double fn_6 = std::pow((Q_i(5) - cs_surf(5)), 2) / std::pow(a_ev(n, 5), 2);
    double dfn_31 = 2.0 * (-1.0/ a_ev(n, 2)) * (br + 1.0 - (Q_i(2) - cs_surf(2)) / a_ev(n, 2));
    double dfn_32 = 2.0 * (1.0 / a_ev(n, 2)) * (((Q_i(2) - cs_surf(2)) / a_ev(n, 2)) - 1.0);
    Vector DQn(6);
    DQn(0) = 2.0 * fn_31 * ((Q_i(0) - cs_surf(0)) / std::pow(a_ev(n, 0), 2));
    DQn(1) = 2.0 * fn_31 * ((Q_i(1) - cs_surf(1)) / std::pow(a_ev(n, 1), 2));
    DQn(2) = (fn_1 + fn_2 + fn_4 + fn_5 + fn_6) * dfn_31 + dfn_32;
    DQn(3) = 2.0 * fn_31 * ((Q_i(3) - cs_surf(3)) / std::pow(a_ev(n, 3), 2));
    DQn(4) = 2.0 * fn_31 * ((Q_i(4) - cs_surf(4)) / std::pow(a_ev(n, 4), 2));
    DQn(5) = 2.0 * fn_31 * ((Q_i(5) - cs_surf(5)) / std::pow(a_ev(n, 5), 2));
    return DQn;
}

Vector TIMgUP::getDCs(const Vector& Q_i, const Vector& cs_surf, const Matrix& a_ev, double br, int n) const {
    Vector DCsn = getDQ(Q_i, cs_surf, a_ev, br, n);
    DCsn *= -1.0;
    return DCsn;
}

Vector TIMgUP::bisectionQ(const Vector& Q_i, const Vector& v, const Vector& cs_surf, const Matrix& a_ev, double br, int n, double tol) const {
    double lam_min = 0.0;
    double lam_max = 1.0;

    Vector Q_mid = Q_i;
    double f_mid = yieldFunction(Q_mid, cs_surf, a_ev, br, n);

    while (f_mid < -tol || f_mid > tol) {
        double f_min = yieldFunction(Q_i + v * lam_min, cs_surf, a_ev, br, n);
        double f_max = yieldFunction(Q_i + v * lam_max, cs_surf, a_ev, br, n);

        while (f_min * f_max > 0.0) {
            lam_max *= 2.0;
            f_max = yieldFunction(Q_i + v * lam_max, cs_surf, a_ev, br, n);
        }

        double lam_mid = 0.5 * (lam_min + lam_max);
        Q_mid = Q_i + v * lam_mid;
        f_mid = yieldFunction(Q_mid, cs_surf, a_ev, br, n);

        if (f_min * f_mid < 0.0) {
            lam_max = lam_mid;
            f_max = f_mid;
        }
        else {
            lam_min = lam_mid;
            f_min = f_mid;
        }
    }

    return Q_mid;
}

int TIMgUP::setTrialStrain(const Vector& v) {
    for (int i = 0; i < 6; i++) trialq(i) = v(i) * scale(i);
    trialq(2) = -trialq(2);

    Vector dq = trialq - commq;

    if (dq.Norm() < 1e-12) {
        trialQ = commQ;
        trialTangent = commTangent;
        trialcs = commcs;
        trialqp = commqp;
        trial_qpU = comm_qpU;
        trialU = commU;
        trialC = commC;
        trialdC = commdC;
        return 0;
    }

    trialdq = dq;

    double dt_val = ops_Dt;
    if (dt_val <= 0.0) dt_val = dt_default;

    // Drainage coefficient
    double dq_norm = dq.Norm();
    double Tf_val = dq_norm / (k_ev * dt_val);
    double Tf0 = -1.0 / a_U;

    if (Tf_val <= Tf_min) trialU = 1.0;
    else if (Tf_val >= Tf_max) trialU = 0.0;
    else {
        trialU = 1.0 / (a_U * (Tf_val - Tf0));
        if (trialU - commU > 0.0008) trialU = commU + 0.0008;
        else if (trialU - commU < -0.0008) trialU = commU - 0.0008;
    }

    Matrix a_ev = var_UC(a_D, a_U0, a_UC, trialU, commC);
    Vector Q_i = commQ;
    Matrix H_a = commTangent;
    Vector cs_a = commcs;
    Vector sum_qpUa = comm_qpU;
    Vector O_a = commO;
    trialpf_flag = commpf_flag;
    Vector pf_a(N_surf); pf_a.Zero();

    double Tol = 1e-8;

    // 2. ADJUSTMENT OF SURFACES DUE TO U OR C VARIATION
    if (std::abs(trialU - commU) > 0.0 || commdC > 0.0) {
        // Ultimate surface adjustment
        O_a((N_surf - 1) * 6 + 2) = a_ev(N_surf - 1, 2);
        Vector cs_last(6); for (int k = 0; k < 6; k++) cs_last(k) = cs_a((N_surf - 1) * 6 + k);
        if (yieldFunction(Q_i, cs_last, a_ev, br, N_surf - 1) > -Tol) {
            pf_a(N_surf - 1) = 1;
            Vector v_cl(6);
            for (int k = 0; k < 6; k++) v_cl(k) = Q_i(k) - O_a((N_surf - 1) * 6 + k);
            Vector Q_start = Q_i - v_cl * 0.5;
            Q_i = bisectionQ(Q_start, v_cl, cs_last, a_ev, br, N_surf - 1, Tol);
        }

        // First surface adjustment
        if (trialpf_flag(0) > 0) {
            cs_a(2) = O_a(2) - a_ev(0, 2);
            if (cs_a(2) < 0.0) { cs_a(2) = 0.0; O_a(2) = a_ev(0, 2); }
            if ((O_a(2) + a_ev(0, 2)) > 2.0 * a_ev(N_surf - 1, 2)) {
                O_a(2) = 2.0 * a_ev(N_surf - 1, 2) - a_ev(0, 2);
                cs_a(2) = O_a(2) - a_ev(0, 2);
            }
            Vector cs0(6); for (int k = 0; k < 6; k++) cs0(k) = cs_a(k);
            if (yieldFunction(Q_i, cs0, a_ev, br, 0) > -Tol) {
                pf_a(0) = 1;
                Vector v_cl(6); for (int k = 0; k < 6; k++) v_cl(k) = Q_i(k) - O_a(k);
                Vector Q_mid = bisectionQ(Q_i - v_cl * 0.5, v_cl, cs0, a_ev, br, 0, Tol);
                for (int k = 0; k < 6; k++) {
                    cs_a(k) += (Q_i(k) - Q_mid(k));
                    O_a(k) = cs_a(k);
                    trialQint(k) = Q_i(k);
                }
                O_a(2) += a_ev(0, 2);
            }
        }

        // Intermediate surfaces adjustment
        for (int n = 1; n < N_surf - 1; n++) {
            if (trialpf_flag(n) > 0) {
                cs_a(n * 6 + 2) = O_a(n * 6 + 2) - a_ev(n, 2);
                if (cs_a(n * 6 + 2) < 0.0) { cs_a(n * 6 + 2) = 0.0; O_a(n * 6 + 2) = a_ev(n, 2); }
                if ((O_a(n * 6 + 2) + a_ev(n, 2)) < (O_a((n - 1) * 6 + 2) + a_ev(n - 1, 2))) {
                    O_a(n * 6 + 2) = O_a((n - 1) * 6 + 2) + a_ev(n - 1, 2) - a_ev(n, 2);
                    cs_a(n * 6 + 2) = O_a(n * 6 + 2) - a_ev(n, 2);
                }
                if (cs_a(n * 6 + 2) > cs_a((n - 1) * 6 + 2)) {
                    cs_a(n * 6 + 2) = cs_a((n - 1) * 6 + 2);
                    O_a(n * 6 + 2) = cs_a(n * 6 + 2) + a_ev(n, 2);
                }
                if ((O_a(n * 6 + 2) + a_ev(n, 2)) > 2.0 * a_ev(N_surf - 1, 2)) {
                    O_a(n * 6 + 2) = 2.0 * a_ev(N_surf - 1, 2) - a_ev(n, 2);
                    cs_a(n * 6 + 2) = O_a(n * 6 + 2) - a_ev(n, 2);
                }
                Vector csn(6); for (int k = 0; k < 6; k++) csn(k) = cs_a(n * 6 + k);
                if (yieldFunction(Q_i, csn, a_ev, br, n) > -Tol) {
                    pf_a(n) = 1;
                    Vector v_cl(6); for (int k = 0; k < 6; k++) v_cl(k) = Q_i(k) - O_a(n * 6 + k);
                    Vector Q_mid = bisectionQ(Q_i - v_cl * 0.5, v_cl, csn, a_ev, br, n, Tol);
                    for (int k = 0; k < 6; k++) {
                        cs_a(n * 6 + k) += (Q_i(k) - Q_mid(k));
                        O_a(n * 6 + k) = cs_a(n * 6 + k);
                        trialQint(n * 6 + k) = Q_i(k);
                    }
                    O_a(n * 6 + 2) += a_ev(n, 2);
                }
            }
            else {
                O_a(n * 6 + 2) = a_ev(n, 2);
                cs_a(n * 6 + 2) = 0.0;
                if ((O_a(n * 6 + 2) + a_ev(n, 2)) < (O_a((n - 1) * 6 + 2) + a_ev(n - 1, 2))) {
                    O_a(n * 6 + 2) = O_a((n - 1) * 6 + 2) + a_ev(n - 1, 2) - a_ev(n, 2);
                    cs_a(n * 6 + 2) = O_a(n * 6 + 2) - a_ev(n, 2);
                }
                Vector csn(6); for (int k = 0; k < 6; k++) csn(k) = cs_a(n * 6 + k);
                if (yieldFunction(Q_i, csn, a_ev, br, n) > -Tol) {
                    pf_a(n) = 1;
                    Vector v_cl(6); for (int k = 0; k < 6; k++) v_cl(k) = Q_i(k) - O_a(n * 6 + k);
                    Vector Q_mid = bisectionQ(Q_i - v_cl * 0.5, v_cl, csn, a_ev, br, n, Tol);
                    for (int k = 0; k < 6; k++) {
                        cs_a(n * 6 + k) += (Q_i(k) - Q_mid(k));
                        O_a(n * 6 + k) = cs_a(n * 6 + k);
                        trialQint(n * 6 + k) = Q_i(k);
                    }
                    O_a(n * 6 + 2) += a_ev(n, 2);
                }
            }
        }

        H_a = var_UC(H_0D, H_0U0, H_0UC, trialU, commC);
        Matrix I(6, 6); I.Zero(); for (int k = 0; k < 6; k++) I(k, k) = 1.0;
        Matrix C_a(6, 6); H_a.Solve(I, C_a);

        Vector dqU(6); dqU.Zero();
        Vector dqC(6); dqC.Zero();

        Matrix He_ev = var_UC(H_eD, H_eU0, H_eUC, trialU, commC);

        for (int n = 0; n < N_surf; n++) {
            if (pf_a(n) == 1) {
                Vector csn(6); for (int k = 0; k < 6; k++) csn(k) = cs_a(n * 6 + k);
                Matrix He(6, 6); He.Zero(); for (int k = 0; k < 6; k++) He(k, k) = He_ev(n * 6 + k, k);
                Vector DQn = getDQ(Q_i, csn, a_ev, br, n);
                Vector DCsn = getDCs(Q_i, csn, a_ev, br, n);
                double DCsnHeDQn = dotProduct(DCsn, He * DQn);

                double du_val = DU_func(Q_i, cs_a, a_ev, a_D, a_U0, a_UC, commC, br, n);
                double dc_val = DC_func(Q_i, cs_a, a_ev, a_U0, a_UC, trialU, br, n);

                dqU += DQn * (du_val * (trialU - commU) / DCsnHeDQn);
                dqC += DQn * (dc_val * commdC / DCsnHeDQn);

                for (int r = 0; r < 6; r++) {
                    for (int c = 0; c < 6; c++) C_a(r, c) -= (DQn(r) * DQn(c)) / DCsnHeDQn;
                }
            }
        }

        C_a.Solve(I, H_a);
        trialdQU = H_a * dqU;
        trialdQC = H_a * dqC;

        if (pf_a(N_surf - 1) == 0) {
            Q_i += (trialdQU + trialdQC);
            for (int n = 0; n < N_surf - 1; n++) {
                Vector csn(6); for (int k = 0; k < 6; k++) csn(k) = cs_a(n * 6 + k);
                if (yieldFunction(Q_i, csn, a_ev, br, n) > -Tol) {
                    for (int k = 0; k < 6; k++) {
                        cs_a(n * 6 + k) += (trialdQU(k) + trialdQC(k));
                        O_a(n * 6 + k) += (trialdQU(k) + trialdQC(k));
                    }
                }
            }
        }

    }

    if (dotProduct(dq, commdq) < 0) {
        H_a = var_UC(H_0D, H_0U0, H_0UC, trialU, commC);
    }

    Vector Q_a = Q_i + H_a * dq;
    Vector cs_last(6); for (int k = 0; k < 6; k++) cs_last(k) = cs_a((N_surf - 1) * 6 + k);

    double H_a25 = H_a(2, 4);
    double H_a52 = H_a(4, 2);
    double Q5 = Q_a(4);
    double dq5 = dq(4);
    double aaaa = 3;

    if (yieldFunction(Q_a, cs_last, a_ev, br, N_surf - 1) > -Tol) {
        Vector v_pcN = 2.0 * (var_UC(H_0D, H_0U0, H_0UC, trialU, commC) * dq);
        Q_a = bisectionQ(Q_i, v_pcN, cs_last, a_ev, br, N_surf - 1, Tol);
    }

    H_a = var_UC(H_0D, H_0U0, H_0UC, trialU, commC);
    Matrix I(6, 6); I.Zero(); for (int k = 0; k < 6; k++) I(k, k) = 1.0;
    Matrix C_a(6, 6); H_a.Solve(I, C_a);

    trialqp = commqp;
    trial_qpU = comm_qpU;

    Matrix He_ev = var_UC(H_eD, H_eU0, H_eUC, trialU, commC);

    for (int n = 0; n < N_surf; n++) {
        Vector csn(6); for (int k = 0; k < 6; k++) csn(k) = cs_a(n * 6 + k);
        double f_check = yieldFunction(Q_a, csn, a_ev, br, n);
        if (f_check > -Tol) {
            pf_a(n) = 1;
            Vector v_int = Q_a - Q_i;
            Vector Q_out = bisectionQ(Q_i, v_int, csn, a_ev, br, n, Tol);
            for (int k = 0; k < 6; k++) trialQint(n * 6 + k) = Q_out(k);

            Matrix He(6, 6); He.Zero(); for (int k = 0; k < 6; k++) He(k, k) = He_ev(n * 6 + k, k);
            Vector DQn = getDQ(Q_out, csn, a_ev, br, n);
            Vector DCsn = getDCs(Q_out, csn, a_ev, br, n);
            double DCsnHeDQn = dotProduct(DCsn, He * DQn);
            
            for (int r = 0; r < 6; r++) {
                 for (int c = 0; c < 6; c++) C_a(r, c) -= (DQn(r) * DQn(c)) / DCsnHeDQn;
            }

            double lam_n = -dotProduct(DQn, Q_a - Q_i) / DCsnHeDQn;
            for (int k = 0; k < 6; k++) {
                trialqp(n * 6 + k) += lam_n * DQn(k);
                trial_qpU(k) += lam_n * DQn(k) * (1 - trialU) / B_foot;
                
            }
        }
    }

    C_a.Solve(I, trialTangent);
    trialQ = Q_a;

    // 4. DEGRADATION DYNAMIC UPDATE (C)
    trialMq = trial_qpU.Norm();

    double a0 = dC0 * (std::pow(Mq_min, 2) - 2 * Mq_min * Mq_max) / std::pow(Mq_max - Mq_min, 2);
    double a1 = 2 * dC0 * Mq_max / std::pow(Mq_max - Mq_min, 2);
    double a2 = - dC0 / std::pow(Mq_max - Mq_min, 2);

    if (trialMq <= Mq_min) trialC = commC;
    else if (trialMq >= Mq_max) {
        double dC_deg = (a0 + a1 * Mq_max + a2 * std::pow(Mq_max, 2)) * (1.0 - 2 * commC + std::pow(commC, 2));
        trialC = commC + dC_deg;
    }
    else {
        double dC_deg = (a0 + a1 * trialMq + a2 * std::pow(trialMq, 2)) * (1.0 - 2 * commC + std::pow(commC, 2));
        trialC = commC + dC_deg;
        
    }
    if (trialC > 1.0) trialC = 1.0;

    trialdC = trialC - commC;

    for (int n = 0; n < N_surf - 1; n++) {
        if (pf_a(n) == 1) {
            for (int k = 0; k < 6; k++) {
                cs_a(n * 6 + k) += (Q_a(k) - trialQint(n * 6 + k));
            }
        }
        for (int k = 0; k < 6; k++) {
            trialcs(n * 6 + k) = cs_a(n * 6 + k);
            trialO(n * 6 + k) = cs_a(n * 6 + k);
        }
        trialO(n * 6 + 2) += a_ev(n, 2);
    }

    for (int n = 0; n < N_surf; n++) trialpf_flag(n) += pf_a(n);

    return 0;
}

int TIMgUP::setTrialStrain(const Vector& v, const Vector& r) { return setTrialStrain(v); }
int TIMgUP::setTrialStrainIncr(const Vector& v) { return setTrialStrain(commq + v); }
int TIMgUP::setTrialStrainIncr(const Vector& v, const Vector& r) { return setTrialStrainIncr(v); }

int TIMgUP::setParameter(const char** argv, int argc, Parameter& param) {
    if (argc < 1) return -1;

    int argIndex = 0;
    // If argv[0] is "material", skip it to inspect the parameter name in argv[1]
    if (strcmp(argv[0], "material") == 0) {
        argIndex = 1;
        if (argc < 2) return -1; // Ensure a second argument exists
    }

    if (strcmp(argv[argIndex], "k_ev") == 0 ||
        strcmp(argv[argIndex], "k_ev_val") == 0 ||
        strcmp(argv[argIndex], "permeability") == 0 ||
        strcmp(argv[argIndex], "k") == 0) {

        return param.addObject(1, this);
    }
    return -1;
}

int TIMgUP::updateParameter(int responseID, Information& info) {
    if (responseID == 1) {
        this->k_ev = info.theDouble;

        // Output message when parameter value is updated
        opserr << "TIMgUP - Parameter 'k_ev' updated to " << this->k_ev << "\n";

        return 0;
    }
    return -1;
}

const Matrix& TIMgUP::getTangent(void) {
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) Tangent(i, j) = trialTangent(i, j) * (scale(i) * scale(j));
    }
    return Tangent;
}

const Matrix& TIMgUP::getInitialTangent(void) {
    Matrix H_0D_out = H_0D;
    for (int i = 0; i < 6; i++) {
        for (int j = 0; j < 6; j++) H_0D_out(i, j) = H_0D(i, j) * (scale(i) * scale(j));
    }
    return H_0D_out;
}

const Vector& TIMgUP::getStress(void) {
    for (int i = 0; i < 6; i++) Stress(i) = trialQ(i) * scale(i);
    Stress(2) = -Stress(2);
    return Stress;
}

const Vector& TIMgUP::getStrain(void) {
    for (int i = 0; i < 6; i++) Strain(i) = trialq(i) / scale(i);
    Strain(2) = -Strain(2);
    return Strain;
}

int TIMgUP::commitState(void) {
    commq = trialq; commQ = trialQ; commTangent = trialTangent;
    commqp = trialqp; commcs = trialcs; commO = trialO; commQint = trialQint;
    commpf_flag = trialpf_flag; comm_qpU = trial_qpU;
    commU = trialU; commC = trialC; commdC = trialdC; commMq = trialMq;
    commdQU = trialdQU; commdQC = trialdQC; commdq = trialdq;
    return 0;
}

int TIMgUP::revertToLastCommit(void) {
    trialq = commq; trialQ = commQ; trialTangent = commTangent;
    trialqp = commqp; trialcs = commcs; trialO = commO; trialQint = commQint;
    trialpf_flag = commpf_flag; trial_qpU = comm_qpU;
    trialU = commU; trialC = commC; trialdC = commdC; trialMq = commMq;
    trialdQU = commdQU; trialdQC = commdQC; trialdq = commdq;
    return 0;
}

int TIMgUP::revertToStart(void) {
    trialq.Zero(); trialQ.Zero(); trialTangent = H_0D;
    trialqp.resize(6 * N_surf); trialqp.Zero();
    trialcs.resize(6 * N_surf); trialcs.Zero();
    trialO.resize(6 * N_surf); trialO.Zero();
    trialQint.resize(6 * N_surf); trialQint.Zero();
    trialpf_flag.resize(N_surf); trialpf_flag.Zero();
    trial_qpU.Zero(); trialdQU.Zero(); trialdQC.Zero(); commdq.Zero();
    trialU = 1.0; trialC = 0.0; trialdC = 0.0; trialMq = 0.0;

    Matrix a_ev0 = var_UC(a_D, a_U0, a_UC, 1.0, 0.0);
    for (int n = 0; n < N_surf; n++) trialO(n * 6 + 2) = a_ev0(n, 2);

    this->commitState();
    return 0;
}

NDMaterial* TIMgUP::getCopy(void) {
    Vector HinD_vec(6), HinU0_vec(6), QlimD_vec(6);
    for (int i = 0; i < 6; i++) {
        HinD_vec(i) = Hin0_D(i, i) * (scale(i) * scale(i));
        HinU0_vec(i) = Hin0_U0(i, i) * (scale(i) * scale(i));
        QlimD_vec(i) = Qlim_D(i) * scale(i);

    }

    TIMgUP* theCopy = new TIMgUP(this->getTag(), N_surf, B_foot, L_foot, tf_foot, Q3_act,
        HinD_vec, HinU0_vec, QlimD_vec, br, Qlim_U0(2), Qlim_UC(2), k_ev, dt_default, Tf_min, Tf_max, a_U,
        Mq_min, Mq_max, dC0, r_disp_D, r_disp_U0, gamma_ratio);
    return theCopy;
}

NDMaterial* TIMgUP::getCopy(const char* type) {
    if (strcmp(type, "6D") == 0 || strcmp(type, "ThreeDimensional") == 0) return getCopy();
    return nullptr;
}

int TIMgUP::sendSelf(int commitTag, Channel& theChannel) { return 0; }
int TIMgUP::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) { return 0; }

void TIMgUP::Print(OPS_Stream& s, int flag) {
    s << "TIMgUP Material, tag: " << this->getTag() << "\n";
    s << "  B: " << B_foot << ", L: " << L_foot << ", tf: " << tf_foot << "\n";
    s << "  k_ev: " << k_ev << ", dt_default: " << dt_default << "\n";
    s << "  Current U: " << trialU << ", Current C: " << trialC << "\n";
}