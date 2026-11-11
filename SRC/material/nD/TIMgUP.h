/*
Developed by: Giuseppe Lombardi (giuseppe.lombardi@uniroma1.it), Davide Noè Gorini (davidenoe.gorini@unitn.it)
Date: 2026
Description: this file contains the implementation of the TIM-UP simulating the three-axial force-deformation response at the foundation-superstructure contact of a shallow foundation and the soil interacting with it
Reference: Gorini, D.N., Lombardi, G., Callisto, L. (2025): "A thermodynamic-based macroelement framework for soil-foundation systems with coupled hydro-mechanical behaviour", Computers and Geotechnics, https://doi.org/10.1016/j.compgeo.2025.107681
*/

#ifndef TIMgUP_h
#define TIMgUP_h

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class TIMgUP : public NDMaterial {
public:
    TIMgUP(int tag,
        int N,
        double B,
        double L,
        double tf,
        double Q3act,
        const Vector& Hin_0D,
        const Vector& Hin_0U0,
        const Vector& Qlim_D,
        double br,
        double Q3lim_U0,
        double Q3lim_UC,
        double k_ev,
        double dt_val,
        double Tf_min,
        double Tf_max,
        double a_U,
        double Mq_min,
        double Mq_max,
        double dC0,
        const Vector& r_D = Vector(),
        const Vector& r_U0 = Vector(),
        double gamma = 0.10);

    TIMgUP();
    ~TIMgUP();

    const char* getType(void) const { return "6D"; }
    int getOrder(void) const { return 6; }

    int setTrialStrain(const Vector& v);
    int setTrialStrain(const Vector& v, const Vector& r);
    int setTrialStrainIncr(const Vector& v);
    int setTrialStrainIncr(const Vector& v, const Vector& r);

    const Matrix& getTangent(void);
    const Matrix& getInitialTangent(void);
    const Vector& getStress(void);
    const Vector& getStrain(void);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    int setParameter(const char** argv, int argc, Parameter& param);
    int updateParameter(int responseID, Information& info);

    NDMaterial* getCopy(void);
    NDMaterial* getCopy(const char* code);

    int sendSelf(int commitTag, Channel& theChannel);
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker);

    void Print(OPS_Stream& s, int flag = 0);

private:
    int N_surf;
    double B_foot, L_foot, tf_foot, Q3_act, gamma_ratio;
    double br;

    double k_ev;
    double dt_default;
    double Tf_min, Tf_max, a_U;
    double Mq_min, Mq_max, dC0;

    Vector Qlim_D, Qlim_U0, Qlim_UC;
    Matrix Hin0_D, Hin0_U0, Hin0_UC;
    Vector r_disp_D, r_disp_U0;
    Vector qlim_D, qlim_U0;
    Vector scale;

    Matrix a_D, a_U0, a_UC;
    Matrix H_0D, H_0U0, H_0UC;
    Matrix H_eD, H_eU0, H_eUC;
    Matrix H_pD, H_pU0, H_pUC;

    double commU, trialU;
    double commC, trialC;
    double trialdC, commdC;
    double commMq, trialMq;

    Vector trialq, trialQ, trialqp, trialcs, trialO, trialQint;
    Matrix trialTangent;

    Vector commq, commQ, commqp, commcs, commO, commQint;
    Matrix commTangent;

    Vector commpf_flag, trialpf_flag;
    Vector comm_qpU, trial_qpU;
    Vector trialdQU, trialdQC, commdQU, commdQC;
    Vector trialdq, commdq;

    static Matrix Tangent;
    static Vector Stress;
    static Vector Strain;

    void calibrateModel(void);
    void calibrateHelper(const Vector& Qlim_in, const Matrix& Hin_in, const Vector& qlim_in,
        Matrix& a_out, Matrix& H0_out, Matrix& Hp_out, Matrix& He_out,
        bool is_undrained, const Matrix& a_D_ref);

    Matrix var_UC(const Matrix& x1, const Matrix& x2, const Matrix& x3, double U, double C) const;
    Vector var_UC_vec(const Vector& x1, const Vector& x2, const Vector& x3, double U, double C) const;
    double var_UC_scalar(double x1, double x2, double x3, double U, double C) const;

    double DU_func(const Vector& Q_i, const Vector& cs_a, const Matrix& a_ev, const Matrix& a_D_mat, const Matrix& a_U0_mat, const Matrix& a_UC_mat, double C_val, double br, int n) const;
    double DC_func(const Vector& Q_i, const Vector& cs_a, const Matrix& a_ev, const Matrix& a_U0_mat, const Matrix& a_UC_mat, double U_val, double br, int n) const;

    double yieldFunction(const Vector& Q_i, const Vector& cs_surf, const Matrix& a_ev, double br, int n_surf) const;
    Vector getDQ(const Vector& Q_i, const Vector& cs_surf, const Matrix& a_ev, double br, int n_surf) const;
    Vector getDCs(const Vector& Q_i, const Vector& cs_surf, const Matrix& a_ev, double br, int n_surf) const;
    Vector bisectionQ(const Vector& Q_i, const Vector& v, const Vector& cs_surf, const Matrix& a_ev, double br, int n_surf, double tol = 1e-8) const;
};

#ifdef __cplusplus
extern "C" {
#endif
    void* OPS_TIMgUP(void);
#ifdef __cplusplus
}
#endif

#endif

