/*
Developed by: Davide Noè Gorini (davidenoe.gorini@unitn.it), Giuseppe Lombardi (giuseppe.lombardi@uniroma1.it)				  
Date: 2025
Description: this file contains the implementation of the TIM-UP simulating the three-axial force-deformation response at the foundation-superstructure contact of a shallow foundation and the soil interacting with it
Degrees of Freedom covered by the TIM-UP: 1 - horizontal translation, 3 - vertical translation, 5 - rotation around axis 2
Reference: Gorini, D.N., Lombardi, G., Callisto, L. (2025): "A thermodynamic-based macroelement framework for soil-foundation systems with coupled hydro-mechanical behaviour", Computers and Geotechnics, https://doi.org/10.1016/j.compgeo.2025.107681
*/

#ifndef TIMUP3D_h
#define TIMUP3D_h

#include <NDMaterial.h>

#include <T2Vector.h>
#include <Matrix.h>
#include <Vector.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>

#include <Information.h>
#include <Parameter.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <string.h>

#include <elementAPI.h>

class TIMUP3D : public NDMaterial
{
public:

	TIMUP3D(int tag, int nYield,
		double H11el_D, double H33el_D, double H55el_D, double mHyperEl,
		double Q3min_D, double Q3max_D, 
		double a1_DN, double a5_DN, double b_DN, double r_DN,
		double H11el_U0, double H33el_U0, double H55el_U0, 
		double Q3min_U0, double Q3max_U0,
		double H11el_UC, double H33el_UC, double H55el_UC, 
		double Q3min_UC, double Q3max_UC,
		double perm, double Lfound, double beta,
		double Tf_min, double Tf_max, double a_U, double Mq_min, double Mq_max, double dC0,
		double tol, double deltaT, 
		double Q3_maxD_cal, double Q3_maxU0_cal, double Q3_maxUC_cal, double Q3_Du_cal, double Q3_U0u_cal, double Q3_UCu_cal,
		int niter,
		double b_11D = 70.1, double b_33D = 1058.2, double b_55D = 50.4, double b_11U0 = 21.4, double b_33U0 = 224.3, double b_55U0 = 17.9);

	TIMUP3D();
	virtual ~TIMUP3D();

	const char* getClassType(void) const { return "TIMUP3D"; };
	const char* getType(void) const { return "ThreeDimensional"; };
	int setTrialStrain(const Vector& strain);
	int setTrialRelVel(const Vector& RelVelNodes);
	int setTrialStrain(const Vector& v, const Vector& r);
	int setTrialStrainIncr(const Vector& v);
	int setTrialStrainIncr(const Vector& v, const Vector& r);
	int getOrder() const;

	// Calculates current tangent stiffness
	const Matrix& getTangent(void);
	const Matrix& getInitialTangent(void);
	
	// Calculates the corresponding stress increment (rate), for a given strain increment
	const Vector& getStress(void);
	const Vector& getStrain(void);
	const Vector& getCommittedStress(void);
	const Vector& getCommittedStrain(void);


	int commitState(void);
	int revertToLastCommit(void);
	int revertToStart(void);

	NDMaterial* getCopy(void);

	NDMaterial* getCopy(const char* code);


	int sendSelf(int commitTag, Channel& theChannel);
	int recvSelf(int commitTag, Channel& theChannel,
		FEM_ObjectBroker& theBroker);

	Response* setResponse(const char** argv, int argc, OPS_Stream& s);
	int getResponse(int responseID, Information& matInformation);
	void Print(OPS_Stream& s, int flag = 0);

	int setParameter(const char** argv, int argc, Parameter& param);
	int updateParameter(int responseID, Information& eleInformation);

protected:

private:

	int ComputeResponse();

	int nYield;
	int niter;
	double H11el;
	double H33el;
	double H55el;

	double H11el_D;
	double H33el_D;
	double H55el_D;
	double H11el_U0;
	double H33el_U0;
	double H55el_U0;
	double H11el_UC;
	double H33el_UC;
	double H55el_UC;
	double Q3min_D;
	double Q3max_D;
	double Q3max_U0;
	double Q3min_U0;
	double Q3max_UC;
	double Q3min_UC;

	double mHyperEl;

	double a1_DN;
	double a5_DN;
	double b_DN;
	double r_DN;
	double a1_U0N;
	double a5_U0N;
	double b_U0N;
	double r_U0N;
	double a1_UCN;
	double a5_UCN;
	double b_UCN;
	double r_UCN;

	double fq3_DN;
	double fq3_U0N;
	double fq3_UCN;

	double Q30_DN;
	double Q30_U0N;
	double Q30_UCN;

	double perm;
	double Lfound;
	double beta;

	Matrix a1;
	Matrix a3;
	Matrix a5;
	Matrix b;
	Matrix r;
	Matrix c10;
	Matrix c30;
	Matrix c50;
	Matrix Q3min;
	Matrix Q3max;
	Matrix fq3;
	Matrix Q30;

	Vector m;

	Matrix H_0D;
	Matrix H_0U0;
	Matrix H_0UC;

	Matrix H_eD;
	Matrix H_eU0;
	Matrix H_eUC;


	double Tf_min;
	double Tf_max;
	double a_U;
	double Tf0;
	double Mq_min;
	double Mq_max;
	double dC0;

	double Q3_maxD_cal;
	double Q3_maxU0_cal;
	double Q3_maxUC_cal;
	double Q3_Du_cal;
	double Q3_U0u_cal;
	double Q3_UCu_cal;

	double b_11D;
	double b_33D;
	double b_55D;
	double b_11U0;
	double b_33U0;
	double b_55U0;

	double deltaT;

	double H_r0D33;
	double H_r0U033;
	double H_r0UC33;
	
	Vector TRelVel;



	double m1;
	double m2;
	double m3;

	Vector H11;
	Vector H33;
	Vector H55;

	double mh1;
	double mh2;
	double mh3;

	double tol;
	double delta;

	double q1el;
	double q3el;
	double q5el;

	Vector t13;
	Vector t23;
	Vector t12;
	Vector t13_2;
	Vector t23_2;
	Vector t12_2;

	Vector Sum_pf;

	// --- trial variables ----------------
	Vector Tstress;
	Vector Tstrain;
	Vector cs;
	Vector O;

	int setFlows;

	int TrialStressElasticFlag;

	Vector Tc1;
	Vector Tc3;
	Vector Tc5;
	Vector Tc1Comm;
	Vector Tc3Comm;
	Vector Tc5Comm;

	Vector Qy1_iter;
	Vector Qy3_iter;
	Vector Qy5_iter;
	Vector Qy1_iter_2;
	Vector Qy3_iter_2;
	Vector Qy5_iter_2;

	int checkUnload;
	int unload;
	int numUnload;

	double U;
	double C;
	double dC;
	double dCC;
	int checkCC;

	Vector a_C;
	Vector b_C;

	Vector O_conv;

	Vector sum_q_p;

	Matrix Hs;
	Vector pf;
	Vector q_p;

	int NumYielded;
	// ------------------------------------
	

	// ---- committed variables -----------
	Vector Cstress;
	Vector Cstrain;
	Vector CO;
	Vector Ccs;
	Vector Csum_q_p;
	Vector CCsum_q_p;
	Matrix CHs;

	Vector Cpf;
	Vector CSum_pf;

	double CU;
	double CC;
	double CCC;
	double CMq;

	Vector Cq_p;

	int CnFlows;
	int PrnFlows;

	double PStr1;
	double PStr3;
	double PStr5;
	double PStr1_2;
	double PStr3_2;
	double PStr5_2;
	double TStr1_2;
	double TStr3_2;
	double TStr5_2;
	double CStr1_2;
	double CStr3_2;
	double CStr5_2;
	double Peps1;
	double Peps3;
	double Peps5;

	Vector CQy1;
	Vector CQy3;
	Vector CQy5;
	Vector CQy1_2;
	Vector CQy3_2;
	Vector CQy5_2;

	int PcheckUnload;
	int Punload;
	int PnumUnload;

	Vector PCQy1;
	Vector PCQy3;
	Vector PCQy5;
	Vector PCQy1_2;
	Vector PCQy3_2;
	Vector PCQy5_2;

	Vector Cc1;
	Vector Cc3;
	Vector Cc5;

	Matrix PKt;

	int CommittedStressElasticFlag;

	// ------------------------------------

	Matrix theTangent;

};

#endif
