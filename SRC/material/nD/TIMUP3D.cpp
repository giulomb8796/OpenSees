/*
Developed by: Davide Noè Gorini (davidenoe.gorini@unitn.it), Giuseppe Lombardi (giuseppe.lombardi@uniroma1.it)
Date: 2025
Description: this file contains the implementation of the TIM-UP simulating the three-axial force-deformation response at the foundation-superstructure contact of a shallow foundation and the soil interacting with it
Degrees of Freedom covered by the TIM-UP: 1 - horizontal translation, 3 - vertical translation, 5 - rotation around axis 2
Reference: Gorini, D.N., Lombardi, G., Callisto, L. (2025): "A thermodynamic-based macroelement framework for soil-foundation systems with coupled hydro-mechanical behaviour", Computers and Geotechnics, https://doi.org/10.1016/j.compgeo.2025.107681
*/

#include <math.h>
#include <stdlib.h>
#include <TIMUP3D.h>
#include <Information.h>
#include <ID.h>
#include <MaterialResponse.h>
#include <Parameter.h>

#include <elementAPI.h>
#include <string.h>
#include <algorithm>



void* OPS_TIMUP3D(void) {

	int tag, nYield, niter;
	double H11el_D, H33el_D, H55el_D, H11el_U0, H33el_U0, H55el_U0,
		H11el_UC, H33el_UC, H55el_UC, mHyperEl,	Q3min_D, Q3max_D, Q3min_U0,
		Q3max_U0, Q3min_UC, Q3max_UC, a1_DN, a5_DN, b_DN, r_DN,
		perm, Lfound, beta, Tf_min, Tf_max, a_U, Mq_min, Mq_max, dC0, tol, deltaT, 
		Q3_maxD_cal, Q3_maxU0_cal, Q3_maxUC_cal, Q3_Du_cal, Q3_U0u_cal, Q3_UCu_cal,
		b_11D = 70.1, b_33D = 1058.2, b_55D = 50.4, b_11U0 = 21.4, b_33U0 = 224.3, b_55U0 = 17.9;
		
	int numArgs = OPS_GetNumRemainingInputArgs();
	if (numArgs != 40 && numArgs != 46) {
		opserr << "ndMaterial TIMUP3D incorrect num args: want tag nYield "
			<< "H11el_D H33el_D H55el_D mHyperEl Q3min_D Q3max_D a1_DN a5_DN b_DN r_DN "
			<< "H11el_U0 H33el_U0 H55el_U0 Q3min_U0 Q3max_U0 "
			<< "H11el_UC H33el_UC H55el_UC Q3min_UC Q3max_UC "
			<< "perm Lfound beta Tf_min Tf_max a_U Mq_min Mq_max dC0 tol deltaT Q3_maxD_cal Q3_maxU0_cal Q3_maxUC_cal Q3_Du_cal Q3_U0u_cal Q3_UCu_cal niter"
			<< "<optional: b_11D b_33D b_55D b_11U0 b_33U0 b_55U0>\n";
		return 0;
	}

	int iData[2];
	double dData[37];

	int numData = 2;
	if (OPS_GetInt(&numData, iData) != 0) {
		opserr << "WARNING invalid integer values: nDMaterial TIMUP3D \n";
		return 0;
	}
	tag = iData[0];
	nYield = iData[1];

	if (nYield > 40) {
		opserr << "WARNING too many yield surfaces: nDMaterial TIMUP3D; nYield must be <=40 \n";
		return 0;
	}
	
	numData = 37;
	if (OPS_GetDouble(&numData, dData) != 0) {
		opserr << "WARNING invalid double values: nDMaterial TIMUP3D " << tag << endln;
		return 0;
	}
	H11el_D = dData[0];
	H33el_D = dData[1];
	H55el_D = dData[2];
	mHyperEl = dData[3];
	Q3min_D = dData[4];
	Q3max_D = dData[5];
	a1_DN = dData[6];
	a5_DN = dData[7];
	b_DN = dData[8];
	r_DN = dData[9];
    
	H11el_U0 = dData[10];
	H33el_U0 = dData[11];
	H55el_U0 = dData[12];
	Q3min_U0 = dData[13];
	Q3max_U0 = dData[14];

	H11el_UC = dData[15];
	H33el_UC = dData[16];
	H55el_UC = dData[17];
	Q3min_UC = dData[18];
	Q3max_UC = dData[19];

	perm = dData[20];
	Lfound = dData[21];
	beta = dData[22];

	Tf_min = dData[23];
	Tf_max = dData[24];
	a_U = dData[25];
	Mq_min = dData[26];
	Mq_max = dData[27];
	dC0 = dData[28];
	tol = dData[29];
	deltaT = dData[30];

	Q3_maxD_cal = dData[31];
	Q3_maxU0_cal = dData[32];
	Q3_maxUC_cal = dData[33];
	Q3_Du_cal = dData[34];
	Q3_U0u_cal = dData[35];
	Q3_UCu_cal = dData[36];

	// Optional parameters b111 and b333 (default values)
	b_11D = 70.1;
	b_33D = 1058.2;
	b_55D = 50.4;
	b_11U0 = 21.4;
	b_33U0 = 224.3;
	b_55U0 = 17.9;

	// Check if additional arguments are provided
	int numRemaining = OPS_GetNumRemainingInputArgs();

	if (numRemaining != 1 && numRemaining != 7) {
		opserr << "WARNING: optional parameters b_11D b_33D b_55D b_11U0 b_33U0 b_55U0 must be all provided or omitted\n";
		return nullptr;
	}

	if (numRemaining == 7) {
		double optData[6];
		int numDataOpt = 6;

		if (OPS_GetDouble(&numDataOpt, optData) != 0) {
			opserr << "WARNING: invalid optional inputs\n";
			return nullptr;
		}

		b_11D = optData[0];
		b_33D = optData[1];
		b_55D = optData[2];
		b_11U0 = optData[3];
		b_33U0 = optData[4];
		b_55U0 = optData[5];
	}

	numData = 1;
	int iterData[1];
	if (OPS_GetInt(&numData, iterData) != 0) {
		opserr << "WARNING invalid integer values: nDMaterial TIMUP3D \n";
		return 0;
	}
	niter = iterData[0];

	NDMaterial* theMaterial = new TIMUP3D(tag,
		nYield, H11el_D, H33el_D, H55el_D, mHyperEl,
		Q3min_D, Q3max_D, a1_DN, a5_DN, b_DN, r_DN,
		H11el_U0, H33el_U0, H55el_U0, Q3min_U0, Q3max_U0,
		H11el_UC, H33el_UC, H55el_UC, Q3min_UC, Q3max_UC,
		perm, Lfound, beta,	Tf_min, Tf_max, a_U, Mq_min, Mq_max, dC0, tol, deltaT,
		Q3_maxD_cal, Q3_maxU0_cal, Q3_maxUC_cal, Q3_Du_cal, Q3_U0u_cal, Q3_UCu_cal,
		niter,
		b_11D, b_33D, b_55D, b_11U0, b_33U0, b_55U0);
	
	return theMaterial;
}

TIMUP3D::TIMUP3D(int ptag, int pnYield,
	double pH11el_D, double pH33el_D, double pH55el_D, double pmHyperEl,
	double pQ3min_D, double pQ3max_D, double pa1_DN, double pa5_DN, double pb_DN, double pr_DN,
	double pH11el_U0, double pH33el_U0, double pH55el_U0,
	double pQ3min_U0, double pQ3max_U0,	double pH11el_UC, double pH33el_UC, double pH55el_UC,
	double pQ3min_UC, double pQ3max_UC, double pperm, double pLfound, double pbeta,
	double pTf_min, double pTf_max, double pa_U, double pMq_min, double pMq_max, double pdC0,
	double ptol, double pdeltaT, 
	double pQ3_maxD_cal, double pQ3_maxU0_cal, double pQ3_maxUC_cal, double pQ3_Du_cal, double pQ3_U0u_cal, double pQ3_UCu_cal,
	int pniter,
	double pb_11D, double pb_33D, double pb_55D, double pb_11U0, double pb_33U0, double pb_55U0)
	: NDMaterial(ptag, ND_TAG_TIMUP3D), Tstress(6),
	Tstrain(6), Cstress(6), Cstrain(6), theTangent(6, 6), TRelVel(6), Sum_pf(pnYield),
	a_C(3), b_C(2), Csum_q_p(6), CCsum_q_p(6), CHs(6,6), Cpf(pnYield), Cq_p(6*pnYield), Ccs(6*pnYield),
	a1(pnYield,3), a5(pnYield, 3), b(pnYield, 3), r(pnYield, 3), c10(pnYield, 3), Q3min(pnYield, 3), Q3max(pnYield, 3),
	c50(pnYield, 3), c30(pnYield, 3), fq3(pnYield, 3), Q30(pnYield, 3), m(pnYield), H_0D(6,6), H_0U0(6, 6), H_0UC(6, 6),
	H_eD(6*pnYield, 6), H_eU0(6*pnYield, 6), H_eUC(6*pnYield, 6), cs(6*pnYield), Hs(6,6), pf(pnYield),
	O_conv(pnYield * 6), CO(pnYield * 6), sum_q_p(6), q_p(6 * pnYield), CSum_pf(pnYield)
{
	this->nYield = pnYield;
	this->niter = pniter;
	this->tol = ptol;
	this->H11el_D = pH11el_D;
	this->H33el_D = pH33el_D;
	this->H55el_D = pH55el_D;
	this->Q3min_D = pQ3min_D;
	this->Q3max_D = pQ3max_D;
	this->a1_DN = pa1_DN;
	this->a5_DN = pa5_DN;
	this->b_DN = pb_DN;
	this->r_DN = pr_DN;

	this->Q3_maxD_cal = pQ3_maxD_cal;
	this->Q3_maxU0_cal = pQ3_maxU0_cal;
	this->Q3_maxUC_cal = pQ3_maxUC_cal;
	this->Q3_Du_cal = pQ3_Du_cal;
	this->Q3_U0u_cal = pQ3_U0u_cal;
	this->Q3_UCu_cal = pQ3_UCu_cal;

	this->b_11D = pb_11D;
	this->b_33D = pb_33D;
	this->b_55D = pb_55D;
	this->b_11U0 = pb_11U0;
	this->b_33U0 = pb_33U0;
	this->b_55U0 = pb_55U0;

	this->perm = pperm;
	this->Lfound = pLfound;
	this->beta = pbeta;

	this->H11el_U0 = pH11el_U0;
	this->H33el_U0 = pH33el_U0;
	this->H55el_U0 = pH55el_U0;
	this->Q3min_U0 = pQ3min_U0;
	this->Q3max_U0 = pQ3max_U0;
	this->a1_U0N = a1_DN*Q3max_U0 / Q3max_D;									// control the amplitude along the Q1 axis
	this->a5_U0N = a5_DN*Q3max_U0 / Q3max_D;									// control the amplitude along the Qr2 axis
	this->b_U0N = b_DN*Q3max_U0 / Q3max_D;										// dimension parameter (b / r control the shape)
	this->r_U0N = r_DN*Q3max_U0 / Q3max_D;										// dimension parameter (b / r control the shape)

	this->H11el_UC = pH11el_UC;
	this->H33el_UC = pH33el_UC;
	this->H55el_UC = pH55el_UC;
	this->Q3min_UC = pQ3min_UC;
	this->Q3max_UC = pQ3max_UC;
	this->a1_UCN = a1_U0N*Q3max_UC / Q3max_U0;									// control the amplitude along the Q1 axis
	this->a5_UCN = a5_U0N*Q3max_UC / Q3max_U0;									// control the amplitude along the Qr2 axis
	this->b_UCN = b_U0N*Q3max_UC / Q3max_U0;									// dimension parameter (b / r control the shape)
	this->r_UCN = r_U0N*Q3max_UC / Q3max_U0;									// dimension parameter (b / r control the shape)

	this->mHyperEl = pmHyperEl;

	this->NumYielded = 0;

	this->CnFlows = 0;
	this->PrnFlows = 0;
	this->PStr1 = 0.0;
	this->PStr3 = 0.0;
	this->PStr5 = 0.0;
	this->Peps1 = 0.0;
	this->Peps3 = 0.0;
	this->Peps5 = 0.0;
	this->PStr1_2 = 0.0;
	this->PStr3_2 = 0.0;
	this->PStr5_2 = 0.0;
	this->TStr1_2 = 0.0;
	this->TStr3_2 = 0.0;
	this->TStr5_2 = 0.0;
	this->CStr1_2 = 0.0;
	this->CStr3_2 = 0.0;
	this->CStr5_2 = 0.0;

	this->U = 1.;
	this->C = 0.0;
	this->dC = 0.;
	this->dCC = 0.;
	this->checkCC = 0;

	if (perm >= 1.) {
		this->U = 1.;
	}
	else if (perm <= 1.e-10) {
		this->U = 0.;
	}

	this->CU = U;

	double Uehi = U;

	PcheckUnload = 0;
	Punload = 0;
	PnumUnload = 0;

	Tstress.Zero();
	Tstrain.Zero();
	
	TRelVel.Zero();

	H_0UC.Zero();
	H_0U0.Zero();

	Cstress.Zero();
	Cstrain.Zero();

	CSum_pf.Zero();
	Csum_q_p.Zero();
	CCsum_q_p.Zero();
	Cpf.Zero();
	Cq_p.Zero();
	q_p.Zero();
	Ccs.Zero();
	cs.Zero();
	pf.Zero();
	O_conv.Zero();
	CO.Zero();
	sum_q_p.Zero();
	
	H_eD.Zero();
	H_eU0.Zero();
	H_eUC.Zero();
	
	this->CC = 0;
	this->CCC = 0;

	this->fq3_DN = (Q3max_D - Q3min_D) / (2. * r_DN);								// scale factor along Q3 axis
	this->fq3_U0N = (Q3max_U0 - Q3min_U0) / (2. * r_U0N);							// scale factor along Q3 axis
	this->fq3_UCN = (Q3max_UC - Q3min_UC) / (2. * r_UCN);							// scale factor along Q3 axis

	this->Q30_DN = Q3max_D - (Q3max_D - Q3min_D)*(b_DN / r_DN + 1.) / 2.;			// traslation along Q3 axis
	this->Q30_U0N = Q3max_U0 - (Q3max_U0 - Q3min_U0)*(b_U0N / r_U0N + 1.) / 2.;		// traslation along Q3 axis
	this->Q30_UCN = Q3max_UC - (Q3max_UC - Q3min_UC)*(b_UCN / r_UCN + 1.) / 2.;		// traslation along Q3 axis

	double c10_DN = 0;																// traslation respect to the origin(Q1 axis)
	double c50_DN = 0;																// traslation respect to the origin(Qr2 axis)
	double c30_DN = 0;																// traslation respect to the origin(Q3 axis)
	double c10_U0N = 0;																// traslation respect to the origin(Q1 axis)
	double c50_U0N = 0;																// traslation respect to the origin(Qr2 axis)
	double c30_U0N = 0;																// traslation respect to the origin(Q3 axis)
	double c10_UCN = 0;																// traslation respect to the origin(Q1 axis)
	double c50_UCN = 0;																// traslation respect to the origin(Qr2 axis)
	double c30_UCN = 0;																// traslation respect to the origin(Q3 axis)

	this->a1(nYield - 1, 0) = a1_DN;
	this->a5(nYield - 1, 0) = a5_DN;
	this->b(nYield - 1, 0) = b_DN;
	this->r(nYield - 1, 0) = r_DN;
	this->c10(nYield - 1, 0) = c10_DN;
	this->Q3min(nYield - 1, 0) = Q3min_D;
	this->Q3max(nYield - 1, 0) = Q3max_D;
	this->c50(nYield - 1, 0) = c50_DN;
	this->c30(nYield - 1, 0) = c30_DN;
	this->fq3(nYield - 1, 0) = fq3_DN;
	this->Q30(nYield - 1, 0) = Q30_DN;

	this->a1(nYield - 1, 1) = a1_U0N;
	this->a5(nYield - 1, 1) = a5_U0N;
	this->b(nYield - 1, 1) = b_U0N;
	this->r(nYield - 1, 1) = r_U0N;
	this->c10(nYield - 1, 1) = c10_U0N;
	this->Q3min(nYield - 1, 1) = Q3min_U0;
	this->Q3max(nYield - 1, 1) = Q3max_U0;
	this->c50(nYield - 1, 1) = c50_U0N;
	this->c30(nYield - 1, 1) = c30_U0N;
	this->fq3(nYield - 1, 1) = fq3_U0N;
	this->Q30(nYield - 1, 1) = Q30_U0N;

	this->a1(nYield - 1, 2) = a1_UCN;
	this->a5(nYield - 1, 2) = a5_UCN;
	this->b(nYield - 1, 2) = b_UCN;
	this->r(nYield - 1, 2) = r_UCN;
	this->c10(nYield - 1, 2) = c10_UCN;
	this->Q3min(nYield - 1, 2) = Q3min_UC;
	this->Q3max(nYield - 1, 2) = Q3max_UC;
	this->c50(nYield - 1, 2) = c50_UCN;
	this->c30(nYield - 1, 2) = c30_UCN;
	this->fq3(nYield - 1, 2) = fq3_UCN;
	this->Q30(nYield - 1, 2) = Q30_UCN;

	// reduction coefficent
	double del = 0.1; // reduction factor of the size of the first yield surface compared to the Ultimate Surface size

	this->a1(0, 0) = a1_DN*del;
	this->a5(0, 0) = a5_DN*del;
	this->b(0, 0) = b_DN*del;
	this->r(0, 0) = r_DN*del;
	this->c10(0, 0) = c10_DN*del;
	this->Q3min(0, 0) = Q3min_D*del;
	this->Q3max(0, 0) = Q3max_D*del;
	this->c50(0, 0) = c50_DN*del;
	this->c30(0, 0) = c30_DN*del;
	this->fq3(0, 0) = (Q3max(0, 0) - Q3min(0, 0)) / (2. * r(0, 0));
	this->Q30(0, 0) = Q3max(0, 0) - (Q3max(0, 0) - Q3min(0, 0))*(b(0, 0) / r(0, 0) + 1.) / 2.;

	this->a1(0, 1) = a1_U0N*del;
	this->a5(0, 1) = a5_U0N*del;
	this->b(0, 1) = b_U0N*del;
	this->r(0, 1) = r_U0N*del;
	this->c10(0, 1) = c10_U0N*del;
	this->Q3min(0, 1) = Q3min_U0*del;
	this->Q3max(0, 1) = Q3max_U0*del;
	this->c50(0, 1) = c50_U0N*del;
	this->c30(0, 1) = c30_U0N*del;
	this->fq3(0, 1) = (Q3max(0, 1) - Q3min(0, 1)) / (2. * r(0, 1));
	this->Q30(0, 1) = Q3max(0, 1) - (Q3max(0, 1) - Q3min(0, 1))*(b(0, 1) / r(0, 1) + 1.) / 2.;

	this->a1(0, 2) = a1_UCN*del;
	this->a5(0, 2) = a5_UCN*del;
	this->b(0, 2) = b_UCN*del;
	this->r(0, 2) = r_UCN*del;
	this->c10(0, 2) = c10_UCN*del;
	this->Q3min(0, 2) = Q3min_UC*del;
	this->Q3max(0, 2) = Q3max_UC*del;
	this->c50(0, 2) = c50_UCN*del;
	this->c30(0, 2) = c30_UCN*del;
	this->fq3(0, 2) = (Q3max(0, 2) - Q3min(0, 2)) / (2. * r(0, 2));
	this->Q30(0, 2) = Q3max(0, 2) - (Q3max(0, 2) - Q3min(0, 2))*(b(0, 2) / r(0, 2) + 1.) / 2.;

	this->m(0) = del;
	for (int j = 0; j < 3; j++)
	{
		for (int jj = 1; jj < nYield; jj++)
		{
			this->m(jj) = del + jj / (nYield - 1.) * (1. - del);
			this->a1(jj, j) = m(jj)*a1(nYield - 1, j);
			this->a5(jj, j) = m(jj)*a5(nYield - 1, j);
			this->b(jj, j) = m(jj)*b(nYield - 1, j);
			this->r(jj, j) = m(jj)*r(nYield - 1, j);
			this->c10(jj, j) = m(jj)*c10(nYield - 1, j);
			this->Q3min(jj, j) = m(jj)*Q3min(nYield - 1, j);
			this->Q3max(jj, j) = m(jj)*Q3max(nYield - 1, j);
			this->c50(jj, j) = m(jj)*c50(nYield - 1, j);
			this->c30(jj, j) = m(jj)*c30(nYield - 1, j);
			this->fq3(jj, j) = (Q3max(jj, j) - Q3min(jj, j)) / (2 * r(jj, j));
			this->Q30(jj, j) = Q3max(jj, j) - (Q3max(jj, j) - Q3min(jj, j))*(b(jj, j) / r(jj, j) + 1.) / 2.;
		}
	}

	// iperbolic law to describe the load - displacement relationahip for the selected path:
	// Q_i = q_i / (a_i + q_i / b_i)
	double a_11D = 1. / H11el_D;
	double a_33D = 1. / H33el_D;
	double a_55D = 1. / H55el_D;
	double a_33U0 = 1. / H33el_U0;
	double a_11U0 = 1. / H11el_U0;
	double a_55U0 = 1. / H55el_U0;
	double a_33UC = 1. / H33el_UC;
	double b_33UC = b_33U0;
	double a_11UC = 1. / H33el_UC;
	double b_11UC = b_11U0;
	double a_55UC = 1. / H55el_UC;
	double b_55UC = b_55U0;

	// elastic stiffness
	// evaluation of the generalized force for the first surface
	double chi3_cal = 0.;
	if (Q3_maxD_cal > Q3max(0, 0)) {
		chi3_cal = Q3max(0, 0) - Q3_Du_cal;
	}
	else {
		chi3_cal = Q3_maxD_cal - Q3_Du_cal;
	}
	// loads used for the calibration
	Vector Q1_calD(nYield), q1_calD(nYield),
		Q3_calD(nYield), q3_calD(nYield),
		Q5_calD(nYield), q5_calD(nYield);

	double ao = a1(0, 0);
	double Q1_calD1 = (pow(r(0, 0), 2.0) - pow((chi3_cal - Q30(0, 0)) / fq3(0, 0) - b(0, 0), 2.0));
	double Q1_calD2 = pow(2. * b(0, 0) - (chi3_cal - Q30(0, 0)) / fq3(0, 0), 2.);
	double Q1_calD3 = Q1_calD1 / Q1_calD2;
	Q1_calD(0) = ao * pow(Q1_calD3, 0.5);
	q1_calD(0) = Q1_calD(0)*a_11D / (1. - Q1_calD(0) / b_11D);
	Q3_calD(0) = Q3max(0, 0);
	q3_calD(0) = Q3_calD(0)*a_33D / (1 - Q3_calD(0) / b_33D);
	Q5_calD(0) = a5(0, 0)*pow((pow(r(0, 0), 2.0) - pow((chi3_cal - Q30(0, 0)) / fq3(0, 0) - b(0, 0), 2.0)) / pow(2. * b(0, 0) - (chi3_cal - Q30(0, 0)) / fq3(0, 0), 2.), 0.5);
	q5_calD(0) = Q5_calD(0)*a_55D / (1 - Q5_calD(0) / b_55D);

	H_0D.Zero();
	H_0D(0, 0) = H11el_D;
	H_0D(2, 2) = H33el_D;
	H_0D(4, 4) = H55el_D;

	// undrained not degraded case
	// surfaces ceners inizialization
	Vector O_3U0(nYield);
	// evaluation of the generalized force for the first surface
	if (Q3_maxU0_cal > Q3max(0, 0)) {
		O_3U0(0) = Q3_maxU0_cal - Q3max(0, 0) / 2.;
		if (O_3U0(0) + Q3max(0, 1) / 2. > Q3_maxU0_cal - Q3_U0u_cal) {
			chi3_cal = Q3_maxU0_cal - Q3_U0u_cal - (O_3U0(0) - Q3max(0, 1) / 2.);
		}
		else if (O_3U0(0) + Q3max(0, 1) / 2. <= Q3_maxU0_cal - Q3_U0u_cal) {
			chi3_cal = 0.999*Q3max(0, 1);
		}
	}
	else if (Q3_maxU0_cal < Q3max(0, 0) && Q3_maxU0_cal - Q3_U0u_cal > Q3max(0, 1)) {
		O_3U0(0) = Q3_maxU0_cal - Q3_U0u_cal - Q3max(0, 0) / 2.;
		chi3_cal = 0.999*Q3max(0, 1);
	}
	else {
		chi3_cal = Q3_maxU0_cal - Q3_U0u_cal;
	}

	// loads used for the calibration
	Vector Q1_calU0(nYield), q1_calU0(nYield),
		Q3_calU0(nYield), q3_calU0(nYield),
		Q5_calU0(nYield), q5_calU0(nYield);

	Q1_calU0(0) = a1(0, 1)*pow((pow(r(0, 1), 2.0) - pow((chi3_cal - Q30(0, 1)) / fq3(0, 1) - b(0, 1), 2.0)) / pow(2. * b(0, 1) - (chi3_cal - Q30(0, 1)) / fq3(0, 1), 2.), 0.5);
	q1_calU0(0) = Q1_calU0(0)*a_11U0 / (1. - Q1_calU0(0) / b_11U0);
	Q3_calU0(0) = Q3max(0, 1);
	q3_calU0(0) = Q3_calU0(0)*a_33U0 / (1 - Q3_calU0(0) / b_33U0);
	Q5_calU0(0) = a5(0, 1)*pow((pow(r(0, 1), 2.0) - pow((chi3_cal - Q30(0, 1)) / fq3(0, 1) - b(0, 1), 2.0)) / pow(2. * b(0, 1) - (chi3_cal - Q30(0, 1)) / fq3(0, 1), 2.), 0.5);
	q5_calU0(0) = Q5_calU0(0)*a_55U0 / (1 - Q5_calU0(0) / b_55U0);

	H_0U0(0, 0) = H11el_U0;
	H_0U0(2, 2) = H33el_U0;
	H_0U0(4, 4) = H55el_U0;
	H_0UC(0, 0) = H11el_U0;
	H_0UC(2, 2) = H33el_U0;
	H_0UC(4, 4) = H55el_U0;

	// edits for undrained
	Hs = CHs = H_0D*U + (1 - U)*H_0U0;
	// edits for undrained

	// elastic vertical stifness evolution with the load
	double factorRot = 2.2; // Gorini et al. (2025) - https://doi.org/10.1016/j.compgeo.2025.107681
	this->H_r0D33 = H33el_D*factorRot;
	this->H_r0U033 = H33el_U0*factorRot;
	this->H_r0UC33 = H33el_UC*factorRot;
	
	// nth stiffness matrix 
	Matrix H_pD(6 * nYield, 6), H_pU0(6 * nYield, 6), H_pUC(6 * nYield, 6);
	H_pD.Zero();
	H_pU0.Zero();
	H_pUC.Zero();
	// elastic-plastic tangent stiffness matrix
	double factorUlt = 0.03; // Gorini et al. (2025) - https://doi.org/10.1016/j.compgeo.2025.107681
	for (int jj = 0; jj < 6; jj++)
	{
		H_pD(6 * nYield - (5 - jj) - 1, jj) = factorUlt*H_0D(jj, jj);
		H_pU0(6 * nYield - (5 - jj) - 1, jj) = factorUlt*H_0U0(jj, jj);
		H_pUC(6 * nYield - (5 - jj) - 1, jj) = factorUlt*H_0UC(jj, jj);
	}

	for (int n = 1; n < nYield; n++)
	{
		// evaluation of the dissipative force vector associated with the nth yield surface
		if (Q3_maxD_cal > Q3max(n,0)) {
			chi3_cal = Q3max(n, 0) - Q3_Du_cal;
		}
		else {
			chi3_cal = Q3_maxD_cal - Q3_Du_cal;
		}
		// forces and deformations used for the initialisation of the stiffness matrices associated with the yield surfaces
		Q1_calD(n) = a1(n, 0)*pow((pow(r(n, 0), 2.0) - pow((chi3_cal - Q30(n, 0)) / fq3(n, 0) - b(n, 0), 2.0)) / pow(2. * b(n, 0) - (chi3_cal - Q30(n, 0)) / fq3(n, 0), 2.), 0.5);
		q1_calD(n) = Q1_calD(n)*a_11D / (1. - Q1_calD(n) / b_11D);
		Q3_calD(n) = Q3max(n, 0);
		q3_calD(n) = Q3_calD(n)*a_33D / (1 - Q3_calD(n) / b_33D);
		Q5_calD(n) = a5(n, 0)*pow((pow(r(n, 0), 2.0) - pow((chi3_cal - Q30(n, 0)) / fq3(n, 0) - b(n, 0), 2.0)) / pow(2. * b(n, 0) - (chi3_cal - Q30(n, 0)) / fq3(n, 0), 2.), 0.5);
		q5_calD(n) = Q5_calD(n)*a_55D / (1 - Q5_calD(n) / b_55D);

		// elastic-plastic tangent stiffness matrix
		H_pD(6 * n - 5 - 1, 0) = (Q1_calD(n) - Q1_calD(n-1)) / (q1_calD(n) - q1_calD(n-1));
		H_pD(6 * n - 3 - 1, 2) = (Q3_calD(n) - Q3_calD(n-1)) / (q3_calD(n) - q3_calD(n-1));
		H_pD(6 * n - 1 - 1, 4) = (Q5_calD(n) - Q5_calD(n-1)) / (q5_calD(n) - q5_calD(n-1));
		// nth stiffness matrix
		Matrix H_eD_Q(6,6);
		H_eD_Q = H_0D;
		H_eD_Q(2, 2) = ((pow((2. * mHyperEl - 1),2.0)) / (2 * mHyperEl))*(pow((2. * H_r0D33),(mHyperEl / (2. * mHyperEl - 1))))*pow(Q3max(n, 0),((2 * mHyperEl - 2) / (2 * mHyperEl - 1)));

		// nth compliant matrix
		Matrix C_eD(6, 6);
		C_eD.Zero();
		C_eD(0, 0) = 1./H_eD_Q(0, 0);
		C_eD(2, 2) = 1./H_eD_Q(2, 2);
		C_eD(4, 4) = 1./H_eD_Q(4, 4);

		if (n > 1) {
			for (int j = 1; j < n; j++)
			{
				C_eD(0, 0) = C_eD(0, 0) + 1. /(H_eD(6 * j - 5 - 1, 0));
				C_eD(2, 2) = C_eD(2, 2) + 1. /(H_eD(6 * j - 3 - 1, 2));
				C_eD(4, 4) = C_eD(4, 4) + 1. /(H_eD(6 * j - 1 - 1, 4));
			}
		}
		C_eD(0, 0) = 1. /(H_pD(6 * n - 5 - 1, 0)) - C_eD(0, 0);
		C_eD(2, 2) = 1. /(H_pD(6 * n - 3 - 1, 2)) - C_eD(2, 2);
		C_eD(4, 4) = 1. /(H_pD(6 * n - 1 - 1, 4)) - C_eD(4, 4);
        // nth stiffness matrix
		H_eD(6 * n - 5 - 1, 0) = 1./C_eD(0, 0);
		H_eD(6 * n - 3 - 1, 2) = 1./C_eD(2, 2);
		H_eD(6 * n - 1 - 1, 4) = 1./C_eD(4, 4);

		// undrained, non-degraded condition
		// evaluation of the dissipative force for the nth yield surface
		if (Q3_maxU0_cal >= Q3max(n, 0)) {
			O_3U0(n) = Q3_maxU0_cal - Q3max(n, 0) / 2.;
		    if (O_3U0(n) + Q3max(n, 1) / 2. < Q3_maxU0_cal - Q3_U0u_cal) {
			   O_3U0(n) = O_3U0(0) + Q3max(0, 1) / 2. - Q3max(n, 1) / 2.;
			}
		    chi3_cal = (Q3_maxU0_cal - Q3_U0u_cal) - (O_3U0(n) - Q3max(n, 1) / 2.);
		}
		else if ((Q3_maxU0_cal < Q3max(n + 1, 0)) && ((Q3_maxU0_cal - Q3_U0u_cal) > Q3max(n, 1))) {
			O_3U0(n) = O_3U0(n-1) + Q3max(n-1, 1) / 2. - Q3max(n, 1) / 2.;
		    chi3_cal = (Q3_maxU0_cal - Q3_U0u_cal) - (O_3U0(n) - Q3max(n, 1) / 2.);
		}
		else {
			chi3_cal = Q3max(n, 1) - (Q3_maxU0_cal - Q3_U0u_cal);
		}

		// forces and deformations used for the initialisation of the stiffness matrices associated with the yield surfaces
        Q1_calU0(n) = a1(n, 1)*pow((pow(r(n, 1), 2.0) - pow((chi3_cal - Q30(n, 1)) / fq3(n, 1) - b(n, 1), 2.0)) / pow(2. * b(n, 1) - (chi3_cal - Q30(n, 1)) / fq3(n, 1), 2.), 0.5);
		q1_calU0(n) = Q1_calU0(n)*a_11U0 / (1. - Q1_calU0(n) / b_11U0);
		Q3_calU0(n) = Q3max(n, 1);
		q3_calU0(n) = Q3_calU0(n)*a_33U0 / (1 - Q3_calU0(n) / b_33U0);
		Q5_calU0(n) = a5(n, 1)*pow((pow(r(n, 1), 2.0) - pow((chi3_cal - Q30(n, 1)) / fq3(n, 1) - b(n, 1), 2.0)) / pow(2. * b(n, 1) - (chi3_cal - Q30(n, 1)) / fq3(n, 1), 2.), 0.5);
		q5_calU0(n) = Q5_calU0(n)*a_55U0 / (1 - Q5_calU0(n) / b_55U0);

		// elastic-plastic tangent stiffness matrix
		H_pU0(6 * n - 5 - 1, 0) = (Q1_calU0(n) - Q1_calU0(n-1)) / (q1_calU0(n) - q1_calU0(n-1));      
		H_pU0(6 * n - 3 - 1, 2) = (Q3_calU0(n) - Q3_calU0(n-1)) / (q3_calU0(n) - q3_calU0(n-1));
		H_pU0(6 * n - 1 - 1, 4) = (Q5_calU0(n) - Q5_calU0(n-1)) / (q5_calU0(n) - q5_calU0(n-1));
	
        // nth stiffness matrix
		Matrix H_eU0_Q;
		H_eU0_Q = H_0U0;
		H_eU0_Q(2, 2) = ((pow((2 * mHyperEl - 1),2.0)) / (2. * mHyperEl))*(pow((2. * H_r0U033),(mHyperEl / (2. * mHyperEl - 1))))*pow(Q3max(n, 1),((2. * mHyperEl - 2) / (2. * mHyperEl - 1)));

		// nth compliant matrix
		Matrix C_eU0(6,6);
		C_eU0.Zero();
		C_eU0(0, 0) = 1. /H_eU0_Q(0, 0);
		C_eU0(2, 2) = 1. /H_eU0_Q(2, 2);
		C_eU0(4, 4) = 1. /H_eU0_Q(4, 4);
		if (n > 1) {
			for (int j = 1; j < n; j++)
			{
				C_eU0(0, 0) = C_eU0(0, 0) + 1. /(H_eU0(6 * j - 5 - 1, 0));
				C_eU0(2, 2) = C_eU0(2, 2) + 1. /(H_eU0(6 * j - 3 - 1, 2));
				C_eU0(4, 4) = C_eU0(4, 4) + 1. /(H_eU0(6 * j - 1 - 1, 4));
			}
		}
		C_eU0(0, 0) = 1. /(H_pU0(6 * n - 5 - 1, 0)) - C_eU0(0, 0);
		C_eU0(2, 2) = 1. /(H_pU0(6 * n - 3 - 1, 2)) - C_eU0(2, 2);
		C_eU0(4, 4) = 1. /(H_pU0(6 * n - 1 - 1, 4)) - C_eU0(4, 4);

        // nth stiffness matrix
		H_eU0(6 * n - 5 - 1, 0) = 1. /C_eU0(0, 0);
		H_eU0(6 * n - 3 - 1, 2) = 1. /C_eU0(2, 2);
		H_eU0(6 * n - 1 - 1, 4) = 1. /C_eU0(4, 4);

		// undrained, degraded condition
		H_eUC(6 * n - 5 - 1, 0) = H_eU0(6 * n - 5 - 1, 0);
		H_eUC(6 * n - 3 - 1, 2) = H_eU0(6 * n - 3 - 1, 2);
		H_eUC(6 * n - 1 - 1, 4) = H_eU0(6 * n - 1 - 1, 4);	
	}
	// last nth stiffness matrix (low % of the elastic stiffness)
	for (int jj = 0; jj < 6; jj++)
	{
		H_eD(6 * nYield - (5 - jj) - 1, jj) = factorUlt*H_0D(jj, jj);
		H_eU0(6 * nYield - (5 - jj) - 1, jj) = factorUlt*H_0U0(jj, jj);
		H_eUC(6 * nYield - (5 - jj) - 1, jj) = factorUlt*H_0UC(jj, jj);
	}

	// drainage coefficent
	// C coefficent
	// Bc's:
	// dC = dC0 if Mq = Mq_min, C = 0
	// dC = 0 if Mq = Mq_max, C = 0
	// ddC/dMq = 0 if Mq = Mq_max,C = 0
	// dC = 0 if Mq = Mq_min, C = 0
	// ddC/dC = 0 if Mq = Mq_min, C = 1
	Mq_min = pMq_min;
	Mq_max = pMq_max;
	dC0 = pdC0;
	// equation system matrix
	Matrix M_C(3,3);
	M_C(0, 0) = 1;
	M_C(0, 1) = Mq_min;
	M_C(0, 2) = pow(Mq_min,2.);
	M_C(1, 0) = 1;
	M_C(1, 1) = Mq_max;
	M_C(1, 2) = pow(Mq_max,2.);
	M_C(2, 2) = 0;
	M_C(2, 2) = 1;
	M_C(2, 2) = 2. * Mq_max;
	// coefficents vectors
	// Vector a_C(3), b_C(2);
	a_C.Zero();
	b_C.Zero();

	a_C(0) = dC0 * (pow(Mq_min, 2.) - 2. * Mq_min * Mq_max) / pow(Mq_max - Mq_min,2.);
	a_C(1) = 2. * dC0 * Mq_max / pow(Mq_max - Mq_min, 2.);
	a_C(2) = -dC0 / pow(Mq_max - Mq_min, 2.);
	b_C(0) = 2.;
	b_C(1) = -1.;

	// U coefficent
	// Bc's: U = 0 if Tf = 0; U = U* if Tf = Tf*
	Tf_min = pTf_min;
	Tf_max = pTf_max;
	// coefficents vector
	a_U = 2.997e-5; // Gorini et al. (2025) - https://doi.org/10.1016/j.compgeo.2025.107681
	this->Tf0;
	Tf0 = -1. / a_U;

	this->deltaT = pdeltaT;

	// initial elastic-plastic stiffness
	Matrix H0(6,6), Cs(6, 6), Hs_dummy(3, 3), Cs_dummy(3, 3);
	H0.Zero();
	Hs.Zero();
	Cs.Zero();
	Hs_dummy.Zero();
	Cs_dummy.Zero();

	Hs = H_0D;

	Hs_dummy(0,0) = H_0D(0, 0);
	Hs_dummy(1,1) = H_0D(2,2);
	Hs_dummy(2,2) = H_0D(4, 4);

	// initial elastic - plastic compliance of the system
	Hs_dummy.Invert(Cs_dummy);

	Cs(0, 0) = Cs_dummy(0, 0);
	Cs(2, 2) = Cs_dummy(1, 1);
	Cs(4, 4) = Cs_dummy(2, 2);

	// initial surface coefficents
	Vector a1_ev(nYield), ar2_ev(nYield), b_ev(nYield), r_ev(nYield), fq3_ev(nYield), Q30_ev(nYield);

	for (int n = 0; n < nYield; n++)
	{
	a1_ev(n) = a1(n, 0);
	ar2_ev(n) = a5(n, 0);
	b_ev(n) = b(n, 0);
	r_ev(n) = r(n, 0);
	fq3_ev(n) = fq3(n, 0);
	Q30_ev(n) = Q30(n, 0);
	}

	this->beta = 0.0;

	theTangent.Zero();
	theTangent(0, 0) = H11el_D;
	theTangent(2, 2) = H33el_D;
	theTangent(4, 4) = H55el_D;

	this->TrialStressElasticFlag = 0;
	this->CommittedStressElasticFlag = 0;

	Sum_pf.Zero();

}

// destructor
TIMUP3D::TIMUP3D() : NDMaterial(0, ND_TAG_TIMUP3D), 
Tstress(6), Tstrain(6), Cstress(6), Cstrain(6), theTangent(6,6), TRelVel(6), Sum_pf(nYield),
a_C(3), b_C(2), Csum_q_p(6), CCsum_q_p(6), CHs(6, 6), Cpf(nYield), Cq_p(6 * nYield), Ccs(6 * nYield),
a1(nYield, 3), a5(nYield, 3), b(nYield, 3), r(nYield, 3), c10(nYield, 3), Q3min(nYield, 3), Q3max(nYield, 3),
c50(nYield, 3), c30(nYield, 3), fq3(nYield, 3), Q30(nYield, 3), m(nYield), H_0D(6, 6), H_0U0(6, 6), H_0UC(6, 6),
H_eD(6 * nYield, 6), H_eU0(6 * nYield, 6), H_eUC(6 * nYield, 6), cs(6 * nYield), Hs(6,6), pf(nYield),
O_conv(nYield * 6), CO(nYield * 6), sum_q_p(6), q_p(6 * nYield), CSum_pf(nYield)
{
	this->nYield = 0.0;
	this->niter = 0.0;
	this->tol = 0.0;
	this->H11el_D = 0.0;
	this->H33el_D = 0.0;
	this->H55el_D = 0.0;
	this->Q3min_D = 0.0;
	this->Q3max_D = 0.0;
	this->a1_DN = 0.0;
	this->a5_DN = 0.0;
	this->b_DN = 0.0;
	this->r_DN = 0.0;

	this->Q3_maxD_cal = 0.;
	this->Q3_maxU0_cal = 0.;
	this->Q3_maxUC_cal = 0.;
	this->Q3_Du_cal = 0.;
	this->Q3_U0u_cal = 0.;
	this->Q3_UCu_cal = 0.;

	this->b_11D = 0.;
	this->b_33D = 0.;
	this->b_55D = 0.;
	this->b_11U0 = 0.;
	this->b_33U0 = 0.;
	this->b_55U0 = 0.;

	this->NumYielded = 0;

	this->perm = 0.;
	this->Lfound = 0.;

	this->H11el_U0 = 0.0;
	this->H33el_U0 = 0.0;
	this->H55el_U0 = 0.0;
	this->Q3min_U0 = 0.0;
	this->Q3max_U0 = 0.0;
	this->a1_U0N = 0.0;
	this->a5_U0N = 0.0;
	this->b_U0N = 0.0;
	this->r_U0N = 0.0;

	this->H11el_UC = 0.0;
	this->H33el_UC = 0.0;
	this->H55el_UC = 0.0;
	this->Q3min_UC = 0.0;
	this->Q3max_UC = 0.0;
	this->a1_UCN = 0.0;
	this->a5_UCN = 0.0;
	this->b_UCN = 0.0;
	this->r_UCN = 0.0;

	this->mHyperEl = 0.0;

	this->CnFlows = 0;
	this->PrnFlows = 0;
	this->PStr1 = 0.0;
	this->PStr3 = 0.0;
	this->PStr5 = 0.0;
	this->Peps1 = 0.0;
	this->Peps3 = 0.0;
	this->Peps5 = 0.0;
	this->PStr1_2 = 0.0;
	this->PStr3_2 = 0.0;
	this->PStr5_2 = 0.0;
	this->TStr1_2 = 0.0;
	this->TStr3_2 = 0.0;
	this->TStr5_2 = 0.0;
	this->CStr1_2 = 0.0;
	this->CStr3_2 = 0.0;
	this->CStr5_2 = 0.0;

	this->U = 1.;
	this->C = 0.;
	this->dC = 0.;
	this->dCC = 0.;
	this->checkCC = 0;

	this->CU = 1.;

	CQy1.Zero();
	CQy3.Zero();
	CQy5.Zero();
	CQy1_2.Zero();
	CQy3_2.Zero();
	CQy5_2.Zero();

	CSum_pf.Zero();

	PcheckUnload = 0;
	Punload = 0;
	PnumUnload = 0;

	Qy1_iter.Zero();
	Qy3_iter.Zero();
	Qy5_iter.Zero();
	Qy1_iter_2.Zero();
	Qy3_iter_2.Zero();
	Qy5_iter_2.Zero();
	this->setFlows = 0;

	PCQy1.Zero();
	PCQy3.Zero();
	PCQy5.Zero();
	PCQy1_2.Zero();
	PCQy3_2.Zero();
	PCQy5_2.Zero();

	H_eD.Zero();
	H_eU0.Zero();
	H_eUC.Zero();

	Tstress.Zero();
	Tstrain.Zero();

	Cstress.Zero();
	Cstrain.Zero();

	TRelVel.Zero();

	Hs.Zero();

	Csum_q_p.Zero();
	CCsum_q_p.Zero();
	CHs.Zero();
	Cpf.Zero();
	Cq_p.Zero();
	q_p.Zero();
	Ccs.Zero();
	cs.Zero();
	pf.Zero();
	O_conv.Zero();
	CO.Zero();
	sum_q_p.Zero();

	this->fq3_DN = 0.0;
	this->fq3_U0N = 0.0;
	this->fq3_UCN = 0.0;

	this->Q30_DN = 0.0;
	this->Q30_U0N = 0.0;
	this->Q30_UCN = 0.0;

	double c10_DN = 0;
	double c50_DN = 0;
	double c30_DN = 0;
	double c10_U0N = 0;
	double c50_U0N = 0;
	double c30_U0N = 0;
	double c10_UCN = 0;
	double c50_UCN = 0;
	double c30_UCN = 0;

	this->a1(nYield - 1, 0) = 0.0;
	this->a5(nYield - 1, 0) = 0.0;
	this->b(nYield - 1, 0) = 0.0;
	this->r(nYield - 1, 0) = 0.0;
	this->c10(nYield - 1, 0) = 0.0;
	this->Q3min(nYield - 1, 0) = 0.0;
	this->Q3max(nYield - 1, 0) = 0.0;
	this->c50(nYield - 1, 0) = 0.0;
	this->c30(nYield - 1, 0) = 0.0;
	this->fq3(nYield - 1, 0) = 0.0;
	this->Q30(nYield - 1, 0) = 0.0;

	this->a1(nYield - 1, 1) = 0.0;
	this->a5(nYield - 1, 1) = 0.0;
	this->b(nYield - 1, 1) = 0.0;
	this->r(nYield - 1, 1) = 0.0;
	this->c10(nYield - 1, 1) = 0.0;
	this->Q3min(nYield - 1, 1) = 0.0;
	this->Q3max(nYield - 1, 1) = 0.0;
	this->c50(nYield - 1, 1) = 0.0;
	this->c30(nYield - 1, 1) = 0.0;
	this->fq3(nYield - 1, 1) = 0.0;
	this->Q30(nYield - 1, 1) = 0.0;

	this->a1(nYield - 1, 2) = 0.0;
	this->a5(nYield - 1, 2) = 0.0;
	this->b(nYield - 1, 2) = 0.0;
	this->r(nYield - 1, 2) = 0.0;
	this->c10(nYield - 1, 2) = 0.0;
	this->Q3min(nYield - 1, 2) = 0.0;
	this->Q3max(nYield - 1, 2) = 0.0;
	this->c50(nYield - 1, 2) = 0.0;
	this->c30(nYield - 1, 2) = 0.0;
	this->fq3(nYield - 1, 2) = 0.0;
	this->Q30(nYield - 1, 2) = 0.0;

	// reduction coefficent
	double del = 0.0;
	this->a1(0, 0) = 0.0;
	this->a5(0, 0) = 0.0;
	this->b(0, 0) = 0.0;
	this->r(0, 0) = 0.0;
	this->c10(0, 0) = 0.0;
	this->Q3min(0, 0) = 0.0;
	this->Q3max(0, 0) = 0.0;
	this->c50(0, 0) = 0.0;
	this->c30(0, 0) = 0.0;
	this->fq3(0, 0) = 0.0;
	this->Q30(0, 0) = 0.0;

	this->a1(0, 1) = 0.0;
	this->a5(0, 1) = 0.0;
	this->b(0, 1) = 0.0;
	this->r(0, 1) = 0.0;
	this->c10(0, 1) = 0.0;
	this->Q3min(0, 1) = 0.0;
	this->Q3max(0, 1) = 0.0;
	this->c50(0, 1) = 0.0;
	this->c30(0, 1) = 0.0;
	this->fq3(0, 1) = 0.0;
	this->Q30(0, 1) = 0.0;

	this->a1(0, 2) = 0.0;
	this->a5(0, 2) = 0.0;
	this->b(0, 2) = 0.0;
	this->r(0, 2) = 0.0;
	this->c10(0, 2) = 0.0;
	this->Q3min(0, 2) = 0.0;
	this->Q3max(0, 2) = 0.0;
	this->c50(0, 2) = 0.0;
	this->c30(0, 2) = 0.0;
	this->fq3(0, 2) = 0.0;
	this->Q30(0, 2) = 0.0;

	this->m(0) = 0.0;
	for (int j = 0; j < 3; j++)
	{
		for (int jj = 1; jj < nYield; jj++)
		{
			this->m(jj) = 0.0;
			this->a1(jj, j) = 0.0;
			this->a5(jj, j) = 0.0;
			this->b(jj, j) = 0.0;
			this->r(jj, j) = 0.0;
			this->c10(jj, j) = 0.0;
			this->Q3min(jj, j) = 0.0;
			this->Q3max(jj, j) = 0.0;
			this->c50(jj, j) = 0.0;
			this->c30(jj, j) = 0.0;
			this->fq3(jj, j) = 0.0;
			this->Q30(jj, j) = 0.0;
		}
	}

	double a_11D = 0.0;
	double a_33D = 0.0;
	double a_55D = 0.0;
	double a_33U0 = 0.0;
	double a_11U0 = 0.0;
	double a_55U0 = 0.0;
	double a_33UC = 0.0;
	double b_33UC = 0.0;
	double a_11UC = 0.0;
	double b_11UC = 0.0;
	double a_55UC = 0.0;
	double b_55UC = 0.0;

	double chi3_cal = 0.;
	if (Q3_maxD_cal > Q3max(0, 0)) {
		chi3_cal = 0.0;
	}
	else {
		chi3_cal = 0.0;
	}
	
	Vector Q1_calD(nYield), q1_calD(nYield), Q3_calD(nYield), q3_calD(nYield), Q5_calD(nYield), q5_calD(nYield);

	Q1_calD(0) = 0.0;
	q1_calD(0) = 0.0;
	Q3_calD(0) = 0.0;
	q3_calD(0) = 0.0;
	Q5_calD(0) = 0.0;
	q5_calD(0) = 0.0;

	H_0D.Zero();
	H_0D(0, 0) = 0.0;
	H_0D(2, 2) = 0.0;
	H_0D(4, 4) = 0.0;

	Vector O_3U0(nYield);
	if (Q3_maxU0_cal > Q3max(0, 0)) {
		O_3U0(0) = 0.0;
		chi3_cal = 0.0;
	}
	else if (Q3_maxU0_cal < Q3max(0, 0) && Q3_maxU0_cal - Q3_U0u_cal > Q3max(0, 1)) {
		O_3U0(0) = 0.0;
		chi3_cal = 0.0;
	}
	else {
		chi3_cal = 0.0;
	}

	Vector Q1_calU0, q1_calU0, Q3_calU0, q3_calU0, Q5_calU0, q5_calU0;

	Q1_calU0(0) = 0.0;
	q1_calU0(0) = 0.0;
	Q3_calU0(0) = 0.0;
	q3_calU0(0) = 0.0;
	Q5_calU0(0) = 0.0;
	q5_calU0(0) = 0.0;

	H_0U0.Zero();
	H_0UC.Zero();

	this->H_r0D33 = 0.0;
	this->H_r0U033 = 0.0;
	this->H_r0UC33 = 0.0;

	Matrix H_pD(6 * nYield, 6), H_pU0(6 * nYield, 6), H_pUC(6 * nYield, 6);
	H_pD.Zero();
	H_pU0.Zero();
	H_pUC.Zero();

	for (int n = 1; n < nYield; n++)
	{
		if (Q3_maxD_cal > Q3max(n, 0)) {
			chi3_cal = 0.0;
		}
		else {
			chi3_cal = 0.0;
		}
		Q1_calD(n) = 0.0;
		q1_calD(n) = 0.0;
		Q3_calD(n) = 0.0;
		q3_calD(n) = 0.0;
		Q5_calD(n) = 0.0;
		q5_calD(n) = 0.0;

		Matrix H_eD_Q(6,6);
		H_eD_Q = H_0D;
		
		Matrix C_eD(6, 6);
		C_eD.Zero();

		chi3_cal = 0.0;

		Q1_calU0(n) = 0.0;
		q1_calU0(n) = 0.0;
		Q3_calU0(n) = 0.0;
		q3_calU0(n) = 0.0;
		Q5_calU0(n) = 0.0;
		q5_calU0(n) = 0.0;

		Matrix H_eU0_Q(6, 6);
		H_eU0_Q = H_0U0;
		
		Matrix C_eU0(6, 6);
		C_eU0.Zero();

	Mq_min = 0.0;
	Mq_max = 0.0;
	dC0 = 0.0;
	Vector a_C(3), b_C(2);
	a_C(0) = 0.;
	a_C(1) = 0.;
	a_C(2) = 0.;
	b_C(0) = 0.;
	b_C(1) = 0.;

	Tf_min = 0.;
	Tf_max = 0.;
	a_U = 0.;
	this->Tf0 = 0.;

	deltaT = 1.;

	Vector a1_ev(nYield), ar2_ev(nYield), b_ev(nYield), r_ev(nYield), fq3_ev(nYield), Q30_ev(nYield);

	for (int n = 0; n < nYield; n++)
	{
		a1_ev(n) = 0.;
		ar2_ev(n) = 0.;
		b_ev(n) = 0.;
		r_ev(n) = 0.;
		fq3_ev(n) = 0.;
		Q30_ev(n) = 0.;
	}
	}

	Matrix Hs_dummy(3,3), Cs(6, 6), Cs_dummy(3, 3);
	Hs.Zero();
	Hs_dummy.Zero();
	Cs.Zero();
	Cs_dummy.Zero();

	this->beta=0.;

	theTangent.Zero();

	PKt.Zero();

	this->TrialStressElasticFlag = 0;
	this->CommittedStressElasticFlag = 0;

	Sum_pf.Zero();

}


TIMUP3D::~TIMUP3D() {

	return;
};


int TIMUP3D::setTrialStrain(const Vector& pStrain) {

	Tstrain = pStrain;

	this->ComputeResponse();

	return 0;

};

int TIMUP3D::setTrialStrain(const Vector& v, const Vector& r) {

	return this->setTrialStrain(v);

};

int TIMUP3D::setTrialStrainIncr(const Vector& v) {

	Tstrain = Cstrain + v;

	this->ComputeResponse();
	return 0;

};

int TIMUP3D::getOrder() const
{
	return 3;
}

int TIMUP3D::setTrialStrainIncr(const Vector& v, const Vector& r) {

	return this->setTrialStrainIncr(v);


};


const Matrix& TIMUP3D::getTangent() {
	
	return this->theTangent;
};


const Matrix& TIMUP3D::getInitialTangent() {

	Matrix initTangent(6, 6);
	initTangent.Zero();
	initTangent(0, 0) = H11el_D;
	initTangent(2, 2) = H33el_D;
	initTangent(4, 4) = H55el_D;
	return initTangent;
};


const Vector& TIMUP3D::getStress(void) {

	return Tstress;

};

const Vector& TIMUP3D::getStrain(void) {

	return Tstrain;
};

const Vector& TIMUP3D::getCommittedStress(void) {

	return Cstress;
};

const Vector& TIMUP3D::getCommittedStrain(void) {

	return Cstrain;

};


int TIMUP3D::commitState(void) {

	double TstrainNorm_c = pow(pow(Tstrain(0), 2.) + pow(Tstrain(1), 2.) + pow(Tstrain(2), 2.) + pow(Tstrain(3), 2.) + pow(Tstrain(4), 2.) + pow(Tstrain(5), 2.), .5);
	double CstrainNorm_c = pow(pow(Cstrain(0), 2.) + pow(Cstrain(1), 2.) + pow(Cstrain(2), 2.) + pow(Cstrain(3), 2.) + pow(Cstrain(4), 2.) + pow(Cstrain(5), 2.), .5);

	if (abs(CstrainNorm_c - TstrainNorm_c) > 1.e-14) {
		Cq_p = q_p;
	}

	Cstress = Tstress;
	Cstrain = Tstrain;

	CnFlows = NumYielded;
	
	CO = O_conv;

	CCsum_q_p = Csum_q_p;
	Csum_q_p = sum_q_p;

	CHs = Hs;
	Cpf = pf;
	CSum_pf = Sum_pf;

	Ccs = cs;

	dCC = C - CC;
	CCC = CC;
	CC = C;

	CU = U;

	return 0;

};

int TIMUP3D::revertToLastCommit(void) {

	Tstress = Cstress;

	Tstrain = Cstrain;

	cs = Ccs;

	TrialStressElasticFlag = CommittedStressElasticFlag;

	return 0;
};

int TIMUP3D::revertToStart(void) {
	Tstress.Zero();

	Tstrain.Zero();

	cs.Zero();
	
	TrialStressElasticFlag = 0;

	return 0;
}



NDMaterial* TIMUP3D::getCopy(void) {

	TIMUP3D* theSA = new TIMUP3D(*this);
	return theSA;

};



NDMaterial* TIMUP3D::getCopy(const char* code)
{
	if (strcmp(code, "ThreeDimensional") == 0) {
		TIMUP3D* theSA = new TIMUP3D(*this);
		return theSA;
	}

	return 0;
}


int TIMUP3D::sendSelf(int commitTag, Channel& theChannel) {
	
	int res = 0;
	
	static Vector data(135 + 24 * nYield + 9 + 1);

	data(0) = this->getTag();
	data(1) = nYield;
	data(2) = niter;
	data(3) = H11el_D;
	data(4) = H33el_D;
	data(5) = H55el_D;
	data(6) = mHyperEl;
	data(7) = H11el_U0;
	data(8) = H33el_U0;
	data(9) = H55el_U0;
	data(10) = H11el_UC;
	data(11) = H33el_UC;
	data(12) = H55el_UC;
    data(13) = Q3max_D;
	data(14) = Q3min_D;
	data(15) = Q3max_U0;
	data(16) = Q3min_U0;
	data(17) = Q3max_UC;
	data(18) = Q3min_UC;
	data(19) = a1_DN;
	data(20) = a5_DN;
	data(21) = b_DN;
	data(22) = r_DN;
	data(23) = Tf_min;
	data(24) = Tf_max;
	data(25) = a_U;
	data(26) = Mq_min;
	data(27) = Mq_max;
	data(28) = dC0;
	data(29) = tol;

	data(30) = b_11D;
	data(31) = b_33D;
	data(32) = b_55D;
	data(33) = b_11U0;
	data(34) = b_33U0;
	data(35) = b_55U0;

	data(36) = Q3_maxD_cal;
	data(37) = Q3_maxU0_cal;
	data(38) = Q3_maxUC_cal;
	data(39) = Q3_Du_cal;
	data(40) = Q3_U0u_cal;
	data(41) = Q3_UCu_cal;

	int idata = 41;
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			idata = idata + 1;
			data(idata) = theTangent(i, j);
		}
	}
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			idata = idata + 1;
			data(idata) = PKt(i, j);
		}
	}
	
	idata = idata + 1;
	data(idata) = a1_U0N;
	idata = idata + 1;
	data(idata) = a5_U0N;
	idata = idata + 1;
	data(idata) = b_U0N;
	idata = idata + 1;
	data(idata) = r_U0N;

	idata = idata + 1;
	data(idata) = a1_UCN;
	idata = idata + 1;
	data(idata) = a5_UCN;
	idata = idata + 1;
	data(idata) = b_UCN;
	idata = idata + 1;
	data(idata) = r_UCN;

	idata = idata + 1;
	data(idata) = Tstress(0);
	idata = idata + 1;
	data(idata) = Tstress(1);
	idata = idata + 1;
	data(idata) = Tstress(2);
	idata = idata + 1;
	data(idata) = Tstress(3);
	idata = idata + 1;
	data(idata) = Tstress(4);
	idata = idata + 1;
	data(idata) = Tstress(5);

	idata = idata + 1;
	data(idata) = Tstrain(0);
	idata = idata + 1;
	data(idata) = Tstrain(1);
	idata = idata + 1;
	data(idata) = Tstrain(2);
	idata = idata + 1;
	data(idata) = Tstrain(3);
	idata = idata + 1;
	data(idata) = Tstrain(4);
	idata = idata + 1;
	data(idata) = Tstrain(5);

	idata = idata + 1;
	data(idata) = Cstress(0);
	idata = idata + 1;
	data(idata) = Cstress(1);
	idata = idata + 1;
	data(idata) = Cstress(2);
	idata = idata + 1;
	data(idata) = Cstress(3);
	idata = idata + 1;
	data(idata) = Cstress(4);
	idata = idata + 1;
	data(idata) = Cstress(5);

	idata = idata + 1;
	data(idata) = Cstrain(0);
	idata = idata + 1;
	data(idata) = Cstrain(1);
	idata = idata + 1;
	data(idata) = Cstrain(2);
	idata = idata + 1;
	data(idata) = Cstrain(3);
	idata = idata + 1;
	data(idata) = Cstrain(4);
	idata = idata + 1;
	data(idata) = Cstrain(5);

	idata = idata + 1;
	data(idata) = fq3_DN;
	idata = idata + 1;
	data(idata) = fq3_U0N;
	idata = idata + 1;
	data(idata) = fq3_UCN;

	idata = idata + 1;
	data(idata) = Q30_DN;
	idata = idata + 1;
	data(idata) = Q30_U0N;
	idata = idata + 1;
	data(idata) = Q30_UCN;

	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
		idata = idata + 1;
		data(idata) = a1(i,ii);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			idata = idata + 1;
			data(idata) = a5(i, ii);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			idata = idata + 1;
			data(idata) = b(i, ii);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			idata = idata + 1;
			data(idata) = r(i, ii);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			idata = idata + 1;
			data(idata) = c10(i, ii);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			idata = idata + 1;
			data(idata) = Q3min(i, ii);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			idata = idata + 1;
			data(idata) = Q3max(i, ii);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			idata = idata + 1;
			data(idata) = c50(i, ii);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			idata = idata + 1;
			data(idata) = c30(i, ii);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			idata = idata + 1;
			data(idata) = fq3(i, ii);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			idata = idata + 1;
			data(idata) = Q30(i, ii);
		}
	}

	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			idata = idata + 1;
			data(idata) = H_0D(i, j);
		}
	}
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			idata = idata + 1;
			data(idata) = H_0U0(i, j);
		}
	}
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			idata = idata + 1;
			data(idata) = H_0UC(i, j);
		}
	}

	idata = idata + 1;
	data(idata) = H_r0D33;
	idata = idata + 1;
	data(idata) = H_r0U033;
	idata = idata + 1;
	data(idata) = H_r0UC33;

	for (int i = 0; i < 6 + nYield; i++) {
		for (int j = 0; j < 6; j++) {
			idata = idata + 1;
			data(idata) = H_eD(i, j);
		}
	}
	for (int i = 0; i < 6 + nYield; i++) {
		for (int j = 0; j < 6; j++) {
			idata = idata + 1;
			data(idata) = H_eU0(i, j);
		}
	}
	for (int i = 0; i < 6 + nYield; i++) {
		for (int j = 0; j < 6; j++) {
			idata = idata + 1;
			data(idata) = H_eUC(i, j);
		}
	}

	idata = idata + 1;
	data(idata) = Mq_min;
	idata = idata + 1;
	data(idata) = Mq_max;
	idata = idata + 1;
	data(idata) = dC0;

	idata = idata + 1;
	data(idata) = Tf_min;
	idata = idata + 1;
	data(idata) = Tf_max;
	idata = idata + 1;
	idata = a_U;
	idata = idata + 1;
	idata = Tf0;

	idata = idata + 1;
	data(idata) = perm;

	idata = idata + 1;
	data(idata) = deltaT;

	res = theChannel.sendVector(this->getDbTag(), commitTag, data);

	return 0;
};

int TIMUP3D::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) {
	
	static Vector data(135 + 24 * nYield + 9 + 1);
	int res = theChannel.recvVector(0, commitTag, data);
	this->setTag(data(0));
	nYield = data(1);
	niter = data(2);
	H11el_D = data(3);
	H33el_D = data(4);
	H55el_D = data(5);
	mHyperEl = data(6);
	H11el_U0 = data(7);
	H33el_U0 = data(8);
	H55el_U0 = data(9);
	H11el_UC = data(10);
	H33el_UC = data(11);
	H55el_UC = data(12);
    Q3max_D = data(13);
	Q3min_D = data(14);
	Q3max_U0 = data(15);
	Q3min_U0 = data(16);
	Q3max_UC = data(17);
	Q3min_UC = data(18);
	a1_DN = data(19);
	a5_DN = data(20);
	b_DN = data(21);
	r_DN = data(22);
	Tf_min = data(23);
	Tf_max = data(24);
	a_U = data(25);
	Mq_min = data(26);
	Mq_max = data(27);
	dC0 = data(28);
	tol = data(29);

    b_11D = data(30);
	b_33D = data(31);
	b_55D = data(32);
	b_11U0 = data(33);
	b_33U0 = data(34);
	b_55U0 = data(35);

	Q3_maxD_cal = data(36);
	Q3_maxU0_cal = data(37);
	Q3_maxUC_cal = data(38);
	Q3_Du_cal = data(39);
	Q3_U0u_cal = data(40);
	Q3_UCu_cal = data(41);

	int rdata = 41;
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			rdata = rdata + 1;
			theTangent(i, j) = data(rdata);
		}
	}
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			rdata = rdata + 1;
			PKt(i, j) = data(rdata);
		}
	}

	rdata = rdata + 1;
	a1_U0N = data(rdata);
	rdata = rdata + 1;
	a5_U0N = data(rdata);
	rdata = rdata + 1;
	b_U0N = data(rdata);
	rdata = rdata + 1;
	r_U0N = data(rdata);

	rdata = rdata + 1;
	a1_UCN = data(rdata);
	rdata = rdata + 1;
	a5_UCN = data(rdata);
	rdata = rdata + 1;
	b_UCN = data(rdata);
	rdata = rdata + 1;
	r_UCN = data(rdata);

	rdata = rdata + 1;
	Tstress(0) = data(rdata);
	rdata = rdata + 1;
	Tstress(1) = data(rdata);
	rdata = rdata + 1;
	Tstress(2) = data(rdata);
	rdata = rdata + 1;
	Tstress(3) = data(rdata);
	rdata = rdata + 1;
	Tstress(4) = data(rdata);
	rdata = rdata + 1;
	Tstress(5) = data(rdata);

	rdata = rdata + 1;
	Tstrain(0) = data(rdata);
	rdata = rdata + 1;
	Tstrain(1) = data(rdata);
	rdata = rdata + 1;
	Tstrain(2) = data(rdata);
	rdata = rdata + 1;
	Tstrain(3) = data(rdata);
	rdata = rdata + 1;
	Tstrain(4) = data(rdata);
	rdata = rdata + 1;
	Tstrain(5) = data(rdata);

	rdata = rdata + 1;
	Cstress(0) = data(rdata);
	rdata = rdata + 1;
	Cstress(1) = data(rdata);
	rdata = rdata + 1;
	Cstress(2) = data(rdata);
	rdata = rdata + 1;
	Cstress(3) = data(rdata);
	rdata = rdata + 1;
	Cstress(4) = data(rdata);
	rdata = rdata + 1;
	Cstress(5) = data(rdata);

	rdata = rdata + 1;
	Cstrain(0) = data(rdata);
	rdata = rdata + 1;
	Cstrain(1) = data(rdata);
	rdata = rdata + 1;
	Cstrain(2) = data(rdata);
	rdata = rdata + 1;
	Cstrain(3) = data(rdata);
	rdata = rdata + 1;
	Cstrain(4) = data(rdata);
	rdata = rdata + 1;
	Cstrain(5) = data(rdata);

	rdata = rdata + 1;
	fq3_DN = data(rdata);
	rdata = rdata + 1;
	fq3_U0N = data(rdata);
	rdata = rdata + 1;
	fq3_UCN = data(rdata);

	rdata = rdata + 1;
	Q30_DN = data(rdata);
	rdata = rdata + 1;
	Q30_U0N = data(rdata);
	rdata = rdata + 1;
	Q30_UCN = data(rdata);

	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			rdata = rdata + 1;
			a1(i, ii) = data(rdata);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			rdata = rdata + 1;
			a5(i, ii) = data(rdata);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			rdata = rdata + 1;
			b(i, ii) = data(rdata);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			rdata = rdata + 1;
			r(i, ii) = data(rdata);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			rdata = rdata + 1;
			c10(i, ii) = data(rdata);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			rdata = rdata + 1;
			Q3min(i, ii) = data(rdata);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			rdata = rdata + 1;
			Q3max(i, ii) = data(rdata);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			rdata = rdata + 1;
			c50(i, ii) = data(rdata);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			rdata = rdata + 1;
			c30(i, ii) = data(rdata);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			rdata = rdata + 1;
			fq3(i, ii) = data(rdata);
		}
	}
	for (int i = 0; i < nYield; i++) {
		for (int ii = 0; ii < 3; i++) {
			rdata = rdata + 1;
			Q30(i, ii) = data(rdata);
		}
	}

	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			rdata = rdata + 1;
			H_0D(i, j) = data(rdata);
		}
	}
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			rdata = rdata + 1;
			H_0U0(i, j) = data(rdata);
		}
	}
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < 6; j++) {
			rdata = rdata + 1;
			H_0UC(i, j) = data(rdata);
		}
	}

	rdata = rdata + 1;
	H_r0D33 = data(rdata);
	rdata = rdata + 1;
	H_r0U033 = data(rdata);
	rdata = rdata + 1;
	H_r0UC33 = data(rdata);

	for (int i = 0; i < 6 + nYield; i++) {
		for (int j = 0; j < 6; j++) {
			rdata = rdata + 1;
			H_eD(i, j) = data(rdata);
		}
	}
	for (int i = 0; i < 6 + nYield; i++) {
		for (int j = 0; j < 6; j++) {
			rdata = rdata + 1;
			H_eU0(i, j) = data(rdata);
		}
	}
	for (int i = 0; i < 6 + nYield; i++) {
		for (int j = 0; j < 6; j++) {
			rdata = rdata + 1;
			H_eUC(i, j) = data(rdata);
		}
	}

	rdata = rdata + 1;
	Mq_min = data(rdata);
	rdata = rdata + 1;
	Mq_max = data(rdata);
	rdata = rdata + 1;
	dC0 = data(rdata);

	rdata = rdata + 1;
	Tf_min = data(rdata);
	rdata = rdata + 1;
	Tf_max = data(rdata);
	rdata = rdata + 1;
	a_U = rdata;
	rdata = rdata + 1;
	Tf0 = rdata;

	rdata = rdata + 1;
	perm = data(rdata);

	rdata = rdata + 1;
	deltaT = data(rdata);

	return 0;
};

Response*
TIMUP3D::setResponse(const char** argv, int argc, OPS_Stream& s)
{
	if (strcmp(argv[0], "stress") == 0 || strcmp(argv[0], "stresses") == 0)
		return new MaterialResponse(this, 1, this->getStress());

	else if (strcmp(argv[0], "strain") == 0 || strcmp(argv[0], "strains") == 0)
		return new MaterialResponse(this, 2, this->getStrain());

	else if (strcmp(argv[0], "tangent") == 0 || strcmp(argv[0], "Tangent") == 0)
		return new MaterialResponse(this, 3, this->getTangent());
	else
		return 0;
}

int TIMUP3D::getResponse(int responseID, Information& matInfo) {

	switch (responseID) {
	case -1:
		return -1;
	case 1:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStress();
		return 0;

	case 2:
		if (matInfo.theVector != 0)
			*(matInfo.theVector) = getStrain();
		return 0;

	case 3:
		if (matInfo.theMatrix != 0)
			*(matInfo.theMatrix) = getTangent();
		return 0;
	}

	return 0;
};

void TIMUP3D::Print(OPS_Stream& s, int flag) {
	s << "TIMSoilAbutment3D tag: " << this->getTag() << endln;
	s << "  nYield: " << nYield << endln;
	s << "  H11el: " << H11el << endln;
	s << "  H33el: " << H33el << endln;
	s << "  H55el: " << H55el << endln;
	s << " Yield Centers 1:" << Cc1 << endln;
	s << " Yield Centers 3:" << Cc3 << endln;
	s << " Yield Centers 5:" << Cc5 << endln;
	s << " Active Yield Surfaces:" << CnFlows << endln;
	s << "  stress: " << Tstress << " tangent: " << theTangent << endln;
	return;
};

int TIMUP3D::setParameter(const char** argv, int argc, Parameter& param)
{
	if (argc < 1)
		return -1;

	if (strcmp(argv[0], "deltaT") == 0) {
		param.setValue(deltaT);
		return param.addObject(1, this);
	}

	return -1;
}

int TIMUP3D::updateParameter(int parameterID, Information& info)
{
	switch (parameterID) {

	case 1:
		deltaT = info.theDouble;
		return 0;

	default:
		return -1;
	}
}

double
DoubleDot2_2_Contr(const Vector& v1, const Vector& v2)
// computes doubledot product for vector-vector arguments, both "contravariant"
{
	if ((v1.Size() != 6) || (v2.Size() != 6))
		opserr << "\n ERROR! DoubleDot2_2_Contr requires vector of size(6)!" << endln;

	double result = 0.0;
	for (int i = 0; i < v1.Size(); i++) {
		result += v1(i) * v2(i) + (i>2) * v1(i) * v2(i);
	}

	return result;
}

int TIMUP3D::ComputeResponse() {

	double n_pf, n_i, Mq, Tf;

	Vector Q(6), q_p_int(6 * nYield), q_0(6), q_0_int(6), Q_int(6 * nYield), O(6 * nYield), f_n(nYield), f_n_cl(nYield), f_n_int(nYield), sum_q_p(6), inc_q_p(6), D(6 * nYield), lam(nYield), QU(6), QC(6), dy_dU(nYield),
		dy_dC(nYield), a1_ev(nYield), a5_ev(nYield), b_ev(nYield), r_ev(nYield), fq3_ev(nYield), Q30_ev(nYield);

	Vector RelVelN(6);
	RelVelN.Zero();

	Matrix Ct(6, 6);
	Ct.Zero();
	Matrix Kt(6, 6);
	Kt.Zero();

	double t1 = 0.0;
	double t2 = 0.0;
	double t3 = 0.0;
	double t4 = 0.0;
	double t5 = 0.0;
	double t6 = 0.0;

	Vector q(6);
	q = Tstrain;
	q(2) = -Tstrain(2);
	q(3) = Tstrain(3)*Lfound;
	q(4) = Tstrain(4)*Lfound;
	q(5) = Tstrain(5)*Lfound;

	for (int j = 0; j < 6; j++) {
		RelVelN(j) = (Tstrain(j) - Cstrain(j)) / deltaT;
	}

	Vector CQ;
	CQ = Cstress;
	CQ(2) = - CQ(2);
	CQ(3) = CQ(3) / Lfound;
	CQ(4) = CQ(4) / Lfound;
	CQ(5) = CQ(5) / Lfound;
	Vector Cdef;
	Cdef = Cstrain;
	Cdef(2) = -Cdef(2);
	Cdef(3) = Cdef(3) * Lfound;
	Cdef(4) = Cdef(4) * Lfound;
	Cdef(5) = Cdef(5) * Lfound;

	// velocity
	Vector DqDt(6);
	for (int i = 0; i < 6; i++) {
		DqDt(i) = RelVelN(i);
	}
	double DqDt_norm;
	DqDt_norm = pow(pow(DqDt(0),2.0)+ pow(DqDt(1), 2.0)+ pow(DqDt(2), 2.0)+ pow(DqDt(3), 2.0)+ pow(DqDt(4), 2.0)+ pow(DqDt(5), 2.0),0.5);

	// drainage coefficent
	Tf = DqDt_norm / perm;

	if (Tf <= Tf_min) {
		U = 1.;
	}
	else if ((Tf>Tf_min) && (Tf<Tf_max)) {
		U = 1. / (a_U*(Tf - Tf0));
	}
	else if (Tf >= Tf_max) {
		U = 0.;
	}

	// initial elastic - plastic stiffness matrix of the system
	Matrix Hs_dummy(3, 3), Cs(6, 6), Cs_dummy(3, 3);

	Hs.Zero();
	Hs_dummy.Zero();
	Cs.Zero();
	Cs_dummy.Zero();

	Hs = H_0D*U + (CC*H_0UC + (1. - CC)*H_0U0)*(1. - U);
	Hs_dummy(0, 0) = Hs(0, 0);
	Hs_dummy(1, 1) = Hs(2, 2);
	Hs_dummy(2, 2) = Hs(4, 4);
	// initial elastic-plastic compliant matrix
	Cs_dummy(0, 0) = 1. / Hs_dummy(0, 0);
	Cs_dummy(1, 1) = 1. / Hs_dummy(1, 1);
	Cs_dummy(2, 2) = 1. / Hs_dummy(2, 2);

	Cs(0, 0) = Cs_dummy(0, 0);
	Cs(2, 2) = Cs_dummy(1, 1);
	Cs(4, 4) = Cs_dummy(2, 2);

	// initial origin of the yield surfaces
	for (int n = 1; n < nYield+1; n++) {
		O(6 * n - 3 - 1) = (Q3max(n - 1, 0)*U + (CC*Q3max(n - 1, 2) + (1 -CC)*Q3max(n - 1, 1))*(1 - U) + Q3min(n - 1, 0)*U + (CC*Q3min(n - 1, 2) + (1 - CC)*Q3min(n - 1, 1))*(1 - U)) / 2.;
    }

	Matrix H_a(6,6);
	Vector Q_i(6), cs_a(6*nYield), O_a(6 * nYield), sum_q_p_a(6), QUa(6), QCa(6), Q_N(6), 
		pf_a(nYield), pf_b(nYield);

	double n_pf_b;
	n_pf_b = 0;													// number of active plastic flows at the previous iteration
	QUa.Zero();													// force variation due to a change of U inizialization
	QCa.Zero();				                                    // force variation due to a change of C inizialization
	Q_N.Zero();													// ultimate surface load value
	
	Q_i = CQ;													// current load inizialization
	O_a = CO;													// inizialization of the centers of the yield surfaces
	sum_q_p_a = Csum_q_p;										// inizialization of the total plastic displacement
	inc_q_p = Csum_q_p - CCsum_q_p;								// increment of plastic displacement
	H_a = CHs;													// inizialization of the current stiffness
	
	pf_a.Zero();
	cs_a = Ccs;													// surfaces center evolution inizialization
	Sum_pf = CSum_pf + pf_a;
	
	double TstrainNorm_init0 = pow(pow(Tstrain(0), 2.) + pow(Tstrain(1), 2.) + pow(Tstrain(2), 2.) + pow(Tstrain(3), 2.) + pow(Tstrain(4), 2.) + pow(Tstrain(5), 2.), .5);
	double CstrainNorm_init0 = pow(pow(Cstrain(0), 2.) + pow(Cstrain(1), 2.) + pow(Cstrain(2), 2.) + pow(Cstrain(3), 2.) + pow(Cstrain(4), 2.) + pow(Cstrain(5), 2.), .5);

	// reverse to last commit
	if (abs(CstrainNorm_init0 - TstrainNorm_init0) < 1.e-14) {
		
		q_p = Cq_p;

		Tstress = Cstress;
		Tstrain = Cstrain;

		NumYielded = CnFlows;

		O_conv = CO;

		sum_q_p = Csum_q_p;

		Hs = CHs;
		pf = Cpf;
		Sum_pf = CSum_pf;

		cs = Ccs;
		
		C = CC;

		U = CU;
	}
	else if (abs(CstrainNorm_init0 - TstrainNorm_init0) > 1.e-14) {

	// ----------------------------------------------------------
	pf_b.Zero();
	// update surfaces coefficents according to the current values for U and C
	for (int n = 1; n < nYield+1; n++) {
		a1_ev(n - 1) = a1(n - 1, 0)*U + (a1(n - 1, 1)*(1 - CC) + a1(n - 1, 2)*CC)*(1 - U);
		a5_ev(n - 1) = a5(n - 1, 0)*U + (a5(n - 1, 1)*(1 - CC) + a5(n - 1, 2)*CC)*(1 - U);
		b_ev(n - 1) = b(n - 1, 0)*U + (b(n - 1, 1)*(1 - CC) + b(n - 1, 2)*CC)*(1 - U);
		r_ev(n - 1) = r(n - 1, 0)*U + (r(n - 1, 1)*(1 - CC) + r(n - 1, 2)*CC)*(1 - U);
		fq3_ev(n - 1) = fq3(n - 1, 0)*U + (fq3(n - 1, 1)*(1 - CC) + fq3(n - 1, 2)*CC)*(1 - U);
		Q30_ev(n - 1) = Q30(n - 1, 0)*U + (Q30(n - 1, 1)*(1 - CC) + Q30(n - 1, 2)*CC)*(1 - U);
	}

    if (abs(U - CU) < tol) {
	   U = CU;
	}

	double TstrainNorm_init = pow(pow(Tstrain(0), 2.) + pow(Tstrain(1), 2.) + pow(Tstrain(2), 2.) + pow(Tstrain(3), 2.) + pow(Tstrain(4), 2.) + pow(Tstrain(5), 2.), .5);
	double CstrainNorm_init = pow(pow(Cstrain(0), 2.) + pow(Cstrain(1), 2.) + pow(Cstrain(2), 2.) + pow(Cstrain(3), 2.) + pow(Cstrain(4), 2.) + pow(Cstrain(5), 2.), .5);
		
	if (abs(U - CU) > 0 || (dCC > 0)) {

    	// update of the currently non-attained yield surfaces
	   	for (int n = 1; n < nYield; n++) {
			if (Sum_pf(n - 1) == 0) {
				O_a(6 * n - 3 - 1) = (2.*fq3_ev(n - 1)*r_ev(n - 1)) / 2.;
			}
			if (Sum_pf(n - 1) > 0) {
				O_a(6 * n - 3 - 1) = CO(6 * n - 3 - 1);
			}
	   	}

		double checkQi = O_a(2) + (2.*fq3_ev(0)*r_ev(0)) / 2.;

		// update of the first surface center; check and update of the back force
		if (Q_i(2) > checkQi) {
			O_a(2) = Q_i(2) - (2.*fq3_ev(0)*r_ev(0)) / 2.;
			cs_a(2) = O_a(2) - (2.*fq3_ev(0)*r_ev(0)) / 2.;
		}
		else if (Q_i(2) < O_a(2) - (2.*fq3_ev(0)*r_ev(0)) / 2.) {
			O_a(2) = Q_i(2) + (2.*fq3_ev(0)*r_ev(0)) / 2.;
			cs_a(2) = O_a(2) - (2.*fq3_ev(0)*r_ev(0)) / 2.;
		}
		else {
			cs_a(2) = O_a(2) - 2.*fq3_ev(0)*r_ev(0) / 2.;
		}

		// ultimate surface update
		O_a(6 * nYield - 3 - 1) = 2. * fq3_ev(nYield-1)*r_ev(nYield - 1) / 2.;
		cs_a(6 * nYield - 3 - 1) = 0.;

		// update of the centers of the inner yield surfaces
		for (int n = 1; n < nYield-1; n++) {
			if (O_a(6 * (n+1) - 3 - 1) + (2 * fq3_ev(n)*r_ev(n)) / 2. < O_a(6 * n - 3 - 1) + (2. * fq3_ev(n - 1)*r_ev(n - 1)) / 2.) {
				O_a(6 * (n + 1) - 3 - 1) = O_a(6 * n - 3 - 1) + (2. * fq3_ev(n - 1)*r_ev(n - 1)) / 2. - (2. * fq3_ev(n)*r_ev(n)) / 2.;
			}
			else if (O_a(6 * (n + 1) - 3 - 1) - (2. * fq3_ev(n)*r_ev(n)) / 2 > O_a(6 * n - 3 - 1) - (2. * fq3_ev(n - 1)*r_ev(n - 1)) / 2.) {
				O_a(6 * (n + 1) - 3 - 1) = O_a(6 * n - 3 - 1) - (2 * fq3_ev(n - 1)*r_ev(n - 1)) / 2. + (2. * fq3_ev(n)*r_ev(n)) / 2.;
			}
			else {
				if (Q_i(2) > O_a(6 * (n + 1) - 3 - 1) + (2. * fq3_ev(n)*r_ev(n)) / 2.) {
					O_a(6 * (n + 1) - 3 - 1) = Q_i(2) - (2. * fq3_ev(n)*r_ev(n)) / 2.;
				}
				else if (Q_i(2) < O_a(6 * (n + 1) - 3 - 1) - (2. * fq3_ev(n)*r_ev(n)) / 2.) {
					O_a(6 * (n + 1) - 3 - 1) = Q_i(2) + (2. * fq3_ev(n)*r_ev(n)) / 2.;
				}
			}
			cs_a(6 * (n + 1) - 3 - 1) = O_a(6 * (n + 1) - 3 - 1) - (2. * fq3_ev(n)*r_ev(n)) / 2.;
		}

		// check the attainment of the ultimate surface
		double f_n_check;
		Vector v_pcN(6);
		v_pcN.Zero();
		f_n_check = (pow((Q_i(0) - cs_a(6 * nYield - 5 - 1)) / a1_ev(nYield-1), 2.) + pow((Q_i(4) - cs_a(6 * nYield - 1 - 1)) / a5_ev(nYield - 1),2.)) * pow(2. * b_ev(nYield - 1) - (Q_i(2) - cs_a(6 * nYield - 3 - 1) - Q30_ev(nYield - 1)) / fq3_ev(nYield - 1),2.) + pow((Q_i(2) - cs_a(6 * nYield - 3 - 1) - Q30_ev(nYield-1)) / fq3_ev(nYield - 1) - b_ev(nYield - 1), 2.) - pow(r_ev(nYield - 1), 2.);
		if (f_n_check > -tol) {
			Q_N = CQ;
			for (int j = 0; j < 6; j++) {
				v_pcN(j) = (Q_i(j) - O_a(6*(nYield-1)-5-1+j)) / 2.;
			}
			while  ((f_n_check > tol) || (f_n_check < -tol)) {
				v_pcN = v_pcN / 2.;
				if (f_n_check < -tol) {
					Q_N = Q_N + v_pcN;
				}
				else if (f_n_check > tol) {
					Q_N = Q_N - v_pcN;
				}
				f_n_check = (pow(Q_N(0) / a1_ev(nYield - 1),2.) + pow(Q_N(4) / a5_ev(nYield - 1),2.))*pow(2 * b_ev(nYield - 1) - (Q_N(2) - Q30_ev(nYield - 1)) / fq3_ev(nYield - 1),2.) + pow((Q_N(2) - Q30_ev(nYield - 1)) / fq3_ev(nYield - 1) - b_ev(nYield - 1),2.) - pow(r_ev(nYield - 1),2.);
			}
		}
		
		// check the attainment of the inner yield surfaces
		for (int n = 1; n < nYield; n++) {
			f_n_check = (pow((Q_i(0) - cs_a(6 * n - 5 - 1)) / a1_ev(n - 1),2.) + pow((Q_i(4) - cs_a(6 * n - 1 - 1)) / a5_ev(n - 1),2.))*pow(2. * b_ev(n - 1) - (Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1),2.) + pow((Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1) - b_ev(n - 1),2.) - pow(r_ev(n - 1),2.);
			if (f_n_check > -tol) {
		   		pf_a(n - 1) = 1;
			}
	    	else {
		   		pf_a(n - 1) = 0;
			}
		}

		Vector tag_pf_a(nYield);
		tag_pf_a.Zero();
		for (int j = 0; j < nYield; j++) {
			if (pf_a(j) == 1) {
				tag_pf_a(j) = j;
			}
		}
		double Max_tag_pf_a = tag_pf_a(0);
		for (int j = 1; j < nYield; j++) {
			if (tag_pf_a(j) > tag_pf_a(j-1)) {
				Max_tag_pf_a = tag_pf_a(j);
			}
		}
	
		if (Max_tag_pf_a >= 1) {
			for (int j = 0; j < Max_tag_pf_a; j++) {
				pf_a(j) = 1;
			}
		}

		// check the correct position of the surfaces
		Vector f_n_cl(nYield), v_cl(6);
		f_n_cl.Zero();
		v_cl.Zero();
		for (int n = 1; n < nYield; n++) {
			if (pf_a(n - 1) == 1) {
				// inizialization
				f_n_cl(n - 1) = (pow((Q_i(0) - cs_a(6 * n - 5 - 1)) / a1_ev(n - 1),2.) + pow((Q_i(4) - cs_a(6 * n - 1 - 1)) / a5_ev(n - 1),2.))*pow(2. * b_ev(n - 1) - (Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1),2.) + pow((Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1) - b_ev(n - 1),2.) - pow(r_ev(n - 1),2.);
	    		if (f_n_cl(n - 1) > tol) {
		   			// center - point load versor
					for (int jj = 0; jj < 6; jj++) {
		   				v_cl(jj) = (Q_i(jj) - O_a(6 * n - 5 - 1 + jj)) / pow(pow(Q_i(0) - O_a(6 * n - 5 - 1),2.) + pow(Q_i(2) - O_a(6 * n - 3 - 1),2.) + pow(Q_i(4) - O_a(6 * n - 1 - 1),2.),0.5);
					}
					while (f_n_cl(n - 1) < -tol || f_n_cl(n - 1) > tol) {
		      			if (f_n_cl(n - 1) < -tol) {
			     			v_cl = v_cl / 2.;
							for (int jj = 0; jj < 6; jj++) {
	             				cs_a(6 * n - 5 - 1 + jj) = cs_a(6 * n - 5 - 1 + jj) - v_cl(jj);
							}
			  			}
	          			else if (f_n_cl(n - 1) > tol) {
							for (int jj = 0; jj < 6; jj++) {
								cs_a(6 * n - 5 - 1 + jj) = cs_a(6 * n - 5 - 1 + jj) + v_cl(jj);
							}
						}
		      			f_n_cl(n - 1) = (pow((Q_i(0) - cs_a(6 * n - 5 - 1)) / a1_ev(n - 1),2.) + pow((Q_i(4) - cs_a(6 * n - 1 - 1)) / a5_ev(n - 1),2.))*pow(2. * b_ev(n - 1) - (Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1),2.) + pow((Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1) - b_ev(n - 1),2.) - pow(r_ev(n - 1),2.);
					}
				}
			}
		}

		// ultimate surface check
		if (pf_a(nYield-1) == 1) {
			Q_i = Q_N;
			for (int n = 1; n < nYield; n++) {
				// inizialization
				f_n_cl(n - 1) = (pow((Q_i(0) - cs_a(6 * n - 5 - 1)) / a1_ev(n - 1),2.) + pow((Q_i(4) - cs_a(6 * n - 1 - 1)) / a5_ev(n - 1),2.))*pow(2. * b_ev(n - 1) - (Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1),2.) + pow((Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1) - b_ev(n - 1),2.) - pow(r_ev(n - 1),2.);
				if (f_n_cl(n - 1) > tol) {
					for (int jj = 0; jj < 6; jj++) {
						cs_a(6 * n - 5 - 1 + jj) = cs_a(6 * n - 5 - 1 + jj) + Q_i(jj) - CQ(jj);
					}
				}
			}
		}
		
		// check yielded surfaces
		for (int n = 1; n < nYield; n++) {
			Vector Doubledot_H_Dstrain(6);
			Doubledot_H_Dstrain.Zero();
			for (int j = 0; j < 6; j++) {
				for (int jj = 0; jj < 6; jj++) {
					Doubledot_H_Dstrain(j) = Doubledot_H_Dstrain(j) + H_a(j,jj) * (q(jj) - Cdef(jj));
				}
			}
			f_n_check = (pow((Q_i(0) + Doubledot_H_Dstrain(0) - cs_a(6 * n - 5 - 1)) / a1_ev(n - 1), 2.) + pow((Q_i(4) + Doubledot_H_Dstrain(4) - cs_a(6 * n - 1 - 1)) / a5_ev(n - 1), 2.))*pow(2. * b_ev(n - 1) - (Q_i(2) + Doubledot_H_Dstrain(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1), 2.) + pow((Q_i(2) + Doubledot_H_Dstrain(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1) - b_ev(n - 1), 2.) - pow(r_ev(n - 1), 2.);
			if (f_n_check > -tol) {
				pf_a(n - 1) = 1;
			}
			else {
				pf_a(n - 1) = 0;
			}
		}
	
		// initial stiffness matrix
		if (Sum_pf(0) == 0) {
			H_0D(2 , 2) = H33el_D;
			H_0U0(2 , 2) = H33el_U0;
			H_0UC(2 , 2) = H33el_UC;
		}	
		else {
			H_0D(2 , 2) = ((pow((2. * mHyperEl - 1),2.0)) / (2.* mHyperEl))*(pow((2. * H_r0D33),(mHyperEl / (2. * mHyperEl - 1))))*pow(Q_i(2),((2. * mHyperEl - 2.) / (2. * mHyperEl - 1)));
			H_0U0(2 , 2) = ((pow((2. * mHyperEl - 1),2.0)) / (2. * mHyperEl))*(pow((2. * H_r0U033),(mHyperEl / (2. * mHyperEl - 1))))*pow(Q_i(2),((2. * mHyperEl - 2.) / (2. * mHyperEl - 1)));
			H_0UC(2 , 2) = ((pow((2. * mHyperEl - 1),2.0)) / (2. * mHyperEl))*(pow((2. * H_r0UC33),(mHyperEl / (2. * mHyperEl - 1))))*pow(Q_i(2),((2. * mHyperEl - 2.) / (2. * mHyperEl - 1)));
		}
		Matrix H0(6, 6);
		H0 = H_0D*U + (H_0U0*(1 - CC) + H_0UC*CC)*(1 - U);

		// current compliance
		Matrix H0_dummy(3,3), C_a(6, 6), C_a_dummy(3, 3), H_a_dummy(3, 3);
		H0_dummy.Zero();
		H0_dummy(0, 0) = H0(0, 0);
		H0_dummy(1, 1) = H0(2, 2);

		H0_dummy(2, 2) = H0(4, 4);
		H0_dummy.Invert(C_a_dummy);
		C_a.Zero();
		C_a(0, 0) = C_a_dummy(0, 0);
		C_a(2, 2) = C_a_dummy(1, 1);
		C_a(4, 4) = C_a_dummy(2, 2);

		Vector D(6*nYield), D_Hen_D(nYield);
		Matrix dC_a(6, 6);
		dC_a.Zero();
		D_Hen_D.Zero();

		for (int n = 1; n < nYield; n++) {
			if (pf_a(n-1) == 1) {
				// coefficents
				D(6 * n - 5 - 1) = 2. * ((Q_i(0) - cs_a(6 * n - 5 - 1)) / pow(a1_ev(n-1),2.))*pow(2. * b_ev(n-1) - (Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n-1)) / fq3_ev(n-1),2.);
				D(6 * n - 1 - 1) = 2. * ((Q_i(4) - cs_a(6 * n - 1 - 1)) / pow(a5_ev(n - 1), 2.))*pow(2. * b_ev(n - 1) - (Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1), 2.);
				D(6 * n - 3 - 1) = -2. / fq3_ev(n-1)*((pow((Q_i(0) - cs_a(6 * n - 5 - 1)) / a1_ev(n-1),2.) + pow((Q_i(4) - cs_a(6 * n - 1 - 1)) / a5_ev(n - 1), 2.))*(2. * b_ev(n-1) - (Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n-1)) / fq3_ev(n-1)) - ((Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n-1)) / fq3_ev(n-1) - b_ev(n-1)));
				Vector Doubledot_D_Hen(6);
				Doubledot_D_Hen.Zero();
				for (int j = 0; j < 6; j++) {
					for (int jj = 0; jj < 6; jj++) {
						Doubledot_D_Hen(j) = Doubledot_D_Hen(j) + (H_eD(6 * n - 5 - 1 + jj,j)*U + (H_eUC(6 * n - 5 - 1 + jj, j)*CC + H_eU0(6 * n - 5 - 1 + jj, j)*(1. - CC))*(1. - U)) * D(6 * n - 5 - 1 + jj);
					}
				}

				for (int jj = 0; jj < 6; jj++) {
					D_Hen_D(n - 1) = D_Hen_D(n - 1) + Doubledot_D_Hen(jj) * D(6 * n - 5 - 1 + jj);
				}

				Vector SignCQ(6), Signinc_q_p(6);

				for (int j = 0; j < 6; j++) {
					SignCQ(j) = CQ(j) / abs(CQ(j));
					if (abs(CQ(j)) == 0) {
						SignCQ(j) = 0.;
					}
				}
				for (int j = 0; j < 6; j++) {
					Signinc_q_p(j) = inc_q_p(j) / abs(inc_q_p(j));
					if (abs(inc_q_p(j)) == 0) {
						Signinc_q_p(j) = 0.;
					}
				}

				dC_a(0, 0) = (1. + beta*CC*SignCQ(0)*Signinc_q_p(0))*D(6 * n - 5 - 1)*D(6 * n - 5 - 1) / D_Hen_D(n - 1);
				dC_a(0, 4) = (1. + beta*CC*SignCQ(4)*Signinc_q_p(4))*D(6 * n - 5 - 1)*D(6 * n - 1 - 1) / D_Hen_D(n - 1);
				dC_a(0, 2) = (1. + beta*CC*SignCQ(2)*Signinc_q_p(2))*D(6 * n - 5 - 1)*D(6 * n - 3 - 1) / D_Hen_D(n - 1);
				dC_a(4, 0) = (1. + beta*CC*SignCQ(0)*Signinc_q_p(0))*D(6 * n - 1 - 1)*D(6 * n - 5 - 1) / D_Hen_D(n - 1);
				dC_a(4, 4) = (1. + beta*CC*SignCQ(4)*Signinc_q_p(4))*D(6 * n - 1 - 1)*D(6 * n - 1 - 1) / D_Hen_D(n - 1);
				dC_a(4, 2) = (1. + beta*CC*SignCQ(2)*Signinc_q_p(2))*D(6 * n - 1 - 1)*D(6 * n - 3 - 1) / D_Hen_D(n - 1);
				dC_a(2, 0) = (1. + beta*CC*SignCQ(0)*Signinc_q_p(0))*D(6 * n - 3 - 1)*D(6 * n - 5 - 1) / D_Hen_D(n - 1);
				dC_a(2, 4) = (1. + beta*CC*SignCQ(4)*Signinc_q_p(4))*D(6 * n - 3 - 1)*D(6 * n - 1 - 1) / D_Hen_D(n - 1);
				dC_a(2, 2) = (1. + beta*CC*SignCQ(2)*Signinc_q_p(2))*D(6 * n - 3 - 1)*D(6 * n - 3 - 1) / D_Hen_D(n - 1);

				C_a = C_a + dC_a;
			}
		}		
										
		// current stiffness of the system inizialization
		C_a_dummy(0, 0) = C_a(0, 0);
		C_a_dummy(1, 1) = C_a(2, 2);
		C_a_dummy(2, 2) = C_a(4, 4);
		C_a_dummy(0, 1) = C_a(0, 2);
		C_a_dummy(0, 2) = C_a(0, 4);
		C_a_dummy(1, 0) = C_a(2, 0);
		C_a_dummy(1, 2) = C_a(2, 4);
		C_a_dummy(2, 0) = C_a(4, 0);
		C_a_dummy(2, 1) = C_a(4, 2);

		H_a_dummy.Zero();
		C_a_dummy.Invert(H_a_dummy);
		H_a(0, 0) = H_a_dummy(0, 0);
		H_a(2, 2) = H_a_dummy(1, 1);
		H_a(4, 4) = H_a_dummy(2, 2);
		H_a(0, 2) = H_a_dummy(0, 1);
		H_a(0, 4) = H_a_dummy(0, 2);
		H_a(2, 0) = H_a_dummy(1, 0);
		H_a(2, 4) = H_a_dummy(1, 2);
		H_a(4, 0) = H_a_dummy(2, 0);
		H_a(4, 2) = H_a_dummy(2, 1);

		Vector qUa(6), QUa(6);
		double da1_dU,da5_dU, db_dU, dr_dU, dfq3_dU,
		dQ30_dU, dy_dU1, dy_dU2, dy_dU3,
		dy_dU;
		qUa.Zero();
		QUa.Zero();

		// current force variation due to change in U
		for (int n = 1; n < nYield+1; n++) {
			if (pf_a(n - 1) == 1) {
			da1_dU = a1(n - 1, 0) - (a1(n - 1, 1)*(1 - CC) + a1(n - 1, 2)*CC);
			da5_dU = a5(n - 1, 0) - (a5(n - 1, 1)*(1 - CC) + a5(n - 1, 2)*CC);
			db_dU = b(n - 1, 0) - (b(n - 1, 1)*(1 - CC) + b(n - 1, 2)*CC);
			dr_dU = r(n - 1, 0) - (r(n - 1, 1)*(1 - CC) + r(n - 1, 2)*CC);
			dfq3_dU = fq3(n - 1, 0) - (fq3(n - 1, 1)*(1 - CC) + fq3(n - 1, 2)*CC);
			dQ30_dU = Q30(n - 1, 0) - (Q30(n - 1, 1)*(1 - CC) + Q30(n - 1, 2)*CC);
			dy_dU1 = -2. * (da1_dU*(pow(Q_i(0) - cs_a(6 * n - 5 - 1),2.) / pow(a1_ev(n - 1),3.))+ da5_dU*(pow(Q_i(4) - cs_a(6 * n - 1 - 1),2.) / pow(a5_ev(n - 1),3.))) * pow(2. * b_ev(n - 1) - (Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1),2.);
			dy_dU2 = 2. * (pow((Q_i(0) - cs_a(6 * n - 5 - 1)) / a1_ev(n - 1),2.) + pow((Q_i(4) - cs_a(6 * n - 1 - 1)) / a5_ev(n - 1),2.)) * (2. * b_ev(n - 1) - (Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1)) * (2. * db_dU + (dQ30_dU*fq3_ev(n - 1) + dfq3_dU*(Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1))) / pow(fq3_ev(n - 1),2.));
			dy_dU3 = 2. * (((Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1)) - b_ev(n - 1))*((-dQ30_dU*fq3_ev(n - 1) - dfq3_dU*(Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1))) / pow(fq3_ev(n - 1),2.) - db_dU);
			dy_dU = dy_dU1 + dy_dU2 + dy_dU3 - 2. * r_ev(n - 1)*dr_dU;
			
			qUa(0) = qUa(0) + (dy_dU*(U - CU)*D(6 * n - 5 - 1)) / D_Hen_D(n - 1);
			qUa(2) = qUa(2) + (dy_dU*(U - CU)*D(6 * n - 3 - 1)) / D_Hen_D(n - 1);
			qUa(4) = qUa(4) + (dy_dU*(U - CU)*D(6 * n - 1 - 1)) / D_Hen_D(n - 1);
			}
		}

		for (int j = 0; j < 6; j++) {
			for (int jj = 0; jj < 6; jj++) {
				QUa(j) = QUa(j) + H_a(j, jj)*qUa(jj);
			}
		}

		Vector qCa(6), QCa(6);
		double da1_dC,da5_dC, db_dC, dr_dC, dfq3_dC,
		dQ30_dC, dy_dC1, dy_dC2, dy_dC3,
		dy_dC;
		qCa.Zero();
		QCa.Zero();

		// current force variation due to change in C
		for (int n = 1; n < nYield+1; n++) {
			if (pf_a(n - 1) == 1) {
				da1_dC = (a1(n - 1, 2) - a1(n - 1, 1)) * (1 - U);
				da5_dC = (a5(n - 1, 2) - a5(n - 1, 1)) * (1 - U);
				db_dC = (b(n - 1, 2) - b(n - 1, 1)) * (1 - U);
				dr_dC = (r(n - 1, 2) - r(n - 1, 1)) * (1 - U);
				dfq3_dC = (fq3(n - 1, 2) - fq3(n - 1, 1)) * (1 - U);
				dQ30_dC = (Q30(n - 1, 2) - Q30(n - 1, 1)) * (1 - U);
				dy_dC1 = -2. * (da1_dC*pow((Q_i(0) - cs_a(6 * n - 5 - 1)) / a1_ev(n - 1),3.) + da5_dC*pow((Q_i(4) - cs_a(6 * n - 1 - 1)) / a5_ev(n - 1),3.)) * pow(2. * b_ev(n - 1) - (Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1),2.);
				dy_dC2 = 2. * (pow((Q_i(0) - cs_a(6 * n - 5 - 1)) / a1_ev(n - 1),2.) + pow((Q_i(4) - cs_a(6 * n - 1 - 1)) / a5_ev(n - 1),2.)) * (2 * b_ev(n - 1) - (Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1))*(2 * db_dC + (dQ30_dC*fq3_ev(n - 1) + dfq3_dC*(Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1))) / pow(fq3_ev(n - 1),2.));
				dy_dC3 = 2. * ((Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1) - b_ev(n - 1))*((-dQ30_dC*fq3_ev(n - 1) - dfq3_dC*(Q_i(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1))) / pow(fq3_ev(n - 1),2.) - db_dC);
				dy_dC = dy_dC1 + dy_dC2 + dy_dC3 - 2. * r_ev(n - 1)*dr_dC;

				qCa(0) = qCa(0) + dy_dC*(CC - CCC)*D(6 * n - 5 - 1) / D_Hen_D(n - 1);
				qCa(2) = qCa(2) + dy_dC*(CC - CCC)*D(6 * n - 3 - 1) / D_Hen_D(n - 1);
				qCa(4) = qCa(4) + dy_dC*(CC - CCC)*D(6 * n - 1 - 1) / D_Hen_D(n - 1);

			}
		}

		for (int j = 0; j < 6; j++) {
			for (int jj = 0; jj < 6; jj++) {
				QCa(j) = QCa(j) + H_a(j, jj)*qCa(jj);
			}
		}
		
		// load update due to a variation of U and C
		Vector SignCQ_c(6), Signinc_q_p_c(6), SignQ_i_c(6);

		for (int j = 0; j < 6; j++) {
			SignCQ_c(j) = CQ(j) / abs(CQ(j));
			if (abs(CQ(j)) == 0) {
				SignCQ_c(j) = 0.;
			}
		}
		for (int j = 0; j < 6; j++) {
			Signinc_q_p_c(j) = inc_q_p(j) / abs(inc_q_p(j));
			if (abs(inc_q_p(j)) == 0) {
				Signinc_q_p_c(j) = 0.;
			}
		}
		for (int j = 0; j < 6; j++) {
			SignQ_i_c(j) = Q_i(j) / abs(Q_i(j));
			if (abs(Q_i(j)) == 0) {
				SignQ_i_c(j) = 0.;
			}
		}

			if (pf_a(nYield-1) == 0) {
				Q_i(0) = CQ(0) - (1 + beta*CC*SignQ_i_c(0)*Signinc_q_p_c(0))*(QUa(0) + QCa(0));
				Q_i(2) = CQ(2) - (1 + beta*CC*SignQ_i_c(2)*Signinc_q_p_c(2))*(QUa(2) + QCa(2));
				Q_i(4) = CQ(4) - (1 + beta*CC*SignQ_i_c(4)*Signinc_q_p_c(4))*(QUa(4) + QCa(4));

				// surface center update due to changes in U and/or C
				for (int n = 1; n < nYield; n++) { // 20260210 for (int n = 1; n < nYield-1; n++) {
				if (pf_a(n-1) == 1) {
					for (int j = 0; j < 6; j++) {
						cs_a(6 * n - 5 - 1 + j) = cs_a(6 * n - 5 - 1 + j) + (Q_i(j) - CQ(j));
					}
				}
			}
			}
	}

	Vector Q_b(6);                                           
	Q_b.Zero();
	Vector Q_a(6);                                            
	Q_a = Q_i;
	n_i = 0;
	
	double Norm_Qab = pow(pow(Q_a(0) - Q_b(0),2.) + pow(Q_a(1) - Q_b(1), 2.) + pow(Q_a(2) - Q_b(2), 2.) + pow(Q_a(3) - Q_b(3), 2.) + pow(Q_a(4) - Q_b(4), 2.) + pow(Q_a(5) - Q_b(5), 2.),.5);
	 
	double Norm_pfab = 0.;
	for (int j = 0; j < nYield-1; j++) {
		Norm_pfab += pow(pf_a(j) - pf_b(j), 2.);
	}
	Norm_pfab = pow(Norm_pfab, 0.5);

	while ((Norm_Qab > tol) && (n_i < niter) || (Norm_pfab > 0) && (n_i < niter) || (n_i == 0)) {
		n_i = n_i + 1;
		Q_b = Q_a;
		for (int j = 0; j < 6; j++) {
			Q_a(j) = Q_i(j);
			for (int jj = 0; jj < 6; jj++) {
				Q_a(j) = Q_a(j) + H_a(j, jj)*(q(jj) - Cdef(jj));
			}
		}
		double Norm_QN = 0.;
		for (int j = 0; j < 6; j++) {
			Norm_QN = Norm_QN + pow(Q_N(j), 2.);
		}
		Norm_QN = pow(Norm_QN,.5);
		// ultimate surface check
		double f_n_check = (pow(Q_a(0) / a1_ev(nYield - 1), 2.) + pow(Q_a(4) / a5_ev(nYield - 1), 2.))*pow(2. * b_ev(nYield - 1) - (Q_a(2) - Q30_ev(nYield - 1)) / fq3_ev(nYield - 1), 2.) + pow((Q_a(2) - Q30_ev(nYield - 1)) / fq3_ev(nYield - 1) - b_ev(nYield - 1), 2.) - pow(r_ev(nYield - 1), 2.);
		if ((f_n_check > -tol) && (Norm_QN == 0)) {
		// vector between the preeceding load and the last surface center
			Vector v_pcN(6);
			v_pcN.Zero();
			for (int j = 0; j < 6; j++) {
				for (int jj = 0; jj < 6; jj++) {
					v_pcN(j) = v_pcN(j) + 2. * H_a(j,jj)*(q(jj) - Cdef(jj));
				}
			}
		while  ((f_n_check > tol) || (f_n_check < -tol)) {
			v_pcN = v_pcN / 2.;
			if (f_n_check < -tol) {
				Q_a = Q_a + v_pcN;
			}
			else if (f_n_check > tol) {
				Q_a = Q_a - v_pcN;
			}
			f_n_check = (pow(Q_a(0) / a1_ev(nYield - 1), 2.) + pow(Q_a(4) / a5_ev(nYield - 1), 2.))*pow(2. * b_ev(nYield - 1) - (Q_a(2) - Q30_ev(nYield - 1)) / fq3_ev(nYield - 1), 2.) + pow((Q_a(2) - Q30_ev(nYield - 1)) / fq3_ev(nYield - 1) - b_ev(nYield - 1), 2.) - pow(r_ev(nYield - 1), 2.);
			}
		}
		else if ((f_n_check > -tol) && (Norm_QN > 0)) {
			Q_a = Q_N;
		}

		pf_b = pf_a;
		pf_a.Zero();

		if (Sum_pf(0) > 0) {
			H_0D(2, 2) = ((pow((2. * mHyperEl - 1), 2.0)) / (2.* mHyperEl))*(pow((2. * H_r0D33), (mHyperEl / (2. * mHyperEl - 1))))*pow(Q_a(2), ((2. * mHyperEl - 2.) / (2. * mHyperEl - 1)));
			H_0U0(2, 2) = ((pow((2. * mHyperEl - 1), 2.0)) / (2. * mHyperEl))*(pow((2. * H_r0U033), (mHyperEl / (2. * mHyperEl - 1))))*pow(Q_a(2), ((2. * mHyperEl - 2.) / (2. * mHyperEl - 1)));
			H_0UC(2, 2) = ((pow((2. * mHyperEl - 1), 2.0)) / (2. * mHyperEl))*(pow((2. * H_r0UC33), (mHyperEl / (2. * mHyperEl - 1))))*pow(Q_a(2), ((2. * mHyperEl - 2.) / (2. * mHyperEl - 1)));
		}
		Matrix H0(6, 6), H0_dummy(3, 3), C_a_dummy(3, 3), C_a(6, 6);
		H0.Zero();
		H0_dummy.Zero();
		C_a_dummy.Zero();
		C_a.Zero();
		H_a.Zero();
		H0 = H_0D*U + (H_0U0*(1 - CC) + H_0UC*CC)*(1 - U);

		H0_dummy(0, 0) = H0(0, 0);
		H0_dummy(1, 1) = H0(2, 2);
		H0_dummy(2, 2) = H0(4, 4);

		// current compliant of the system(inside the iteration)

		C_a_dummy(0, 0) = 1. / H0_dummy(0, 0);
		C_a_dummy(1, 1) = 1. / H0_dummy(1, 1);
		C_a_dummy(2, 2) = 1. / H0_dummy(2, 2);

		C_a(0, 0) = C_a_dummy(0, 0);
		C_a(2, 2) = C_a_dummy(1, 1);
		C_a(4, 4) = C_a_dummy(2, 2);

		Vector f_n(nYield);
		for (int n = 1; n < nYield + 1; n++) {
			// nth yield function
			f_n(n-1) = (pow((Q_a(0) - cs_a(6 * n - 5 - 1)) / a1_ev(n - 1), 2.) + pow((Q_a(4) - cs_a(6 * n - 1 - 1)) / a5_ev(n - 1), 2.))*pow(2. * b_ev(n - 1) - (Q_a(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1), 2.) + pow((Q_a(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1) - b_ev(n - 1), 2.) - pow(r_ev(n - 1), 2.);
			
			if (f_n(n - 1) >= -tol) {
				pf_a(n-1) = 1; // stop
				if (n-1 == 1) {
					int breakp = 1;
				}
				// load increment vector
				Vector vers(6);
				vers.Zero();
				vers = Q_a - Q_i;

				// inizialization
				Vector Q_int(6);
				Q_int = Q_i;
				double f_n_int = (pow((Q_int(0) - cs_a(6 * n - 5 - 1)) / a1_ev(n - 1), 2.) + pow((Q_int(4) - cs_a(6 * n - 1 - 1)) / a5_ev(n - 1), 2.))*pow(2. * b_ev(n - 1) - (Q_int(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1), 2.) + pow((Q_int(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1) - b_ev(n - 1), 2.) - pow(r_ev(n - 1), 2.);
				while (f_n_int < -tol) {
					Q_int = Q_int + vers / 100.;
					f_n_int = (pow((Q_int(0) - cs_a(6 * n - 5 - 1)) / a1_ev(n - 1), 2.) + pow((Q_int(4) - cs_a(6 * n - 1 - 1)) / a5_ev(n - 1), 2.))*pow(2. * b_ev(n - 1) - (Q_int(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1), 2.) + pow((Q_int(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1) - b_ev(n - 1), 2.) - pow(r_ev(n - 1), 2.);
				}
				// coefficents
				D(6 * n - 5 - 1) = 2. * ((Q_int(0) - cs_a(6 * n - 5 - 1)) / pow(a1_ev(n - 1), 2.))*pow(2. * b_ev(n - 1) - (Q_int(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1), 2.);
				D(6 * n - 1 - 1) = 2. * ((Q_int(4) - cs_a(6 * n - 1 - 1)) / pow(a5_ev(n - 1), 2.))*pow(2. * b_ev(n - 1) - (Q_int(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1), 2.);
				D(6 * n - 3 - 1) = -2. / fq3_ev(n - 1)*((pow((Q_int(0) - cs_a(6 * n - 5 - 1)) / a1_ev(n - 1), 2.) + pow((Q_int(4) - cs_a(6 * n - 1 - 1)) / a5_ev(n - 1), 2.))*(2. * b_ev(n - 1) - (Q_int(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1)) - ((Q_int(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1) - b_ev(n - 1)));
				
				// term to add in the compliant matrix of the system
				Vector Doubledot_D_Hen(6);
				Doubledot_D_Hen.Zero();

				for (int j = 0; j < 6; j++) {
					for (int jj = 0; jj < 6; jj++) {
						Doubledot_D_Hen(j) = Doubledot_D_Hen(j) + (H_eD(6 * n - 5 - 1 + jj, j)*U + (H_eUC(6 * n - 5 - 1 + jj, j)*CC + H_eU0(6 * n - 5 - 1 + jj, j)*(1. - CC))*(1. - U)) * D(6 * n - 5 - 1 + jj);
					}
				}

				double D_Hen_D = 0;
				for (int jj = 0; jj < 6; jj++) {
					double D2 = D(6 * n - 5 - 1 + jj);
					D_Hen_D = D_Hen_D + Doubledot_D_Hen(jj) * D(6 * n - 5 - 1 + jj);
				}

				Matrix dC_a(6, 6);
				Vector SignCQ(6), Signinc_q_p(6);
				
				for (int j = 0; j < 6; j++) {
					SignCQ(j) = CQ(j) / abs(CQ(j));
					if (abs(CQ(j))==0) {
						SignCQ(j) = 0.;
					}
				}
				for (int j = 0; j < 6; j++) {
					Signinc_q_p(j) = inc_q_p(j) / abs(inc_q_p(j));
					if (abs(inc_q_p(j)) == 0) {
						Signinc_q_p(j) = 0.;
					}
				}

				dC_a(0, 0) = (1. + beta*CC*SignCQ(0)*Signinc_q_p(0))*D(6 * n - 5 - 1)*D(6 * n - 5 - 1) / D_Hen_D;
				dC_a(0, 4) = (1. + beta*CC*SignCQ(4)*Signinc_q_p(4))*D(6 * n - 5 - 1)*D(6 * n - 1 - 1) / D_Hen_D;
				dC_a(0, 2) = (1. + beta*CC*SignCQ(2)*Signinc_q_p(2))*D(6 * n - 5 - 1)*D(6 * n - 3 - 1) / D_Hen_D;
				dC_a(4, 0) = (1. + beta*CC*SignCQ(0)*Signinc_q_p(0))*D(6 * n - 1 - 1)*D(6 * n - 5 - 1) / D_Hen_D;
				dC_a(4, 4) = (1. + beta*CC*SignCQ(4)*Signinc_q_p(4))*D(6 * n - 1 - 1)*D(6 * n - 1 - 1) / D_Hen_D;
				dC_a(4, 2) = (1. + beta*CC*SignCQ(2)*Signinc_q_p(2))*D(6 * n - 1 - 1)*D(6 * n - 3 - 1) / D_Hen_D;
				dC_a(2, 0) = (1. + beta*CC*SignCQ(0)*Signinc_q_p(0))*D(6 * n - 3 - 1)*D(6 * n - 5 - 1) / D_Hen_D;
				dC_a(2, 4) = (1. + beta*CC*SignCQ(4)*Signinc_q_p(4))*D(6 * n - 3 - 1)*D(6 * n - 1 - 1) / D_Hen_D;
				dC_a(2, 2) = (1. + beta*CC*SignCQ(2)*Signinc_q_p(2))*D(6 * n - 3 - 1)*D(6 * n - 3 - 1) / D_Hen_D;
				
				// modified compliant stiffness matrix
				C_a = C_a + dC_a;

				// plastic multiplyer
				double lam1 = D(6 * n - 5 - 1)*(Q_a(0) - Q_i(0)) / D_Hen_D;
				double lam2 = D(6 * n - 3 - 1)*(Q_a(2) - Q_i(2)) / D_Hen_D;
				double lam3 = D(6 * n - 1 - 1)*(Q_a(4) - Q_i(4)) / D_Hen_D;
				double lam = lam1 + lam2 + lam3;

				// plastic displacements
				inc_q_p(6 * n - 5 - 1) = D(6 * n - 5 - 1)*lam;
				inc_q_p(6 * n - 3 - 1) = D(6 * n - 3 - 1)*lam;
				inc_q_p(6 * n - 1 - 1) = D(6 * n - 1 - 1)*lam;
				q_p(6 * n - 5 - 1) = Cq_p(6 * n - 5 - 1) + inc_q_p(6 * n - 5 - 1);
				q_p(6 * n - 3 - 1) = Cq_p(6 * n - 3 - 1) + inc_q_p(6 * n - 3 - 1);
				q_p(6 * n - 1 - 1) = Cq_p(6 * n - 1 - 1) + inc_q_p(6 * n - 1 - 1);
			}
			if (f_n(n - 1) < -tol) {
				q_p(6 * n - 5 - 1) = Cq_p(6 * n - 5 - 1);
				q_p(6 * n - 3 - 1) = Cq_p(6 * n - 3 - 1);
				q_p(6 * n - 1 - 1) = Cq_p(6 * n - 1 - 1);
			}
		}

		for (int j = 0; j < 6; j++) {
			sum_q_p(j) = 0.;
		}

		// total plastic deformation
		for (int n = 1; n < nYield + 1; n++) {
			for (int j = 0; j < 6; j++) {
				sum_q_p(j) = sum_q_p(j) + q_p(6 * n - 5 - 1 + j);
			}
		}

		// current stiffness matrix of the system
		C_a_dummy(0, 0) = C_a(0, 0);
		C_a_dummy(1, 1) = C_a(2, 2);
		C_a_dummy(2, 2) = C_a(4, 4);
		C_a_dummy(0, 1) = C_a(0, 2);
		C_a_dummy(0, 2) = C_a(0, 4);
		C_a_dummy(1, 0) = C_a(2, 0);
		C_a_dummy(1, 2) = C_a(2, 4);
		C_a_dummy(2, 0) = C_a(4, 0);
		C_a_dummy(2, 1) = C_a(4, 2);

		Matrix H_a_dummy(3, 3);
		H_a_dummy.Zero();
		//H_a_dummy.Invert(C_a_dummy);
		C_a_dummy.Invert(H_a_dummy);
		
		H_a(0, 0) = H_a_dummy(0, 0);
		H_a(2, 2) = H_a_dummy(1, 1);
		H_a(4, 4) = H_a_dummy(2, 2);
		H_a(0, 2) = H_a_dummy(0, 1);
		H_a(0, 4) = H_a_dummy(0, 2);
		H_a(2, 0) = H_a_dummy(1, 0);
		H_a(2, 4) = H_a_dummy(1, 2);
		H_a(4, 0) = H_a_dummy(2, 0);
		H_a(4, 2) = H_a_dummy(2, 1);

		Norm_Qab = pow(pow(Q_a(0) - Q_b(0), 2.) + pow(Q_a(1) - Q_b(1), 2.) + pow(Q_a(2) - Q_b(2), 2.) + pow(Q_a(3) - Q_b(3), 2.) + pow(Q_a(4) - Q_b(4), 2.) + pow(Q_a(5) - Q_b(5), 2.), .5);

		}

		//  degradation
		Mq = pow(pow(sum_q_p(0), 2.) + pow(sum_q_p(1), 2.) + pow(sum_q_p(2), 2.) + pow(sum_q_p(3), 2.) + pow(sum_q_p(4), 2.) + pow(sum_q_p(5), 2.), 0.5);
		if (Mq <= CMq) {
			Mq = CMq;
		}

		double TstrainNorm = pow(pow(Tstrain(0), 2.) + pow(Tstrain(1), 2.) + pow(Tstrain(2), 2.) + pow(Tstrain(3), 2.) + pow(Tstrain(4), 2.) + pow(Tstrain(5), 2.), .5);
		double CstrainNorm = pow(pow(Cstrain(0), 2.) + pow(Cstrain(1), 2.) + pow(Cstrain(2), 2.) + pow(Cstrain(3), 2.) + pow(Cstrain(4), 2.) + pow(Cstrain(5), 2.), .5);

		dC = 0.;
		if ((Mq > Mq_min) && (Mq <= Mq_max)) {
			dC = (a_C(0) + a_C(1) * Mq + a_C(2) * pow(Mq, 2.)) * (1. - b_C(0) * CC - b_C(1) * pow(CC, 2.));
		}
		if ((Mq > Mq_max) && CC<1) {
			dC = (a_C(0) + a_C(1) * Mq_max + a_C(2) * pow(Mq_max, 2.)) * (1. - b_C(0) * CC - b_C(1) * pow(CC, 2.));
		}
		if (CC >= 1) {
			dC = 0.;
			CC = 1;
		}
		if (U >= 0.99) {
			dC = 0.;
		}
		
		if (abs(CstrainNorm - TstrainNorm) < 1.e-14) {
			C = CC;
		}
		else if (abs(CstrainNorm - TstrainNorm) > 1.e-14) {
			C = CC + dC;
		}
		if (C >= 1) {
			C = 1.;
		}

	double f_n_cl_check = 0.;
	for (int n = 1; n < nYield; n++) {
		if ((pf_a(n - 1) == 1) && (pf_a(nYield - 1) == 0)) {
			for (int j = 0; j < 6; j++) {
				for (int jj = 0; jj < 6; jj++) {
					cs_a(6 * n - 5 - 1 + j) = cs_a(6 * n - 5 - 1 + j) + (H_eD(6 * n - 5 - 1 + j,jj) * U + (H_eUC(6 * n - 5 - 1 + j,jj) * CC + H_eU0(6 * n - 5 - 1 + j,jj) * (1-CC)) * (1-U))*(inc_q_p(6 * n - 5 - 1 + jj));
				}
			}
		}
		else if ((pf_a(n - 1) == 1) && (pf_a(nYield - 1) == 1)) {
			f_n_cl_check = (pow((Q_a(0) - cs_a(6 * n - 5 - 1)) / a1_ev(n - 1), 2.) + pow((Q_a(4) - cs_a(6 * n - 1 - 1)) / a5_ev(n - 1), 2.))*pow(2. * b_ev(n - 1) - (Q_a(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1), 2.) + pow((Q_a(2) - cs_a(6 * n - 3 - 1) - Q30_ev(n - 1)) / fq3_ev(n - 1) - b_ev(n - 1), 2.) - pow(r_ev(n - 1), 2.);
			if ((f_n_cl_check < -tol) || (f_n_cl_check > tol)) {
				for (int j = 0; j < 6; j++) {
				cs_a(6 * n - 5 - 1 + j) = cs_a(6 * n - 5 - 1 + j)+Q_a(j)-Q_i(j);
				}
			}
		}

		O_conv(6 * n - 5 - 1) = cs_a(6 * n - 5 - 1);
		O_conv(6 * n - 3 - 1) = (2*fq3_ev(n - 1)*r_ev(n - 1)) / 2.+cs_a(6 * n - 3 - 1);
		O_conv(6 * n - 1 - 1) = cs_a(6 * n - 1 - 1);
	}

	O_conv(6 * nYield - 3 - 1) = (2. * fq3_ev(nYield - 1) * r_ev(nYield - 1)) / 2.;

	Tstress(0) = Q_a(0);
	Tstress(1) = Q_a(1);
	Tstress(2) = - Q_a(2);
	Tstress(3) = Q_a(3)*Lfound;
	Tstress(4) = Q_a(4)*Lfound;
	Tstress(5) = Q_a(5)*Lfound;

	dCC = dC;

	cs = cs_a;

	if (abs(CstrainNorm - TstrainNorm) < 1.e-14) {
		Hs = CHs;
		theTangent = Hs;
		pf = Cpf;
	}
	else if (abs(CstrainNorm - TstrainNorm) > 1.e-14) {
		Hs = H_a;
		theTangent = H_a;
		pf = pf_a;
	}

	NumYielded = 0;
	for (int k = 0; k < nYield - 1; k++) {
		if (pf_a(k) == 1) {
		NumYielded = NumYielded + 1;
		}
	}

	sum_q_p = sum_q_p_a;

	Sum_pf = Sum_pf + pf_a;

	}

	return 0;
};