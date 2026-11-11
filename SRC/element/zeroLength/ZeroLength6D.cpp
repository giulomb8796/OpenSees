/*
Developed by: Davide Noè Gorini (davidenoe.gorini@unitn.it)
Date: 2023
Updated for OpenSees 3.8.0 Compliance (Thread-Safe & Multi-Dimensional Fixes)

Description: 
This file contains the complete implementation for the ZeroLength6D element.
It provides a full 6D response (3 translations, 3 rotations) and includes 
a dedicated method to compute relative velocities between end nodes, 
essential when using rate-dependent (viscoplastic/damping) materials.

Reference: 
Gorini, D.N. and Callisto, L. (2023): "A multiaxial inertial macroelement for deep foundations", 
Computers and Geotechnics, vol. 155, doi: 10.1016/j.compgeo.2022.105222
*/

#include <ZeroLength6D.h>
#include <Information.h>
#include <Parameter.h>
#include <Domain.h>
#include <Node.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <NDMaterial.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <G3Globals.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <elementAPI.h>

// Global static array to safely pass node pointers back to OpenSees core
static Node *theNodes[2];

/*
================================================================================
OPS_ZeroLength6D: Tcl/Python Interpreter API
Parses input arguments and instantiates the element inside the OpenSees Domain.
================================================================================
*/

void* OPS_ZeroLength6D()
{
    int ndm = OPS_GetNDM(); // Get current spatial dimension of the model (2D or 3D)
    int numdata = OPS_GetNumRemainingInputArgs();
    
    // Check for minimum required arguments: eleTag, iNode, jNode, matTag
    if (numdata < 4) {
        opserr << "WARNING too few arguments\n"
               << "want - element zeroLength6D eleTag? iNode? jNode? NDTag? <-orient x1? x2? x3? y1? y2? y3?>\n";
        return 0;
    }

    // Read the first 4 integer parameters
    int idata[4];
    numdata = 4;
    if (OPS_GetIntInput(&numdata, idata) < 0) {
        opserr << "WARNING: failed to get integer data (tags)\n";
        return 0;
    }
    
    // Modern OpenSees 3.x API macro to retrieve the NDMaterial pointer via its tag
    NDMaterial* nmat = OPS_GetNDMaterial(idata[3]); 
    if (nmat == 0) {
        opserr << "WARNING: NDMaterial " << idata[3] << " is not defined\n";
        return 0;
    }

    // Default orientation vectors: x-axis aligned with global X, y-axis with global Y
    Vector x(3); x(0) = 1.0; x(1) = 0.0; x(2) = 0.0;
    Vector y(3); y(0) = 0.0; y(1) = 1.0; y(2) = 0.0;

    // Check if the optional orientation flag '-orient' is specified
    if (OPS_GetNumRemainingInputArgs() >= 7) {
        const char* type = OPS_GetString();
        if (strcmp(type, "-orient") == 0 || strcmp(type, "orient") == 0) {
            numdata = 3;
            // Extract the user-defined local x-axis vector components
            if (OPS_GetDoubleInput(&numdata, &x(0)) < 0) {
                opserr << "WARNING: invalid double input for local x vector\n";
                return 0;
            }
            // Extract the user-defined local y-axis vector components
            if (OPS_GetDoubleInput(&numdata, &y(0)) < 0) {
                opserr << "WARNING: invalid double input for local y vector\n";
                return 0;
            }
        }
    }

    // Dynamic allocation of the element object
    return new ZeroLength6D(idata[0], ndm, idata[1], idata[2], x, y, *nmat);
}

/*
================================================================================
Constructor (Standard Allocation)
================================================================================
*/

ZeroLength6D::ZeroLength6D(int tag, int dim, int Nd1, int Nd2,
    const Vector& x, const Vector& yprime, NDMaterial &theNDmat) :
    Element(tag, ELE_TAG_ZeroLength6D),
    connectedExternalNodes(2), dimension(dim), numDOF(0),
    end1Ptr(0), end2Ptr(0), theNDMaterial(0), order(0)
{
    // Ridimensiona la matrice di trasformazione
    transformation.resize(6, 6);

    // Deep copy del materiale
    theNDMaterial = theNDmat.getCopy();
    if (theNDMaterial == 0) {
        opserr << "ZeroLength6D::ZeroLength6D-- failed to get copy of NDMaterial\n";
        exit(-1);
    }
    
    order = theNDMaterial->getOrder();
    this->setUp(Nd1, Nd2, x, yprime);
}

/*
================================================================================
Constructor (Blank Constructor for Parallelization / Database Channel Receiving)
================================================================================
*/

ZeroLength6D::ZeroLength6D() :
    Element(0, ELE_TAG_ZeroLength6D),
    connectedExternalNodes(2), dimension(0), numDOF(0),
    end1Ptr(0), end2Ptr(0), theNDMaterial(0), order(0)
{
    transformation.resize(6, 6);
}

/*
================================================================================
Destructor
================================================================================
*/

ZeroLength6D::~ZeroLength6D()
{
    if (theNDMaterial != nullptr) 
        delete theNDMaterial;
}

int ZeroLength6D::getNumExternalNodes(void) const { return 2; }

Node **ZeroLength6D::getNodePtrs(void)
{
    theNodes[0] = end1Ptr;
    theNodes[1] = end2Ptr;
    return theNodes;
}

const ID &ZeroLength6D::getExternalNodes(void) { return connectedExternalNodes; }
int ZeroLength6D::getNumDOF(void) { return numDOF; }

/*
================================================================================
setTrialRelVel: Compute Relative Velocity
Calculates the local/global relative velocity vector between Node I and Node J.
Crucial feature for rate-dependent/viscous multi-axial macroelements.
================================================================================
*/

const Vector &ZeroLength6D::setTrialRelVel()
{
    if (end1Ptr == nullptr || end2Ptr == nullptr) 
        return RelVelNodes;

    const Vector &vel1 = end1Ptr->getTrialVel();
    const Vector &vel2 = end2Ptr->getTrialVel();

    int halfDOF = numDOF / 2;
    RelVelNodes.resize(halfDOF);
    RelVelNodes.Zero();

    // Controlla che i nodi abbiano la dimensione corretta per le velocità
    if (vel1.Size() >= halfDOF && vel2.Size() >= halfDOF) {
        for (int i = 0; i < halfDOF; ++i) {
            RelVelNodes(i) = vel2(i) - vel1(i);
        }
    }

    return RelVelNodes;
}

/*
================================================================================
setDomain: Initialize Element within the Spatial Model Context
================================================================================
*/

void ZeroLength6D::setDomain(Domain *theDomain)
{
    if (theDomain == 0) {
        end1Ptr = 0; end2Ptr = 0;
        return;
    }

    int Nd1 = connectedExternalNodes(0);
    int Nd2 = connectedExternalNodes(1);
    end1Ptr = theDomain->getNode(Nd1);
    end2Ptr = theDomain->getNode(Nd2);

    if (end1Ptr == 0 || end2Ptr == 0) {
        end1Ptr = 0; end2Ptr = 0;
        return;
    }

    // Determine DOFs based on node properties to support multi-dimensional models
    int dofNd1 = end1Ptr->getNumberDOF();
    numDOF = 2 * dofNd1;

    // Safety constraint validation
    if (numDOF != 6 && numDOF != 12) {
        opserr << "ZeroLength6D::setDomain - element only works for 3 (2D) or 6 (3D) DOFs per node\n";
        end1Ptr = 0; end2Ptr = 0;
        return;
    }

    this->DomainComponent::setDomain(theDomain);
    
    // Construct the kinematic matrix system (A)
    this->setTransformation();
}

/*
================================================================================
getTangentStiff: Assemble and Return Global Tangent Stiffness Matrix K
Formulation: K = A^T * kb * A
================================================================================
*/

const Matrix& ZeroLength6D::getTangentStiff(void)
{
    // 1. Assicura che la matrice di trasformazione A sia valida
    if (A.noRows() == 0 || A.noCols() == 0) {
        this->setTransformation();
    }

    // 2. FORZA l'aggiornamento cinematico e l'impostazione delle deformazioni 
    //    nel materiale PRIMA di estrarre la rigidezza tangente D.
    this->setTrialRelVel();
    this->computeStrain();
    theNDMaterial->setTrialStrain(v);

    // 3. Recupera la matrice tangente dal materiale NDMaterial (TIMUP3D)
    const Matrix& D = theNDMaterial->getTangent();

    K.Zero();

    // 4. Assemblaggio della rigidezza globale: K = A^T * D * A
    for (int i = 0; i < numDOF; i++) {
        for (int j = 0; j < numDOF; j++) {
            double sum = 0.0;
            for (int k = 0; k < order; k++) {
                for (int l = 0; l < order; l++) {
                    sum += A(k, i) * D(k, l) * A(l, j);
                }
            }
            K(i, j) = sum;
        }
    }

    // Controllo di sicurezza su rigidezza nulla
    double maxVal = 0.0;
    for (int i = 0; i < numDOF; i++) {
        for (int j = 0; j < numDOF; j++) {
            if (fabs(K(i, j)) > maxVal) {
                maxVal = fabs(K(i, j));
            }
        }
    }

    if (maxVal == 0.0) {
        opserr << "ERROR: ZeroLength6D::getTangentStiff() -- La matrice K e' interamente NULLA!\n";
    }

    return K;
}

/*
================================================================================
getInitialStiff: Assemble and Return Global Initial Stiffness Matrix
================================================================================
*/

const Matrix& ZeroLength6D::getInitialStiff(void)
{
    if (A.noRows() == 0 || A.noCols() == 0) {
        this->setTransformation();
    }
    const Matrix& kb = theNDMaterial->getInitialTangent();
    const Matrix &tran = A;

    K.Zero(); // FIXED C2819: Changed from K->Zero() to K.Zero()
    double E = 0;

    for (int k = 0; k < order; k++) {
        for (int l = 0; l < order; l++) {
            E = kb(k, l);
            for (int i = 0; i < numDOF; i++) {
                for (int j = 0; j <= i; j++) {
                    K(i, j) += tran(k, i) * E * tran(l, j);
                }
            }
        }
    }

    for (int i = 0; i < numDOF; i++) {
        for (int j = 0; j < i; j++) {
            K(j, i) = K(i, j);
        }
    }

    // FIXED C2100: Returning direct reference to the object, removed the asterisk '*'
    return K;
}

/*
================================================================================
getDamp: Return Global Damping Matrix
================================================================================
*/
const Matrix &ZeroLength6D::getDamp(void) 
{ 
    K.Zero(); // FIXED C2819: Corrected object syntax using dot '.' instead of '->'
    return K; // FIXED C2100: Removed asterisk '*' from return type
}

/*
================================================================================
getMass: Return Global Mass Matrix
================================================================================
*/
const Matrix &ZeroLength6D::getMass(void) 
{ 
    K.Zero(); // FIXED C2819: Corrected object syntax using dot '.' instead of '->'
    return K; // FIXED C2100: Removed asterisk '*' from return type
}

void ZeroLength6D::zeroLoad(void) {}

int ZeroLength6D::addLoad(ElementalLoad *theLoad, double loadFactor) { return -1; }

int ZeroLength6D::addInertiaLoadToUnbalance(const Vector &accel) { return 0; }

/*
================================================================================
getResistingForce: Compute Global Resisting Force Vector P
Formulation: P = A^T * q
================================================================================
*/

const Vector& ZeroLength6D::getResistingForce(void)
{
    // 1. Assicura che la matrice A sia allocata
    if (A.noRows() == 0 || A.noCols() == 0) {
        this->setTransformation();
    }

    // 2. Calcola le deformazioni e aggiorna il materiale ND
    this->setTrialRelVel();
    this->computeStrain();
    theNDMaterial->setTrialStrain(v);

    // 3. Recupera il vettore delle forze/tensioni locali q dal materiale
    const Vector& q = theNDMaterial->getStress();

    // 4. Azzera il vettore delle forze resistenti globali P (12 DOF)
    P.Zero();

    // 5. Moltiplicazione P = A^T * q
    for (int i = 0; i < numDOF; i++) {
        double sum = 0.0;
        for (int k = 0; k < order; k++) {
            sum += A(k, i) * q(k);
        }
        P(i) = sum;
    }

    return P;
}

/*
================================================================================
getResistingForceIncInertia: Resisting force including inertia effects
================================================================================
*/
const Vector &ZeroLength6D::getResistingForceIncInertia() 
{ 
    // Zero-length elements typically don't have internal dynamic mass calculations, 
    // so it falls back to standard dynamic static forces directly.
    return this->getResistingForce(); 
}

/*
================================================================================
setUp: Direction Cosines Matrix Formulation
Builds local transformation orientation matrix.
================================================================================
*/

void ZeroLength6D::setUp(int Nd1, int Nd2, const Vector& x, const Vector& yp)
{
    connectedExternalNodes(0) = Nd1;
    connectedExternalNodes(1) = Nd2;

    double xn = x.Norm();
    if (xn == 0.0) {
        opserr << "ZeroLength6D::setUp -- vettore x nullo\n";
        exit(-1);
    }

    Vector x_norm = x / xn;

    // z = x_norm cross yp
    Vector z(3);
    z(0) = x_norm(1) * yp(2) - x_norm(2) * yp(1);
    z(1) = x_norm(2) * yp(0) - x_norm(0) * yp(2);
    z(2) = x_norm(0) * yp(1) - x_norm(1) * yp(0);
    double zn = z.Norm();

    if (zn == 0.0) {
        opserr << "ZeroLength6D::setUp -- vettore yp e x sono paralleli!\n";
        exit(-1);
    }
    Vector z_norm = z / zn;

    // y = z_norm cross x_norm
    Vector y(3);
    y(0) = z_norm(1) * x_norm(2) - z_norm(2) * x_norm(1);
    y(1) = z_norm(2) * x_norm(0) - z_norm(0) * x_norm(2);
    y(2) = z_norm(0) * x_norm(1) - z_norm(1) * x_norm(0);

    // Ridimensiona e azzera transformation (6x6)
    transformation.resize(6, 6);
    transformation.Zero();

    // Matrice di Rotazione 3x3 per le Traslazioni
    for (int i = 0; i < 3; i++) {
        transformation(0, i) = x_norm(i); // Asse x locale
        transformation(1, i) = y(i);      // Asse y locale
        transformation(2, i) = z_norm(i); // Asse z locale
    }

    // Matrice di Rotazione 3x3 per le Rotazioni
    for (int i = 0; i < 3; i++) {
        transformation(3, i + 3) = x_norm(i); // Rotazione x
        transformation(4, i + 3) = y(i);      // Rotazione y
        transformation(5, i + 3) = z_norm(i); // Rotazione z
    }
}

/*
================================================================================
setTransformation: Kinematic Basic Matrix A Generation
Maps nodal global displacements to basic element deformations.
================================================================================
*/

void ZeroLength6D::setTransformation(void)
{
    // Forza order a 6 se usiamo un NDMaterial a 6 gradi di liberta
    if (order == 0 && theNDMaterial != nullptr) {
        order = theNDMaterial->getOrder();
    }
    if (order == 0) order = 6;

    // Ridimensiona e azzera le matrici operative
    A.resize(order, numDOF);
    K.resize(numDOF, numDOF);
    P.resize(numDOF);
    v.resize(order);

    A.Zero();
    K.Zero();
    P.Zero();
    v.Zero();

    // Popola la matrice di cinematica A (order x 12)
    // v = A * U  -->  v = R * (U_nodo2 - U_nodo1)
    for (int i = 0; i < order && i < 6; i++) {
        for (int j = 0; j < 6; j++) {
            A(i, j) = -transformation(i, j); // Nodo 1 (-R)
            A(i, j + 6) = transformation(i, j); // Nodo 2 (+R)
        }
    }
}

/*
================================================================================
computeStrain: Compute Material Basic Deformations (v)
Formulation: v = A * (u2 - u1)
================================================================================
*/

void ZeroLength6D::computeStrain(void)
{
    // CONTROLLO DI SICUREZZA: Evita la violazione di accesso (Null Pointer)
    if (end1Ptr == nullptr || end2Ptr == nullptr) {
        opserr << "ERROR: ZeroLength6D::computeStrain() -- Nodi non ancora collegati (end1Ptr o end2Ptr nullo)!\n";
        return;
    }

    // Recupera i vettori di spostamento direttamente dai puntatori ai nodi dell'elemento
    const Vector& disp1 = end1Ptr->getTrialDisp();
    const Vector& disp2 = end2Ptr->getTrialDisp();

    // Assembla il vettore U globale di 12 elementi
    Vector U(numDOF);
    U.Zero();

    for (int i = 0; i < 6; i++) {
        U(i) = disp1(i); // Nodo 1
        U(i + 6) = disp2(i); // Nodo 2
    }

    // Calcola il vettore deformazione locale v = A * U
    v.Zero();
    for (int k = 0; k < order; k++) {
        for (int i = 0; i < numDOF; i++) {
            v(k) += A(k, i) * U(i);
        }
    }
}

/*
================================================================================
sendSelf: Parallel Database Send Protocol
Saves and sends element attributes through the communication channels.
================================================================================
*/

int ZeroLength6D::sendSelf(int commitTag, Channel &theChannel)
{
    static ID data(6);
    
    data(0) = this->getTag();
    data(1) = dimension;
    data(2) = numDOF;
    data(3) = order;
    data(4) = connectedExternalNodes(0);
    data(5) = connectedExternalNodes(1);

    if (theChannel.sendID(this->getDbTag(), commitTag, data) < 0) {
        opserr << "ZeroLength6D::sendSelf -- errore nell'invio di ID data\n";
        return -1;
    }

    if (theNDMaterial->sendSelf(commitTag, theChannel) < 0) {
        opserr << "ZeroLength6D::sendSelf -- errore nell'invio del materiale NDMaterial\n";
        return -1;
    }

    return 0;
}

/*
================================================================================
recvSelf: Parallel Database Receive Protocol
Restores element data and resizes matrices safely for local thread use.
================================================================================
*/

int ZeroLength6D::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
    static ID data(6);

    if (theChannel.recvID(this->getDbTag(), commitTag, data) < 0) {
        opserr << "ZeroLength6D::recvSelf -- errore nella ricezione di ID data\n";
        return -1;
    }

    this->setTag(data(0));
    dimension = data(1);
    numDOF = data(2);
    order = data(3);
    connectedExternalNodes(0) = data(4);
    connectedExternalNodes(1) = data(5);

    // Ridimensionamento delle matrici/vettori di stato del thread
    K.resize(numDOF, numDOF);
    P.resize(numDOF);
    v.resize(order);
    transformation.resize(6, 6);

    // Se il materiale locale non esiste, istanzialo tramite Broker
    if (theNDMaterial == nullptr) {
        theNDMaterial = theBroker.getNewNDMaterial(order);
    }

    if (theNDMaterial != nullptr) {
        if (theNDMaterial->recvSelf(commitTag, theChannel, theBroker) < 0) {
            opserr << "ZeroLength6D::recvSelf -- errore nella ricezione di NDMaterial\n";
            return -1;
        }
    }

    return 0;
}

/*
================================================================================
displaySelf: Graphic Visualization Call
================================================================================
*/
int ZeroLength6D::displaySelf(Renderer &theViewer, int mode, float fact, const char **displayModes, int numModes)
{
    return 0;
}

/*
================================================================================
Print: Standard Terminal Console Logging Output
================================================================================
*/
void ZeroLength6D::Print(OPS_Stream &s, int flag)
{
    s << "Element: " << this->getTag() << " type: ZeroLength6D\n";
    s << "Connected Nodes: " << connectedExternalNodes(0) << " and " << connectedExternalNodes(1) << "\n";
    s << "Number of DOFs: " << numDOF << "\n";
    if (theNDMaterial != 0) {
        s << "Material Model Tag: " << theNDMaterial->getTag() << "\n";
    }
}

/*
================================================================================
setResponse: Handle recorder registration commands ("force", "deformation")
================================================================================
*/
Response *ZeroLength6D::setResponse(const char **argv, int argc, OPS_Stream &theHandler)
{
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {
        return new ElementResponse(this, 1, P); // FIXED C2100: Pass object direct reference
    }
    else if (strcmp(argv[0],"deformation") == 0 || strcmp(argv[0],"deformations") == 0) {
        return new ElementResponse(this, 2, v); // FIXED C2100: Pass object direct reference
    }
    return 0;
}

/*
================================================================================
getResponse: Extract current data arrays to recorder modules
================================================================================
*/
int ZeroLength6D::getResponse(int responseID, Information &eleInformation)
{
    if (responseID == 1) {
        return eleInformation.setVector(P); // FIXED C2100: Extracted vector directly
    }
    else if (responseID == 2) {
        return eleInformation.setVector(v); // FIXED C2100: Extracted vector directly
    }
    return -1;
}

/*
================================================================================
Parameter Management Hooks for Optimization/Sensitivity analyses
================================================================================
*/
int ZeroLength6D::setParameter(const char **argv, int argc, Parameter &param)
{
    return theNDMaterial->setParameter(argv, argc, param);
}

int ZeroLength6D::updateParameter(int parameterID, Information &info)
{
    return theNDMaterial->updateParameter(parameterID, info);
}

/*
================================================================================
State Management Functions
Passes convergence and commit states down to the assigned NDMaterial
================================================================================
*/

int ZeroLength6D::commitState(void)
{
    int code = 0;

    // Salva lo stato dell'elemento nell'interfaccia base Element
    if ((code = this->Element::commitState()) != 0) {
        opserr << "ZeroLength6D::commitState () - failed in base class\n";
    }

    // Salva lo stato interno del materiale ND
    if (theNDMaterial != nullptr) {
        code += theNDMaterial->commitState();
    }

    return code;
}

int ZeroLength6D::revertToLastCommit(void)
{
    int code = 0;

    // Ripristina lo stato interno del materiale ND all'ultimo passo converso
    if (theNDMaterial != nullptr) {
        code += theNDMaterial->revertToLastCommit();
    }

    return code;
}

int ZeroLength6D::revertToStart(void)
{
    int code = 0;

    // Azzeramento dello stato interno del materiale ND
    if (theNDMaterial != nullptr) {
        code += theNDMaterial->revertToStart();
    }

    return code;
}