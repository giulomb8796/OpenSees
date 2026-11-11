/*
Developed by: Davide Noè Gorini (davidenoe.gorini@unitn.it)
Date: 2023
Updated for OpenSees 3.8.0 Compliance (Thread-Safe & Multi-Dimensional Fixes)

Description:
This file contains the class declaration for ZeroLength6D,
which extends the standard ZeroLengthND element by providing:
1. Full 6D response - 3 translations and 3 rotations.
2. A dedicated method to compute relative velocity between end nodes,
   essential when pairing the element with rate-dependent (viscoplastic/damping) materials.

Reference:
Gorini, D.N. and Callisto, L. (2023): "A multiaxial inertial macroelement for deep foundations",
Computers and Geotechnics, vol. 155, doi: 10.1016/j.compgeo.2022.105222
*/

#ifndef ZeroLength6D_h
#define ZeroLength6D_h

#include <Element.h>
#include <Matrix.h>
#include <Vector.h>
#include <ID.h>

// Tolerance for zero length of element
#define LENTOL 1.0e-6

// Forward declarations to ensure compilation independence
class Node;
class Channel;
class NDMaterial;
class Response;
class Parameter;
class FEM_ObjectBroker;
class OPS_Stream;
class Domain;
class ElementalLoad;
class Renderer;
class Information;

class ZeroLength6D : public Element
{
public:
    // Costruttore standard (Coerente con il .cpp: accetta NDMaterial &)
    ZeroLength6D(int tag,
                 int dimension,
                 int Nd1, int Nd2,
                 const Vector& x,
                 const Vector& yprime,
                 NDMaterial& theNDMaterial);

    // Costruttore vuoto per FEM_ObjectBroker / Parallelismo
    ZeroLength6D();
    
    // Distruttore
    ~ZeroLength6D() override;

    // Metodi informativi su DOF e connettività    
    int getNumExternalNodes(void) const override;
    const ID& getExternalNodes(void) override;
    Node** getNodePtrs(void) override;

    int getNumDOF(void) override;
    void setDomain(Domain* theDomain) override;

    // Metodi di gestione dello stato (presenti nel .cpp) 
    int commitState(void) override;
    int revertToLastCommit(void) override;
    int revertToStart(void) override;

    // Matrici di rigidezza, smorzamento e massa    
    const Matrix& getTangentStiff(void) override;
    const Matrix& getInitialStiff(void) override;
    const Matrix& getDamp(void) override;
    const Matrix& getMass(void) override;

    // Gestione carichi elementari
    void zeroLoad(void) override;
    int addLoad(ElementalLoad* theLoad, double loadFactor) override;
    int addInertiaLoadToUnbalance(const Vector& accel) override;

    // Vettori di forza resistente
    const Vector& getResistingForce(void) override;
    const Vector& getResistingForceIncInertia(void) override;

    // Calcolo della velocità relativa locale (Metodo custom definito nel .cpp)
    const Vector& setTrialRelVel(void);

    // Serializzazione e Database (Parallel/Database OpenSees)
    int sendSelf(int commitTag, Channel& theChannel) override;
    int recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker) override;
    int displaySelf(Renderer& theViewer, int mode, float fact, const char** displayModes = 0, int numModes = 0) override;
    void Print(OPS_Stream& s, int flag = 0) override;

    // Recorders, Sensibilità e Parametri
    Response* setResponse(const char** argv, int argc, OPS_Stream& theHandler) override;
    int getResponse(int responseID, Information& eleInformation) override;
    int setParameter(const char** argv, int argc, Parameter& param) override;
    int updateParameter(int parameterID, Information& info) override;

protected:

private:
    // Metodi privati cinematici (chiamati nel .cpp)
    void setUp(int Nd1, int Nd2, const Vector& x, const Vector& yprime);
    void setTransformation(void);
    void computeStrain(void);

    // Attributi privati (Tutti utilizzati nel .cpp)
    ID connectedExternalNodes;         // Contiene i tag dei nodi (dimensione 2)
    int dimension;                     // Dimensione del problema (2D o 3D)
    int numDOF;                        // Numero totale di DOFs (6 o 12)
    
    Matrix transformation;             // Matrice di rotazione locale (6x6)
    Matrix A;                          // Matrice cinematica locale (order x numDOF)

    Matrix K;                          // Matrice di rigidezza globale
    Vector P;                          // Vettore delle forze resistenti globali
    Vector v;                          // Vettore deformazione per NDMaterial
    Vector RelVelNodes;                // Vettore velocità relativa tra nodi

    Node* end1Ptr;                     // Puntatore al Nodo I (Nodo 1)
    Node* end2Ptr;                     // Puntatore al Nodo J (Nodo 2)

    NDMaterial* theNDMaterial;         // Copia locale dell'NDMaterial
    int order;                         // Numero di componenti di sforzo/deformazione
};

#endif