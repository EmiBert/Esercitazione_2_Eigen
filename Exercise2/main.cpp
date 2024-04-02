#include <iostream>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


void SolveSystemPALU(const Matrix2d& A, const Vector2d& b, const Vector2d& exactSolution){

    Vector2d x = A.fullPivLu().solve(b); //risoluzione del sistema tramite la decomposizione PALU
    double errRel = (exactSolution - x).norm() / exactSolution.norm();
    cout<<"Risoluzione col metodo PALU: errore relativo= "<<errRel<<endl;

}

void SolveSystemQR(const Matrix2d& A, const Vector2d& b, const Vector2d& exactSolution){

    Vector2d x = A.fullPivHouseholderQr().solve(b); //risoluzione del sistema tramite la decomposizione QR
    double errRel = (exactSolution - x).norm() / exactSolution.norm();
    cout<<"Risoluzione col metodo QR:   errore relativo= "<<errRel<<"\n"<<endl;

}


int main()

{   Matrix2d M[3] = {      //vettore ontenete le 3 matrici "A" assegnate

        Matrix2d {{5.547001962252291e-01,-3.770900990025203e-02},
                 {8.320502943378437e-01,-9.992887623566787e-01}},

        Matrix2d {{5.547001962252291e-01,-5.540607316466765e-01},
                 {8.320502943378437e-01,-8.324762492991313e-01}},

        Matrix2d {{5.547001962252291e-01,-5.547001955851905e-01},
                 {8.320502943378437e-01,-8.320502947645361e-01}}
    };


    Vector2d TermineNoto[3] {       //vettore contente i 3 termini noti "b" assegnati

        Vector2d {-5.169911863249772e-01, 1.672384680188350e-01},

        Vector2d {-6.394645785530173e-04, 4.259549612877223e-04},

        Vector2d {-6.400391328043042e-10, 4.266924591433963e-10}
    };

    const Vector2d sol {-1.0e+0,-1.0e+00}; //soluzione esatta

    int i;

    for (i=0; i<3; i++){
        cout<<"Risoluzione della matrice numero "<<i<<endl;
        double detA = M[i].determinant();
        if (abs(detA) < 1e-16){     //controllo che il sistema ammetta un'unica soluzione
            cout<<"errore: la matrice "<<i<<" e' singolare \n"<<endl;
        }
        else{
            SolveSystemPALU(M[i],TermineNoto[i],sol);
            SolveSystemQR(M[i],TermineNoto[i],sol);
        }

    }

    return 0;
}



