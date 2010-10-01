#include <protomol/force/hessian/ReducedHessGB.h>
#include <protomol/topology/GenericTopology.h>

#include <protomol/base/Report.h>

#include <iomanip>
#include <iostream>

using namespace ProtoMol;
using namespace ProtoMol::Report;
using namespace std;


Matrix3By3 ReducedHessGB::operator()( Real a,
                                          const Vector3D &rij,
                                          const GenericTopology *topo,
                                          int atom1, int atom2,
                                          int numatoms,
                                          Real soluteDielec, Real solventDielec,
                                          ExclusionClass excl) const {




   Real dist = sqrt(a);

   //derivatives of burial term


   //Born radius
   Real bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
   Real bornRad_j = topo->atoms[atom2].myGBSA_T->bornRad;



   //cout << setprecision(10) <<"HessGBACE : Atom "<<atom1<<", Atom2 "<<atom2<<", bornRadiusDerivative_ij "<<bornRadiusDerivative_ij<<", bornRadiusDerivative_ji "<<bornRadiusDerivative_ji<<endl;

   int type1 = topo->atoms[atom1].type;
   int type2 = topo->atoms[atom2].type;

   Real radius_i = topo->atomTypes[type1].vdwR;
   Real radius_j = topo->atomTypes[type2].vdwR;

   //offset radii ({\tilde{\rho}_{j}})
   Real offsetRadius_i = topo->atomTypes[type1].vdwR - topo->atoms[atom1].myGBSA_T->offsetRadius;
   Real offsetRadius_j = topo->atomTypes[type2].vdwR - topo->atoms[atom2].myGBSA_T->offsetRadius;

   Real psi_i = topo->atoms[atom1].myGBSA_T->PsiValue;

   Real tanhparam_i = topo->alphaObc*psi_i - topo->betaObc*psi_i*psi_i + topo->gammaObc*psi_i*psi_i*psi_i;

   Real tanh_i = tanh(tanhparam_i);

   Real tanhparam_i_derv = topo->alphaObc - 2*topo->betaObc*psi_i + 3*topo->gammaObc*psi_i*psi_i;

   //data required for the second derivative of the burial term
   Real Lij = topo->atoms[atom1].myGBSA_T->Lvalues[atom2];
   Real Uij = topo->atoms[atom1].myGBSA_T->Uvalues[atom2];

   Real invLij = 1/Lij;
   Real invUij = 1/Uij;

   //Derivatives for calculation of the derivative of the born radii
   Real dLijdrij, dUijdrij, dCijdrij;

     //Scaling factors
     Real S_i = topo->atoms[atom1].myGBSA_T->scalingFactor;
     Real S_j = topo->atoms[atom2].myGBSA_T->scalingFactor;

  //Check Equations (11-13)
   if (offsetRadius_i <= (dist - S_j*offsetRadius_j)) dLijdrij = 1;
   else dLijdrij = 0;

   if (offsetRadius_i < (dist + S_j*offsetRadius_j)) dUijdrij = 1;
   else dUijdrij = 0;

   if (offsetRadius_i <= (S_j*offsetRadius_j - dist)) dCijdrij = 2*(1/Lij)*(1/Lij)*dLijdrij;
   else dCijdrij = 0;

    //tanhk derivative
     Real tanhparam_derv_i = (topo->alphaObc - 2*topo->betaObc*psi_i + 3*topo->gammaObc*psi_i*psi_i);

     Real S_i_term = (S_i*offsetRadius_i)/dist;
     Real S_j_term = (S_j*offsetRadius_j)/dist;


   //First Derivative of burial term
   Real btderv_ij = -0.5*dLijdrij*(1/(Lij*Lij)) + 0.5*dUijdrij*(1/(Uij*Uij)) + 0.125*((1/(Uij*Uij)) - (1/(Lij*Lij))) + 0.125*dist*((2/(Lij*Lij*Lij))*dLijdrij - (2/(Uij*Uij*Uij))*dUijdrij) - 0.25*(1/(dist*dist))*log(Lij/Uij) + (Uij/(4*dist*Lij))*((1/Uij)*dLijdrij - (Lij/(Uij*Uij))*dUijdrij) - 0.125*power(S_j_term,2)*((1/(Lij*Lij)) - (1/(Uij*Uij))) + 0.25*((S_j*S_j*offsetRadius_j*offsetRadius_j)/(dist*Uij*Uij*Uij))*dUijdrij - 0.25*((S_j*S_j*offsetRadius_j*offsetRadius_j)/(dist*Lij*Lij*Lij))*dLijdrij + dCijdrij;

   Real psiderv_i_ij = offsetRadius_i*btderv_ij;

   Real bornRadiusDerivative_ij = power(bornRad_i,2)*offsetRadius_i*(1-tanh_i*tanh_i)*tanhparam_derv_i*(1/radius_i)*btderv_ij;

   //second derivative of the burial term

   Real d2BTijdrij2 = power(invLij,3)*power(dLijdrij,2) - power(invUij,3)*power(dUijdrij,2) + 0.5*power(invLij,3)*dLijdrij - 0.5*power(invUij,3)*dUijdrij + (dist/4)*(3*power(invUij,4)*power(dUijdrij,2) - 3*power(invLij,4)*power(dLijdrij,2)) + 0.5*(1/power(dist,3))*log(Lij/Uij) + 0.5*invLij*(1/(dist*dist))*(Lij*invUij*dUijdrij - dLijdrij) + 0.5*(1/dist)*invLij*invUij*(Lij*invUij*dUijdrij*dUijdrij - dUijdrij*dLijdrij) + 0.25*(1/dist)*invLij*invLij*dLijdrij*(Lij*invUij*dUijdrij - dLijdrij) - 0.25*(1/dist)*invLij*invUij*dUijdrij*(Lij*invUij*dUijdrij - dLijdrij) + ((power(S_j,2)*power(offsetRadius_j,2))/(4*dist*dist*dist))*(invLij*invLij - invUij*invUij) - ((power(S_j,2)*power(offsetRadius_j,2))/(2*dist*dist))*(invUij*invUij*invUij*dUijdrij - invLij*invLij*invLij*dLijdrij) + ((power(S_j,2)*power(offsetRadius_j,2))/(4*dist))*(3*power(invLij,4)*power(dLijdrij,2) - 3*power(invUij,4)*power(dUijdrij,2));


  Real d2Psi_i_drij2 = d2BTijdrij2*offsetRadius_i;

  Real alpha = topo->alphaObc;
  Real beta = topo->betaObc;
  Real gamma = topo->gammaObc;

  //second derivative of born radius
  Real d2Ridrij2 =  2*(1 - tanh_i*tanh_i)*(1-tanh_i*tanh_i)*power(psiderv_i_ij,2)*power(tanhparam_i_derv,2)*(power(bornRad_i,3)/power(radius_i,2)) -2*power(bornRad_i,2)*tanh_i*(1 - tanh_i*tanh_i)*power(psiderv_i_ij,2)*power(tanhparam_i_derv,2)*(1/radius_i) + (power(bornRad_i,2)/radius_i)*(1 - power(tanh_i,2))*(alpha*d2Psi_i_drij2 - 2*beta*power(psiderv_i_ij,2)- 2*beta*psi_i*d2Psi_i_drij2) + (power(bornRad_i,2)/radius_i)*(1 - power(tanh_i,2))*(6*gamma*psi_i*power(psiderv_i_ij,2) + 3*gamma*power(psi_i,2)*d2Psi_i_drij2);

   //data required for the second derivative of the burial term
   Real Lji = topo->atoms[atom2].myGBSA_T->Lvalues[atom1];
   Real Uji = topo->atoms[atom2].myGBSA_T->Uvalues[atom1];

   Real invLji = 1/Lji;
   Real invUji = 1/Uji;


     //Derivatives for calculation of the derivative of the born radii
      Real dLjidrij, dUjidrij, dCjidrij;

      if (offsetRadius_j <= (dist - S_i*offsetRadius_i)) dLjidrij = 1;
      else dLjidrij = 0;

      if (offsetRadius_j < (dist + S_i*offsetRadius_i)) dUjidrij = 1;
      else dUjidrij = 0;

      if (offsetRadius_j <= (S_i*offsetRadius_i - dist)) dCjidrij = 2*(1/Lji)*(1/Lji)*dLjidrij;
      else dCjidrij = 0;

   Real psi_j = topo->atoms[atom2].myGBSA_T->PsiValue;
   Real tanhparam_j = topo->alphaObc*psi_j - topo->betaObc*psi_j*psi_j + topo->gammaObc*psi_j*psi_j*psi_j;
   Real tanh_j = tanh(tanhparam_j);
   Real tanhparam_j_derv = topo->alphaObc - 2*topo->betaObc*psi_j + 3*topo->gammaObc*psi_j*psi_j;

     //tanhi derivative
     Real tanhparam_derv_j = (topo->alphaObc - 2*topo->betaObc*psi_j + 3*topo->gammaObc*psi_j*psi_j);

   //first derivative of burial term
   Real btderv_ji = -0.5*dLjidrij*(1/(Lji*Lji)) + 0.5*dUjidrij*(1/(Uji*Uji)) + 0.125*((1/(Uji*Uji)) - (1/(Lji*Lji))) + 0.125*dist*((2/(Lji*Lji*Lji))*dLjidrij - (2/(Uji*Uji*Uji))*dUjidrij) - 0.25*(1/(dist*dist))*log(Lji/Uji) + (Uji/(4*dist*Lji))*((1/Uji)*dLjidrij - (Lji/(Uji*Uji))*dUjidrij) - 0.125*power(S_i_term,2)*((1/(Lji*Lji)) - (1/(Uji*Uji))) + 0.25*((S_i*S_i*offsetRadius_i*offsetRadius_i)/(dist*Uji*Uji*Uji))*dUjidrij - 0.25*((S_i*S_i*offsetRadius_i*offsetRadius_i)/(dist*Lji*Lji*Lji))*dLjidrij + dCjidrij;

   Real psiderv_j_ij = offsetRadius_j*btderv_ji;

   Real bornRadiusDerivative_ji = power(bornRad_j,2)*offsetRadius_j*(1-tanh_j*tanh_j)*tanhparam_derv_j*(1/radius_j)*btderv_ji;

   //second derivative of the burial term

   Real d2BTjidrji2 = power(invLji,3)*power(dLjidrij,2) - power(invUji,3)*power(dUjidrij,2) + 0.5*power(invLji,3)*dLjidrij - 0.5*power(invUji,3)*dUjidrij + (dist/4)*(3*power(invUji,4)*power(dUjidrij,2) - 3*power(invLji,4)*power(dLjidrij,2)) + 0.5*(1/power(dist,3))*log(Lji/Uji) + 0.5*invLji*(1/(dist*dist))*(Lji*invUji*dUjidrij - dLjidrij) + 0.5*(1/dist)*invLji*invUji*(Lji*invUji*dUjidrij*dUjidrij - dUjidrij*dLjidrij) + 0.25*(1/dist)*invLji*invLji*dLjidrij*(Lji*invUji*dUjidrij - dLjidrij) - 0.25*(1/dist)*invLji*invUji*dUjidrij*
(Lji*invUji*dUjidrij - dLjidrij) + ((power(S_i,2)*power(offsetRadius_i,2))/(4*dist*dist*dist))*(invLji*invLji - invUji*invUji) - ((power(S_i,2)*power(offsetRadius_i,2))/(2*dist*dist))*(invUji*invUji*invUji*dUjidrij - invLji*invLji*invLji*dLjidrij) + ((power(S_i,2)*power(offsetRadius_i,2))/(4*dist))*(3*power(invLji,4)*power(dLjidrij,2) - 3*power(invUji,4)*power(dUjidrij,2));

   Real d2Psi_j_drij2 = d2BTjidrji2*offsetRadius_j;

  //second derivative of born radius
  Real d2Rjdrij2 =  2*(1 - tanh_j*tanh_j)*(1-tanh_j*tanh_j)*power(psiderv_j_ij,2)*power(tanhparam_j_derv,2)*(power(bornRad_j,3)/power(radius_j,2)) -2*power(bornRad_j,2)*tanh_j*(1 - tanh_j*tanh_j)*power(psiderv_j_ij,2)*power(tanhparam_j_derv,2)*(1/radius_j) + (power(bornRad_j,2)/radius_j)*(1 - power(tanh_j,2))*(alpha*d2Psi_j_drij2 - 2*beta*power(psiderv_j_ij,2)- 2*beta*psi_j*d2Psi_j_drij2) + (power(bornRad_j,2)/radius_j)*(1 - power(tanh_j,2))*(6*gamma*psi_j*power(psiderv_j_ij,2) + 3*gamma*power(psi_j,2)*d2Psi_j_drij2);


   //grab f_ij 
   //Real fGB_ij = topo->atoms[atom1].myGBSA_T->fij[atom2];
     Real expterm_ij = (dist*dist)/(4*bornRad_i*bornRad_j);

     Real exp_w = exp(-expterm_ij);
   Real fGB_ij = sqrt(dist*dist + bornRad_i*bornRad_j*exp_w);


     //Real dfGBdrij = 0.5*(1/fGB_ij)*(exp(-expterm_ij)*(bornRadiusDerivative_ij*bornRad_j + bornRad_i*bornRadiusDerivative_ji-(dist/2) + ((dist*dist)/(4.0*bornRad_i)*bornRadiusDerivative_ij)+((dist*dist)/(4.0*bornRad_j)*bornRadiusDerivative_ji)));
   Real dfGBdrij = 0.5*(1/fGB_ij)*((2*dist - 0.5*exp(-expterm_ij)*dist) + exp(-expterm_ij)*bornRadiusDerivative_ij*(bornRad_j + (dist*dist)/(4*bornRad_i)) + exp(-expterm_ij)*bornRadiusDerivative_ji*(bornRad_i + (dist*dist)/(4*bornRad_j)));

/*
     //second derivative
     Real d2fGBdrij2 = (-1/fGB_ij)*dfGBdrij*dfGBdrij + (1/fGB_ij) + 0.5*(1/fGB_ij)*exp_w*(d2Ridrij2*bornRad_j + 2*bornRadiusDerivative_ij*bornRadiusDerivative_ji + 2*bornRadiusDerivative_ij*bornRad_j*(-(dist/(2*bornRad_i*bornRad_j)) + ((dist*dist)/(4*bornRad_i*bornRad_i*bornRad_j))*bornRadiusDerivative_ij + ((dist*dist)/(4*bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ji) + bornRad_i*d2Rjdrij2 + 2*bornRad_j*bornRadiusDerivative_ij*(-(dist/(2*bornRad_i*bornRad_j)) + ((dist*dist)/(4*bornRad_i*bornRad_i*bornRad_j))*bornRadiusDerivative_ij + ((dist*dist)/(4*bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ji) + bornRad_i*bornRad_j*(-(1/(2*bornRad_i*bornRad_j)) + ((dist)/(4*bornRad_i*bornRad_i*bornRad_j))*bornRadiusDerivative_ij + ((dist)/(4*bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ji - ((dist*dist)/(2*bornRad_i*bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ij*bornRadiusDerivative_ji ) + bornRad_i*bornRad_j*( ((dist*dist)/(4*bornRad_i*bornRad_i*bornRad_j))*d2Ridrij2 + ((dist*dist)/(4*bornRad_i*bornRad_j*bornRad_j))*d2Rjdrij2  - ((dist*dist)/(2*bornRad_i*power(bornRad_j,3)))*power(bornRadiusDerivative_ji,2) - ((dist*dist)/(2*bornRad_j*power(bornRad_i,3)))*power(bornRadiusDerivative_ij,2)) + bornRad_i*bornRad_j*(-(dist/(2*bornRad_i*bornRad_j)) + ((dist*dist)/(4*bornRad_i*bornRad_i*bornRad_j))*bornRadiusDerivative_ij + ((dist*dist)/(4*bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ji));
*/
   Real d2fGBdrij2 = (-1/fGB_ij)*dfGBdrij*dfGBdrij + (1/fGB_ij) + 0.5*(1/fGB_ij)*exp_w*(d2Ridrij2*bornRad_j + 2*bornRadiusDerivative_ij*bornRadiusDerivative_ji + 2*bornRadiusDerivative_ij*bornRad_j*(-(dist/(2*bornRad_i*bornRad_j)) + ((dist*dist)/(4*bornRad_i*bornRad_i*bornRad_j))*bornRadiusDerivative_ij + ((dist*dist)/(4*bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ji) + bornRad_i*d2Rjdrij2 + 2*bornRad_i*bornRadiusDerivative_ji*(-(dist/(2*bornRad_i*bornRad_j)) + ((dist*dist)/(4*bornRad_i*bornRad_i*bornRad_j))*bornRadiusDerivative_ij + ((dist*dist)/(4*bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ji)+ bornRad_i*bornRad_j*(-(1/(2*bornRad_i*bornRad_j)) + ((dist)/(bornRad_i*bornRad_i*bornRad_j))*bornRadiusDerivative_ij + ((dist)/(bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ji - ((dist*dist)/(2*bornRad_i*bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ij*bornRadiusDerivative_ji ) + bornRad_i*bornRad_j*( ((dist*dist)/(4*bornRad_i*bornRad_i*bornRad_j))*d2Ridrij2 + ((dist*dist)/(4*bornRad_i*bornRad_j*bornRad_j))*d2Rjdrij2  - ((dist*dist)/(2*bornRad_i*power(bornRad_j,3)))*power(bornRadiusDerivative_ji,2) - ((dist*dist)/(2*bornRad_j*power(bornRad_i,3)))*power(bornRadiusDerivative_ij,2)) + bornRad_i*bornRad_j*(-(dist/(2*bornRad_i*bornRad_j)) + ((dist*dist)/(4*bornRad_i*bornRad_i*bornRad_j))*bornRadiusDerivative_ij + ((dist*dist)/(4*bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ji )*(-(dist/(2*bornRad_i*bornRad_j)) + ((dist*dist)/(4*bornRad_i*bornRad_i*bornRad_j))*bornRadiusDerivative_ij + ((dist*dist)/(4*bornRad_i*bornRad_j*bornRad_j))*bornRadiusDerivative_ji) );

   //charges on atoms i and j
   Real charge_i = topo->atoms[atom1].scaledCharge;
   Real charge_j = topo->atoms[atom2].scaledCharge;

   Real second_derivative_born = charge_i*charge_j*(1/power(fGB_ij,2))*d2fGBdrij2 - (2/power(fGB_ij,3))*dfGBdrij*dfGBdrij;

   // self terms
   second_derivative_born += power(charge_i, 2) / power(bornRad_i, 2) * (d2Ridrij2 / 2 - power( bornRadiusDerivative_ij, 2) / bornRad_i);
   second_derivative_born += power(charge_j, 2) / power(bornRad_j, 2) * (d2Rjdrij2 / 2 - power( bornRadiusDerivative_ji, 2) / bornRad_j);

   Real first_derivative_born = (charge_i*charge_j)*(1/power(fGB_ij,2))*dfGBdrij;

   //self terms
   first_derivative_born += power(charge_i, 2) * bornRadiusDerivative_ij / (2 * power(bornRad_i, 2));
   first_derivative_born += power(charge_j, 2) * bornRadiusDerivative_ji / (2 * power(bornRad_j, 2));

   Real d2Gikterm = 0;
   Real dGikterm = 0;
   Real d2Gjkterm = 0;
   Real dGjkterm = 0;
   //Real dfGBik, d2fGBik, dfGBjk, d2fGBjk;

   for(int k = 0; k<numatoms; k++) {

     if ((k == atom1) || (k == atom2)) continue;

     d2Gikterm += SecondDerivativeFGB(topo, atom1, k, d2Ridrij2, bornRadiusDerivative_ij);
     d2Gjkterm += SecondDerivativeFGB(topo, atom2, k, d2Rjdrij2, bornRadiusDerivative_ji);

     dGikterm += FirstDerivativeFGB(topo, atom1, k, bornRadiusDerivative_ij);
     dGjkterm += FirstDerivativeFGB(topo, atom1, k, bornRadiusDerivative_ji); 
   }

   second_derivative_born += d2Gikterm + d2Gjkterm;
   second_derivative_born *= ((1/soluteDielec) - (1/solventDielec));

   first_derivative_born += dGikterm + dGjkterm;
   first_derivative_born *= ((1/soluteDielec) - (1/solventDielec))*(1/dist);


   Matrix3By3 I(1, 0 ,0 , 0 , 1, 0, 0 , 0, 1);
   Matrix3By3 vec_rij_ij(rij, rij);

   Matrix3By3 H((I - vec_rij_ij*(1/(dist*dist)))*first_derivative_born);

   H += vec_rij_ij*(1/(dist*dist))*second_derivative_born;

   return H;

}

Real ReducedHessGB::FirstDerivativeFGB(
                          const GenericTopology *topo,
                          const int atom1, const int atom2,
                          const Real bornRadiusDerivative_ij) const {

   //grab r_{ik} / r_{jk}
   Real dist = topo->atoms[atom1].myGBSA_T->distij[atom2];

   //born radius of i and k
   Real bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
   Real bornRad_k = topo->atoms[atom2].myGBSA_T->bornRad;

   Real expfactor = (dist*dist)/(4.0*bornRad_i*bornRad_k);
   Real exp_ik = exp(-expfactor);
   Real fGBik = sqrt(dist*dist + bornRad_i*bornRad_k*exp_ik);

   Real dfGBikdrij = 0.5*(1/fGBik)*exp_ik*bornRadiusDerivative_ij*(bornRad_k + ((dist*dist)/(4*bornRad_i)));

   Real charge_i = topo->atoms[atom1].scaledCharge;
   Real charge_k = topo->atoms[atom2].scaledCharge;

   Real dGBik = (charge_i*charge_k)*(1/power(fGBik,2))*dfGBikdrij;

   return dGBik;
}

Real ReducedHessGB::SecondDerivativeFGB(
                          const GenericTopology *topo,
                          int atom1, int atom2,
                          Real d2Ridrij2, Real bornRadiusDerivative_ij) const {

   //grab r_{ik} / r_{jk}
   Real dist = topo->atoms[atom1].myGBSA_T->distij[atom2];

   //report << plain <<"Seen ReducedHessGB : Atom1 "<<atom1<<", Atom2 "<<atom2<<" dist "<<topo->atoms[atom1].myGBSA_T->distij[atom2]<<" Lij "<<topo->atoms[atom1].myGBSA_T->Lvalues[atom2]<<endr;


   //born radius of i and k
   Real bornRad_i = topo->atoms[atom1].myGBSA_T->bornRad;
   Real bornRad_k = topo->atoms[atom2].myGBSA_T->bornRad;

   Real expfactor = (dist*dist)/(4.0*bornRad_i*bornRad_k);
   Real exp_ik = exp(-expfactor);
   Real fGBik = sqrt(dist*dist + bornRad_i*bornRad_k*exp_ik);

   Real dfGBikdrij = 0.5*(1/fGBik)*exp_ik*bornRadiusDerivative_ij*(bornRad_k + ((dist*dist)/(4*bornRad_i)));

   //Real d2fGBikdrij2 = -(1/fGBik)*dfGBikdrij*dfGBikdrij + (1/(2.0*fGBik))*exp_ik*(d2Ridrij2*(bornRad_k + ((dist*dist)/(4.0*bornRad_i))) + (power(dist,4)/(16*power(bornRad_i,3)*bornRad_k))*dfGBikdrij*dfGBikdrij);
   Real d2fGBikdrij2 = -(1/fGBik)*dfGBikdrij*dfGBikdrij + (1/(2.0*fGBik))*exp_ik*(d2Ridrij2*(bornRad_k + ((dist*dist)/(4.0*bornRad_i))) + (power(dist,4)/(16*power(bornRad_i,3)*bornRad_k))*bornRadiusDerivative_ij*bornRadiusDerivative_ij);

   //charge of i/j
   Real charge_i = topo->atoms[atom1].scaledCharge;
   Real charge_k = topo->atoms[atom2].scaledCharge;

   Real d2GBik = (charge_i*charge_k)*d2fGBikdrij2 - (2/power(fGBik,3))*dfGBikdrij*dfGBikdrij;

   return d2GBik;

}

