/*
 * quad2d.h
 *   Implementation for four-node quad element shapes.
 *
 * Copyright (c) 2007  David Bindel
 * See the file COPYING for copying permissions
 */
#include "quad2d.h"

#define ME Quad2d


ME::~ME()
{
}


void ME::set_node(int nodenum, const double* x)
{
    nodex[2*nodenum+0] = x[0];
    nodex[2*nodenum+1] = x[1];
}


void ME::set_nodes(const double* x)
{
    for (int i = 0; i < 8; ++i)
        nodex[i] = x[i];
}


void ME::eval(const double *XX, double* xx,
              double* N, double* dN,
              double& J) const
{
    double X = XX[0];
    double Y = XX[1];
    double FF[4];

    /* <generator matexpr>
    // Evaluate element functions

    input X, Y;
    input nodex(2,4);

    output N(4);
    output dN(4,2);
    output xx(2);
    output FF(2,2);

    N1x = (1-X)/2;  N1y = (1-Y)/2;
    N2x = (1+X)/2;  N2y = (1+Y)/2;    

    N  = [ N1x*N1y; N2x*N1y; N2x*N2y; N1x*N2y];
    dN = [    -N1y,     N1y,     N2y,    -N2y;
              -N1x,    -N2x,     N2x,     N1x ]'/2;

    xx = nodex*N;
    FF = nodex*dN;
    */
    /* <generated matexpr> */ {
    double tmp1_ = X;
    double tmp2_ = Y;
    double tmp3_ = nodex[0*2+0];
    double tmp4_ = nodex[0*2+1];
    double tmp5_ = nodex[1*2+0];
    double tmp6_ = nodex[1*2+1];
    double tmp7_ = nodex[2*2+0];
    double tmp8_ = nodex[2*2+1];
    double tmp9_ = nodex[3*2+0];
    double tmp10_ = nodex[3*2+1];
    double tmp12_ = 1.0 - tmp1_;
    double tmp14_ = tmp12_ / 2.0;
    double tmp15_ = 1.0 - tmp2_;
    double tmp16_ = tmp15_ / 2.0;
    double tmp17_ = 1.0 + tmp1_;
    double tmp18_ = tmp17_ / 2.0;
    double tmp19_ = 1.0 + tmp2_;
    double tmp20_ = tmp19_ / 2.0;
    double tmp21_ = tmp14_ * tmp16_;
    double tmp22_ = tmp18_ * tmp16_;
    double tmp23_ = tmp18_ * tmp20_;
    double tmp24_ = tmp14_ * tmp20_;
    double tmp25_ = -tmp16_;
    double tmp26_ = -tmp20_;
    double tmp27_ = -tmp14_;
    double tmp28_ = -tmp18_;
    double tmp29_ = tmp25_ / 2.0;
    double tmp30_ = tmp16_ / 2.0;
    double tmp31_ = tmp20_ / 2.0;
    double tmp32_ = tmp26_ / 2.0;
    double tmp33_ = tmp27_ / 2.0;
    double tmp34_ = tmp28_ / 2.0;
    double tmp35_ = tmp18_ / 2.0;
    double tmp36_ = tmp14_ / 2.0;
    double tmp37_ = tmp3_ * tmp21_;
    double tmp38_ = tmp5_ * tmp22_;
    double tmp39_ = tmp37_ + tmp38_;
    double tmp40_ = tmp7_ * tmp23_;
    double tmp41_ = tmp39_ + tmp40_;
    double tmp42_ = tmp9_ * tmp24_;
    double tmp43_ = tmp41_ + tmp42_;
    double tmp44_ = tmp4_ * tmp21_;
    double tmp45_ = tmp6_ * tmp22_;
    double tmp46_ = tmp44_ + tmp45_;
    double tmp47_ = tmp8_ * tmp23_;
    double tmp48_ = tmp46_ + tmp47_;
    double tmp49_ = tmp10_ * tmp24_;
    double tmp50_ = tmp48_ + tmp49_;
    double tmp51_ = tmp3_ * tmp29_;
    double tmp52_ = tmp5_ * tmp30_;
    double tmp53_ = tmp51_ + tmp52_;
    double tmp54_ = tmp7_ * tmp31_;
    double tmp55_ = tmp53_ + tmp54_;
    double tmp56_ = tmp9_ * tmp32_;
    double tmp57_ = tmp55_ + tmp56_;
    double tmp58_ = tmp4_ * tmp29_;
    double tmp59_ = tmp6_ * tmp30_;
    double tmp60_ = tmp58_ + tmp59_;
    double tmp61_ = tmp8_ * tmp31_;
    double tmp62_ = tmp60_ + tmp61_;
    double tmp63_ = tmp10_ * tmp32_;
    double tmp64_ = tmp62_ + tmp63_;
    double tmp65_ = tmp3_ * tmp33_;
    double tmp66_ = tmp5_ * tmp34_;
    double tmp67_ = tmp65_ + tmp66_;
    double tmp68_ = tmp7_ * tmp35_;
    double tmp69_ = tmp67_ + tmp68_;
    double tmp70_ = tmp9_ * tmp36_;
    double tmp71_ = tmp69_ + tmp70_;
    double tmp72_ = tmp4_ * tmp33_;
    double tmp73_ = tmp6_ * tmp34_;
    double tmp74_ = tmp72_ + tmp73_;
    double tmp75_ = tmp8_ * tmp35_;
    double tmp76_ = tmp74_ + tmp75_;
    double tmp77_ = tmp10_ * tmp36_;
    double tmp78_ = tmp76_ + tmp77_;
    N[0*4+0] = tmp21_;
    N[0*4+1] = tmp22_;
    N[0*4+2] = tmp23_;
    N[0*4+3] = tmp24_;
    dN[0*4+0] = tmp29_;
    dN[0*4+1] = tmp30_;
    dN[0*4+2] = tmp31_;
    dN[0*4+3] = tmp32_;
    dN[1*4+0] = tmp33_;
    dN[1*4+1] = tmp34_;
    dN[1*4+2] = tmp35_;
    dN[1*4+3] = tmp36_;
    xx[0*2+0] = tmp43_;
    xx[0*2+1] = tmp50_;
    FF[0*2+0] = tmp57_;
    FF[0*2+1] = tmp64_;
    FF[1*2+0] = tmp71_;
    FF[1*2+1] = tmp78_;
    } /* </generated matexpr> */
    remap_gradients(FF, dN, J);
}


void ME::remap_gradients(double* FF, double* dN, double& J) const
{
    J = (FF[0]*FF[3]-FF[1]*FF[2]);
    double invF[4] = {  FF[3]/J, -FF[1]/J, 
		       -FF[2]/J,  FF[0]/J };

    /* <generator matexpr>
    // Remap gradients

    inout dN(4,2);
    input invF(2,2);
    dN = dN*invF;
    */
    /* <generated matexpr> */ {
    double tmp1_ = dN[0*4+0];
    double tmp2_ = dN[0*4+1];
    double tmp3_ = dN[0*4+2];
    double tmp4_ = dN[0*4+3];
    double tmp5_ = dN[1*4+0];
    double tmp6_ = dN[1*4+1];
    double tmp7_ = dN[1*4+2];
    double tmp8_ = dN[1*4+3];
    double tmp9_ = invF[0*2+0];
    double tmp10_ = invF[0*2+1];
    double tmp11_ = invF[1*2+0];
    double tmp12_ = invF[1*2+1];
    double tmp13_ = tmp1_ * tmp9_;
    double tmp14_ = tmp5_ * tmp10_;
    double tmp15_ = tmp13_ + tmp14_;
    double tmp16_ = tmp2_ * tmp9_;
    double tmp17_ = tmp6_ * tmp10_;
    double tmp18_ = tmp16_ + tmp17_;
    double tmp19_ = tmp3_ * tmp9_;
    double tmp20_ = tmp7_ * tmp10_;
    double tmp21_ = tmp19_ + tmp20_;
    double tmp22_ = tmp4_ * tmp9_;
    double tmp23_ = tmp8_ * tmp10_;
    double tmp24_ = tmp22_ + tmp23_;
    double tmp25_ = tmp1_ * tmp11_;
    double tmp26_ = tmp5_ * tmp12_;
    double tmp27_ = tmp25_ + tmp26_;
    double tmp28_ = tmp2_ * tmp11_;
    double tmp29_ = tmp6_ * tmp12_;
    double tmp30_ = tmp28_ + tmp29_;
    double tmp31_ = tmp3_ * tmp11_;
    double tmp32_ = tmp7_ * tmp12_;
    double tmp33_ = tmp31_ + tmp32_;
    double tmp34_ = tmp4_ * tmp11_;
    double tmp35_ = tmp8_ * tmp12_;
    double tmp36_ = tmp34_ + tmp35_;
    dN[0*4+0] = tmp15_;
    dN[0*4+1] = tmp18_;
    dN[0*4+2] = tmp21_;
    dN[0*4+3] = tmp24_;
    dN[1*4+0] = tmp27_;
    dN[1*4+1] = tmp30_;
    dN[1*4+2] = tmp33_;
    dN[1*4+3] = tmp36_;
    } /* </generated matexpr> */
}
