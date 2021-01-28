/*--
    Fresnel.cpp  

    This file is part of the Cornucopia curve sketching library.
    Copyright (C) 2010 Ilya Baran (baran37@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
The code in this file is adapted from:
Cephes Math Library Release 2.8:  June, 2000
Copyright 1984, 1987, 1989, 2000 by Stephen L. Moshier
*/

#include "Fresnel.h"
#include <vector>
#include <iostream>

using namespace std;
using namespace Eigen;
using namespace internal; //use Eigen's internal namespace for packet math
NAMESPACE_Cornu

//Polynomial evaluation routines
template<typename Scalar>
static Scalar polevl( const Scalar &x, const Matrix<Scalar, Dynamic, 1> &coefVec ) //regular
{
    int i = (int)coefVec.size() - 1;
    const Scalar *coef = &(coefVec[0]);
    Scalar ans = *coef++;

    do
    {
        ans = ans * x  +  *coef++;
    } while( --i );

    return ans;
}

template<typename Scalar>
static Scalar p1evl( const Scalar &x, const Matrix<Scalar, Dynamic, 1> &coefVec ) //leading coef is 1
{
    int i = (int)coefVec.size() - 1;
    const Scalar *coef = &(coefVec[0]);

    Scalar ans = x + *coef++;

    do
    {
        ans = ans * x  + *coef++;
    } while( --i );

    return ans;
}

//==================Coefficients===========================

//double precision rational coefficients for s, c, f, and g
static VectorXd dsn(6), dsd(6);
static VectorXd dcn(6), dcd(7);
static VectorXd dfn(10), dfd(10);
static VectorXd dgn(11), dgd(11);
//single precision polynomial coefficients
static VectorXf ssn(7);
static VectorXf scn(7);
static VectorXf sfn(8);
static VectorXf sgn(8);
//double precision polynomial coefficients
static VectorXd dssn(7);
static VectorXd dscn(7);
static VectorXd dsfn(8);
static VectorXd dsgn(8);

struct InitCoefs
{
    InitCoefs()
    {
        initPolynomial();
        initRational();
        ssn = dssn.cast<float>();
        scn = dscn.cast<float>();
        sfn = dsfn.cast<float>();
        sgn = dsgn.cast<float>();
    }

    void initPolynomial()
    {
        dssn <<
            1.647629463788700E-009,
            -1.522754752581096E-007,
            8.424748808502400E-006,
            -3.120693124703272E-004,
            7.244727626597022E-003,
            -9.228055941124598E-002,
            5.235987735681432E-001;
        dscn <<
            1.416802502367354E-008,
            -1.157231412229871E-006,
            5.387223446683264E-005,
            -1.604381798862293E-003,
            2.818489036795073E-002,
            -2.467398198317899E-001,
            9.999999760004487E-001;
        dsfn <<
            -1.903009855649792E+012,
            1.355942388050252E+011,
            -4.158143148511033E+009,
            7.343848463587323E+007,
            -8.732356681548485E+005,
            8.560515466275470E+003,
            -1.032877601091159E+002,
            2.999401847870011E+000;
        dsgn <<
            -1.860843997624650E+011,
            1.278350673393208E+010,
            -3.779387713202229E+008,
            6.492611570598858E+006,
            -7.787789623358162E+004,
            8.602931494734327E+002,
            -1.493439396592284E+001,
            9.999841934744914E-001;
    }

    void initRational()
    {
        dsn <<
            -2.99181919401019853726E3,
            7.08840045257738576863E5,
            -6.29741486205862506537E7,
            2.54890880573376359104E9,
            -4.42979518059697779103E10,
            3.18016297876567817986E11;
        dsd <<
            /* 1.00000000000000000000E0,*/
            2.81376268889994315696E2,
            4.55847810806532581675E4,
            5.17343888770096400730E6,
            4.19320245898111231129E8,
            2.24411795645340920940E10,
            6.07366389490084639049E11;
        dcn <<
            -4.98843114573573548651E-8,
            9.50428062829859605134E-6,
            -6.45191435683965050962E-4,
            1.88843319396703850064E-2,
            -2.05525900955013891793E-1,
            9.99999999999999998822E-1;
        dcd <<
            3.99982968972495980367E-12,
            9.15439215774657478799E-10,
            1.25001862479598821474E-7,
            1.22262789024179030997E-5,
            8.68029542941784300606E-4,
            4.12142090722199792936E-2,
            1.00000000000000000118E0;
        dfn <<
            4.21543555043677546506E-1,
            1.43407919780758885261E-1,
            1.15220955073585758835E-2,
            3.45017939782574027900E-4,
            4.63613749287867322088E-6,
            3.05568983790257605827E-8,
            1.02304514164907233465E-10,
            1.72010743268161828879E-13,
            1.34283276233062758925E-16,
            3.76329711269987889006E-20;
        dfd <<
            /*  1.00000000000000000000E0,*/
            7.51586398353378947175E-1,
            1.16888925859191382142E-1,
            6.44051526508858611005E-3,
            1.55934409164153020873E-4,
            1.84627567348930545870E-6,
            1.12699224763999035261E-8,
            3.60140029589371370404E-11,
            5.88754533621578410010E-14,
            4.52001434074129701496E-17,
            1.25443237090011264384E-20;
        dgn <<
            5.04442073643383265887E-1,
            1.97102833525523411709E-1,
            1.87648584092575249293E-2,
            6.84079380915393090172E-4,
            1.15138826111884280931E-5,
            9.82852443688422223854E-8,
            4.45344415861750144738E-10,
            1.08268041139020870318E-12,
            1.37555460633261799868E-15,
            8.36354435630677421531E-19,
            1.86958710162783235106E-22;
        dgd <<
            /*  1.00000000000000000000E0,*/
            1.47495759925128324529E0,
            3.37748989120019970451E-1,
            2.53603741420338795122E-2,
            8.14679107184306179049E-4,
            1.27545075667729118702E-5,
            1.04314589657571990585E-7,
            4.60680728146520428211E-10,
            1.10273215066240270757E-12,
            1.38796531259578871258E-15,
            8.39158816283118707363E-19,
            1.86958710162783236342E-22;
    }
} init;

//full double precision accuracy using rational functions
void fresnel( double xxa, double *ssa, double *cca )
{
    double f, g, cc, ss, c, s, t, u;
    double x, x2;

    x = fabs(xxa);
    x2 = x * x;
    if( x2 < 2.5625 )
    {
        t = x2 * x2;
        ss = x * x2 * polevl( t, dsn)/p1evl( t, dsd);
        cc = x * polevl( t, dcn)/polevl(t, dcd);
        goto done;
    }

    if( x > 36974.0 )
    {
        cc = 0.5;
        ss = 0.5;
        goto done;
    }

    /*        Asymptotic power series auxiliary functions
    *        for large argument
    */
    x2 = x * x;
    t = PI * x2;
    u = 1.0/(t * t);
    t = 1.0/t;
    f = 1.0 - u * polevl( u, dfn)/p1evl(u, dfd);
    g = t * polevl( u, dgn)/p1evl(u, dgd);

    t = HALFPI * x2;
    c = cos(t);
    s = sin(t);
    t = PI * x;
    cc = 0.5  +  (f * s  -  g * c)/t;
    ss = 0.5  -  (f * c  +  g * s)/t;

done:
    if( xxa < 0.0 )
    {
        cc = -cc;
        ss = -ss;
    }

    *cca = cc;
    *ssa = ss;
}

//roughly single-precision accuracy, using polynomial approximations
void fresnelApprox( double xxa, double *ssa, double *cca )
{
    double f, g, cc, ss, c, s, t, u;
    double x, x2;

    x = fabs(xxa);
    x2 = x * x;
    if( x2 < 2.5625 )
    {
        t = x2 * x2;
        ss = x * x2 * polevl( t, dssn);
        cc = x * polevl( t, dscn);
        goto done;
    }

    if( x > 36974.0 )
    {
        cc = 0.5;
        ss = 0.5;
        goto done;
    }

    /*        Asymptotic power series auxiliary functions
    *        for large argument
    */
    x2 = x * x;
    t = PI * x2;
    u = 1.0/(t * t);
    t = 1.0/t;
    f = 1.0 - u * polevl( u, dsfn);
    g = t * polevl( u, dsgn);

    t = HALFPI * x2;
    c = cos(t);
    s = sin(t);
    t = PI * x;
    cc = 0.5  +  (f * s  -  g * c)/t;
    ss = 0.5  -  (f * c  +  g * s)/t;

done:
    if( xxa < 0.0 )
    {
        cc = -cc;
        ss = -ss;
    }

    *cca = cc;
    *ssa = ss;
}

//The double precision version is not vectorized because the scalar version is actually faster
void fresnel(const VectorXd &t, VectorXd *s, VectorXd *c)
{
    s->resize(t.size());
    c->resize(t.size());
    for(int i = 0; i < t.size(); ++i)
        fresnel(t[i], &((*s)[i]), &((*c)[i]));
}

//Vectorization stuff
#ifdef EIGEN_VECTORIZE_SSE
typedef Packet4f PSETParam;

//vectorized polynomial evaluation
template<typename Packet>
static Packet vecpolevl( const Packet &x, const Matrix<typename unpacket_traits<Packet>::type, Dynamic, 1> &coefVec )
{
    int i = (int)coefVec.size() - 1;
    const typename unpacket_traits<Packet>::type *coef = &(coefVec[0]);
    Packet ans = pset1<PSETParam>(*coef++);

    do
    {
        ans = padd(pmul(ans, x), pset1<PSETParam>(*coef++));
    } while(--i);

    return ans;
}

template<typename Packet>
static Packet vecp1evl( const Packet &x, const Matrix<typename unpacket_traits<Packet>::type, Dynamic, 1> &coefVec )
{
    int i = (int)coefVec.size() - 1;
    const typename unpacket_traits<Packet>::type *coef = &(coefVec[0]);
    Packet ans = padd(x, pset1<PSETParam>(*coef++));

    do
    {
        ans = padd(pmul(ans, x), pset1<PSETParam>(*coef++));
    } while(--i);

    return ans;
}

//other vectorized utilities
EIGEN_STRONG_INLINE Packet4f packetTransferSign(const Packet4f& to, const Packet4f& from)
{
  const Packet4f mask = _mm_and_ps(from, _mm_castsi128_ps(_mm_setr_epi32(0x80000000,0x80000000,0x80000000,0x80000000)));
  return _mm_or_ps(to, mask);
}
EIGEN_STRONG_INLINE Packet2d packetTransferSign(const Packet2d& to, const Packet2d from)
{
  const Packet2d mask = _mm_and_pd(from, _mm_castsi128_pd(_mm_setr_epi32(0x00000000,0x80000000,0x00000000,0x80000000)));
  return _mm_or_pd(to, mask);
}

EIGEN_STRONG_INLINE Packet4f packetFmod(const Packet4f &a, const Packet4f &b)
{
    Packet4f div = pdiv(a, b);
    div = _mm_cvtepi32_ps(_mm_cvtps_epi32(div));
    return psub(a, pmul(div, b));
}

//vectorized for the low branch of the Fresnel approximation
void fresnelLow( const Packet4f &xxa, Packet4f *ssa, Packet4f *cca )
{
    Packet4f cc, ss, t;
    Packet4f x, x2;

    x = pabs(xxa);
    x2 = pmul(x, x);

    t = pmul(x2, x2);
    ss = pmul(x, pmul(x2, vecpolevl( t, ssn)));
    cc = pmul(x, vecpolevl( t, scn));

    *ssa = packetTransferSign(ss, xxa);
    *cca = packetTransferSign(cc, xxa);
}

//vectorized for the high branch
void fresnelMed( const Packet4f &xxa, Packet4f *ssa, Packet4f *cca )
{
    Packet4f cc, ss, t, u, f, g, c, s;
    Packet4f x, x2;

    x = pabs(xxa);
    x2 = pmul(x, x);

    t = pmul(pset1<PSETParam>(float(PI)), x2);
    t = pdiv(pset1<PSETParam>(float(1.0)), t);
    u = pmul(t, t);
    f = psub(pset1<PSETParam>(float(1.0)), pmul(u, vecpolevl( u, sfn)));
    g = pmul(t, vecpolevl( u, sgn));

    t = pmul(pset1<PSETParam>(float(HALFPI)), x2);

    //The following line is necessary because Eigen's psin and pcos don't handle large
    //inputs well.
    t = packetFmod(t, pset1<PSETParam>(float(TWOPI))); 

    c = pcos(t);
    s = psin(t);

    t = pdiv(pset1<PSETParam>(float(1. / PI)), x);
    cc = padd(pset1<PSETParam>(float(0.5)), pmul(t, psub(pmul(f, s), pmul(g, c))));
    ss = psub(pset1<PSETParam>(float(0.5)), pmul(t, padd(pmul(f, c), pmul(g, s))));

    *ssa = packetTransferSign(ss, xxa);
    *cca = packetTransferSign(cc, xxa);
}

//This version is vectorized
void fresnelApprox(const VectorXd &t, VectorXd *s, VectorXd *c)
{
    const int packetSize = 4;
    typedef Packet4f Packet;
    typedef float Scalar;

    typedef Matrix<float, packetSize, 1> ValVec;
    typedef Matrix<int, packetSize, 1> IdxVec;

    s->resize(t.size());
    c->resize(t.size());

    ValVec lowVal, medVal;
    IdxVec lowIdx, medIdx;
    Packet pval, ps, pc;
    ValVec vs, vc;
    int lowNum = 0, medNum = 0;

    for(int i = 0; i < t.size(); ++i)
    {
        double vsq = t[i] * t[i];
        if(vsq > 1367076676)
        {
            if(t[i] < 0)
                (*s)[i] = (*c)[i] = -0.5;
            else
                (*s)[i] = (*c)[i] = 0.5;

            continue;
        }
        if(vsq < 2.5625) //low
        {
            lowVal[lowNum] = Scalar(t[i]);
            lowIdx[lowNum++] = i;

            if(lowNum == packetSize)
            {
                pval = lowVal.packet<Aligned>(0);
                fresnelLow(pval, &ps, &pc);
                vs.writePacket<Aligned>(0, ps);
                vc.writePacket<Aligned>(0, pc);
                for(int j = 0; j < packetSize; ++j)
                {
                    (*s)[lowIdx[j]] = vs[j];
                    (*c)[lowIdx[j]] = vc[j];
                }
                lowNum = 0;
            }
        }
        else //med
        {
            medVal[medNum] = Scalar(t[i]);
            medIdx[medNum++] = i;

            if(medNum == packetSize)
            {
                pval = medVal.packet<Aligned>(0);
                fresnelMed(pval, &ps, &pc);
                vs.writePacket<Aligned>(0, ps);
                vc.writePacket<Aligned>(0, pc);
                for(int j = 0; j < packetSize; ++j)
                {
                    (*s)[medIdx[j]] = vs[j];
                    (*c)[medIdx[j]] = vc[j];
                }
                medNum = 0;
            }
        }
    }

    //finish up
    for(int i = 0; i < lowNum; ++i)
        fresnelApprox(lowVal[i], &((*s)[lowIdx[i]]), &((*c)[lowIdx[i]]));
    for(int i = 0; i < medNum; ++i)
        fresnelApprox(medVal[i], &((*s)[medIdx[i]]), &((*c)[medIdx[i]]));
}

#else //EIGEN_VECTORIZE_SSE

//The unvectorized version
void fresnelApprox(const VectorXd &t, VectorXd *s, VectorXd *c)
{
    s->resize(t.size());
    c->resize(t.size());
    for(int i = 0; i < t.size(); ++i)
        fresnelApprox(t[i], &((*s)[i]), &((*c)[i]));
}

#endif //EIGEN_VECTORIZE_SSE

END_NAMESPACE_Cornu

