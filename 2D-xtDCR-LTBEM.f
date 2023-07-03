      program main
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      parameter (pi=3.141592654d0)
      parameter (ns=12)
      parameter (ntime=4)
      character*20 outfil,infil,pr
      dimension nsse(10),kode(nto),xsil(10),ysil(10),xsiu(10),ysiu(10),
     &          sqlo1(nto),sqlo2(nto),sqhi1(nto),sqhi2(nto)
      dimension sqm1(nto),sqm2(nto)
      dimension sqn1(nto),sqn2(nto)
      dimension aa(nto,nto),bb(nto),xx(nto),uu(nto),pp(nto),bvalue(nto)
      dimension ttt(ntime)
      dimension xxint(19), yyint(19)
      dimension sum_uu(nto),sum_pp(nto)
      dimension sum_mu(19,19), sum_mu1(19,19), sum_mu2(19,19)
      dimension v(ns)
      common /moduli/ rk11,rk12,rk22
      common /koord/ sqlo1,sqlo2,sqhi1,sqhi2
      common /normal/ sqn1,sqn2
      common /tengah/ sqm1,sqm2
      common /constantD/ CapD
      common /vektorv/ v1, v2
      common /laplace_constant/ s
      common /alpha/ ralpha
      common /lambda/ rlambda
      common /constantKsmall/ SmallK
      data ttt / .5d0, 1.d0, 1.50, 2.00 /
      data xxint /.05d0,.1d0,.15d0,.2d0,.25d0,.3d0,.35d0,.4d0,.45d0,.5d0
     &           ,.55d0,.6d0,.65d0,.7d0,.75d0,.8d0,.85d0,.9d0,.95d0/
      data yyint /.05d0,.1d0,.15d0,.2d0,.25d0,.3d0,.35d0,.4d0,.45d0,.5d0
     &           ,.55d0,.6d0,.65d0,.7d0,.75d0,.8d0,.85d0,.9d0,.95d0/
      call cpu_time ( t1 )
c  INPUT AND OUTPUT FILES
      outfil='2D-xtDCR-LTBEM.out'
      open(unit=1,file=outfil,status='old')
      infil='2D-xtDCR-LTBEM.dat'
      OPEN(unit=5,file=infil,status='old')
      write (1,*) '# ',outfil
      write (6,*) '# ',outfil
      write (1,*) '# N =',ns
      write (6,*) '# N =',ns
      write(1,*)
      write(6,*)
c  DISCRETISATION OF THE BOUNDARY
      nn=0
      read(5,*) nsides
      do i=1,nsides
        read(5,*) nsse(i),xsil(i),ysil(i),xsiu(i),ysiu(i)
 33     format(5x,i2,10x,i3,13x,4f5.2)
        nn=nn+nsse(i)
        dx1=(xsiu(i)-xsil(i))/dfloat(nsse(i))
        dx2=(ysiu(i)-ysil(i))/dfloat(nsse(i))
        do j=1,nsse(i)
          k=k+1
          sqlo1(k)=xsil(i)+dfloat(j-1)*dx1
          sqlo2(k)=ysil(i)+dfloat(j-1)*dx2
          sqhi1(k)=xsil(i)+dfloat(j)*dx1
          sqhi2(k)=ysil(i)+dfloat(j)*dx2
          sqm1(k)=sqlo1(k)+(sqhi1(k)-sqlo1(k))/2.d0
          sqm2(k)=sqlo2(k)+(sqhi2(k)-sqlo2(k))/2.d0
          ddx=sqhi1(k)-sqlo1(k)
          ddy=sqhi2(k)-sqlo2(k)
          dds=dsqrt( ddx*ddx+ddy*ddy )
          sqn1(k)=ddy/dds
          sqn2(k)=-ddx/dsd
        enddo
      enddo
c  Calculate Stehfest coefficients v_m
      ns2=ns/2
      do m=1,ns
        if (m.le.ns2) then
          ks=m
        else
          ks=ns2
        endif
        v(m)=.0d0
        do k=int(m+1)/2,ks
          v(m)=v(m)+k**ns2*fakt(2*k)
     &      /fakt(ns2-k)/fakt(k)/fakt(k-1)/fakt(m-k)/fakt(2*k-m)
        enddo
        v(m)=v(m)*(-1)**(ns2+m)
      enddo
      trmse = 0.0d0
      trmse1 = 0.0d0
      trmse2 = 0.0d0
      tmae = 0.0d0
      tmae1 = 0.0d0
      tmae2 = 0.0d0
      tmape = 0.0d0
      tmape1 = 0.0d0
      tmape2 = 0.0d0
c  TIME t ITERATION
      do kkk=1,ntime
        tau=ttt(kkk)
        do ii=1,19
          do jj=1,19
            sum_mu(ii,jj)  = 0.0d0
            sum_mu1(ii,jj) = 0.0d0
            sum_mu2(ii,jj) = 0.0d0
          enddo
        enddo
        do k=1,nn
          sum_uu(k) = 0.0d0
          sum_pp(k) = 0.0d0
        enddo
c  STEHFEST ITERATION
        do nnnnn=1,ns
          s = real(nnnnn)*dlog(2.d0)/tau
          rk11=1.d0/(1.d0+s)
          rk12=.25d0*rk11
          rk22=.25d0*rk11
          rho = .5d0*dsqrt(pi/s**3.d0)*dexp(.25d0/s)
          v1 = 1.d0/s - .8862269255d0*dexp(-.25d0/s)/(s**1.5d0)
          v2 = 2.d0*v1
          rlambda = .0075d0/(1.d0+s) - .1d0/s
     &      + .08862269255d0*dexp(-.25d0/s)/s**1.5d0
          ralpha = -( .8862269255d0*dexp(.25d0/s)*(s+1.d0)
     &    - 0.06d0*s**1.5d0 - .1d0*s**.5d0 +
     &    .08862269255d0*dexp(-.25d0/s)*(s+1.d0)  )
     &    / (s**2.5d0*(1.d0+s))
          call findroot(rk11,rk12,rk22)
c  THE COEFFICIENT k IN THE GOVERNING EQUATION
          SmallK=rlambda+rho+s*ralpha
          print *,"n =", nnnnn
          call fundpar(rk11,rk12,rk22,SmallK)
c  THE BOUNDARY DATA
          do k=1,nn
            if (k.le.nsse(1)+nsse(2)+nsse(3)) then
              kode(k)=0
              bvalue(k)=Pexs(k,sqm1(k),sqm2(k))
            else
              kode(k)=1
              bvalue(k)=exs(sqm1(k),sqm2(k))
            endif
          enddo
c  SETUP OF THE LINEAR SYSTEM OF EQUATIONS
          do k=1,nn
            bb(k)=0.d0
            spt1=sqm1(k)
            spt2=sqm2(k)
            do l=1,nn
              call bode(99,l,spt1,spt2,smp,smg,smpv,smp1,smg1,smpv1,
     &               smp2,smg2,smpv2)
              if (kode(l).eq.0) then
                aa(k,l)=-(smpv-smg)
                bb(k)=bb(k)+smp*bvalue(l)
              else
                aa(k,l)=-smp
                bb(k)=bb(k)+(smpv-smg)*bvalue(l)
              endif
            enddo
            if (kode(k).eq.1) then
              bb(k)=bb(k)-0.5d0*bvalue(k)*h12(sqm1(k),sqm2(k))
            else
              aa(k,k)=aa(k,k)+0.5d0*h12(sqm1(k),sqm2(k))
            endif
          enddo
c  SOLVING THE LINEAR ALGEBRAIC EQUATION SYSTEM
          call laes(nn,aa,bb,xx)
cc  SOLUTIONS ON THE BOUNDARY
          do k=1,nn
            if (kode(k) .eq. 1) then
              uu(k)=bvalue(k)
              pp(k)=xx(k)
              sum_uu(k) = sum_uu(k) + ( uu(k)*v(nnnnn)*
     &                    dlog(2.0d0)/tau )
              sum_pp(k) = sum_pp(k) + ( pp(k)*v(nnnnn)*
     &                    dlog(2.0d0)/tau )
            else
              pp(k)=bvalue(k)
              uu(k)=xx(k)
              sum_uu(k) = sum_uu(k) + ( uu(k)*v(nnnnn)*
     &                    dlog(2.0d0)/tau )
              sum_pp(k) = sum_pp(k) + ( pp(k)*v(nnnnn)*
     &                    dlog(2.0d0)/tau )
            endif
          enddo
          do ii=1,19
            xint=xxint(ii)
            do jj=1,19
              yint=yyint(jj)
              psi=0.d0
              psi1=0.d0
              psi2=0.d0
              do l=1,nn
                call bode(1,l,xint,yint,smp,smg,smpv,smp1,smg1,smpv1,
     &               smp2,smg2,smpv2)
                psi  = psi + smp*pp(l) + (smpv-smg)*uu(l)
                psi1 = psi1 + smp1*pp(l) + (smpv1-smg1)*uu(l)
                psi2 = psi2 + smp2*pp(l) + (smpv2-smg2)*uu(l)
              enddo
              c=psi/h12(xint,yint)
              dcdx1=(psi1-c*dh12dx1(xint,yint))/h12(xint,yint)
              dcdx2=(psi2-c*dh12dx2(xint,yint))/h12(xint,yint)
              sum_mu(ii,jj)  = sum_mu(ii,jj) + (c *v(nnnnn)
     &                         *dlog(2.0d0)/tau)
              sum_mu1(ii,jj) = sum_mu1(ii,jj)+ (dcdx1*v(nnnnn)
     &                         *dlog(2.0d0)/tau)
              sum_mu2(ii,jj) = sum_mu2(ii,jj)+ (dcdx2*v(nnnnn)
     &                         *dlog(2.0d0)/tau)
            enddo
          enddo
        enddo
        write(6,*)
c  END OF STEHFEST ITERATION
 117    format('#',2x,186('-'))
        write(6,*) '# Solutions :'
        write(1,*) '# Solutions :'
        write(6,117)
        write(1,117)
        write(6,236)
        write(1,236)
 236    format('#',1x,'t (1)',4x,'(x,y) (2-3)',10x,
     &'BEM solutions (4-6)',10x,'Analytical solutions (7-9)',10x,
     &'Error Square (10-12)',14x,'Abs Error (13-15)',
     &    15x,'Abs Rel Error (16-18)')
        write(6,117)
        write(1,117)
        ermse = 0.0d0
        ermse1 = 0.0d0
        ermse2 = 0.0d0
        emae = 0.0d0
        emae1 = 0.0d0
        emae2 = 0.0d0
        emape = 0.0d0
        emape1 = 0.0d0
        emape2 = 0.0d0
        do ii=1,19
          xxi=xxint(ii)
          do jj=1,19
            yyi=yyint(jj)
       write(6,237)
     &   tau,xxi,yyi,sum_mu(ii,jj),sum_mu1(ii,jj),sum_mu2(ii,jj),
     &   ext(xxi,yyi,tau), ext1(xxi,yyi,tau), ext2(xxi,yyi,tau),
     &   (sum_mu(ii,jj)-ext(xxi,yyi,tau))**2.d0,
     &   (sum_mu1(ii,jj)-ext1(xxi,yyi,tau))**2.d0,
     &   (sum_mu2(ii,jj)-ext2(xxi,yyi,tau))**2.d0,
     &   dabs(sum_mu(ii,jj)-ext(xxi,yyi,tau)),
     &   dabs(sum_mu1(ii,jj)-ext1(xxi,yyi,tau)),
     &   dabs(sum_mu2(ii,jj)-ext2(xxi,yyi,tau)),
     &   dabs( (sum_mu(ii,jj)-ext(xxi,yyi,tau)) / ext(xxi,yyi,tau) ),
     &   dabs( (sum_mu1(ii,jj)-ext1(xxi,yyi,tau)) / ext1(xxi,yyi,tau) ),
     &   dabs( (sum_mu2(ii,jj)-ext2(xxi,yyi,tau)) / ext2(xxi,yyi,tau) )
       write(1,237)
     &   tau,xxi,yyi,sum_mu(ii,jj),sum_mu1(ii,jj),sum_mu2(ii,jj),
     &   ext(xxi,yyi,tau), ext1(xxi,yyi,tau), ext2(xxi,yyi,tau),
     &   (sum_mu(ii,jj)-ext(xxi,yyi,tau))**2.d0,
     &   (sum_mu1(ii,jj)-ext1(xxi,yyi,tau))**2.d0,
     &   (sum_mu2(ii,jj)-ext2(xxi,yyi,tau))**2.d0,
     &   dabs(sum_mu(ii,jj)-ext(xxi,yyi,tau)),
     &   dabs(sum_mu1(ii,jj)-ext1(xxi,yyi,tau)),
     &   dabs(sum_mu2(ii,jj)-ext2(xxi,yyi,tau)),
     &   dabs( (sum_mu(ii,jj)-ext(xxi,yyi,tau)) / ext(xxi,yyi,tau) ),
     &   dabs( (sum_mu1(ii,jj)-ext1(xxi,yyi,tau)) / ext1(xxi,yyi,tau) ),
     &   dabs( (sum_mu2(ii,jj)-ext2(xxi,yyi,tau)) / ext2(xxi,yyi,tau) )
 237    format(f7.3,1x,f7.2,1x,f7.2,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,
     &   1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,
     &   1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6)
          ermse = ermse
     &                + (sum_mu(ii,jj)-ext(xxi,yyi,tau))**2.d0
          ermse1 = ermse1
     &                 + (sum_mu1(ii,jj)-ext1(xxi,yyi,tau))**2.d0
          ermse2 = ermse2
     &                 + (sum_mu2(ii,jj)-ext2(xxi,yyi,tau))**2.d0
          emae = emae
     &               + dabs(sum_mu(ii,jj)-ext(xxi,yyi,tau))
          emae1 = emae1
     &                + dabs(sum_mu1(ii,jj)-ext1(xxi,yyi,tau))
          emae2 = emae2
     &                + dabs(sum_mu2(ii,jj)-ext2(xxi,yyi,tau))
          emape = emape
     &    + dabs( (sum_mu(ii,jj)-ext(xxi,yyi,tau)) / ext(xxi,yyi,tau) )
          emape1 = emape1
     &    + dabs( (sum_mu1(ii,jj)-ext1(xxi,yyi,tau))/ext1(xxi,yyi,tau) )
          emape2 = emape2
     &    + dabs( (sum_mu2(ii,jj)-ext2(xxi,yyi,tau))/ext2(xxi,yyi,tau) )

          trmse = trmse
     &                + (sum_mu(ii,jj)-ext(xxi,yyi,tau))**2.d0
          trmse1 = trmse1
     &                + (sum_mu1(ii,jj)-ext1(xxi,yyi,tau))**2.d0
          trmse2 = trmse2
     &                + (sum_mu2(ii,jj)-ext2(xxi,yyi,tau))**2.d0
          tmae = tmae
     &               + dabs(sum_mu(ii,jj)-ext(xxi,yyi,tau))
          tmae1 = tmae1
     &               + dabs(sum_mu1(ii,jj)-ext1(xxi,yyi,tau))
          tmae2 = tmae2
     &               + dabs(sum_mu2(ii,jj)-ext2(xxi,yyi,tau))
          tmape = tmape
     &    + dabs( (sum_mu(ii,jj)-ext(xxi,yyi,tau)) / ext(xxi,yyi,tau) )
          tmape1 = tmape1
     &    + dabs( (sum_mu1(ii,jj)-ext1(xxi,yyi,tau))/ext1(xxi,yyi,tau) )
          tmape2 = tmape2
     &    + dabs( (sum_mu2(ii,jj)-ext2(xxi,yyi,tau))/ext2(xxi,yyi,tau) )
          enddo
        enddo
        write(6,117)
        write(1,117)
      write(1,*)'# For t =', tau
      write(6,*)'# For t =', tau
      write(1,*)'# t (1)               RMSE (2-4)
     & MAE (5-7)                      MAPE (8-10)'
      write(6,*)'# t (1)               RMSE (2-4)
     & MAE (5-7)                      MAPE (8-10)'
        write(1,999) tau,dsqrt(ermse/19.d0/19.d0),
     &                 dsqrt(ermse1/19.d0/19.d0),
     &                 dsqrt(ermse2/19.d0/19.d0),
     &                 emae/19.d0/19.d0,
     &                 emae1/19.d0/19.d0,
     &                 emae2/19.d0/19.d0,
     &                 emape/19.d0/19.d0,
     &                 emape1/19.d0/19.d0,
     &                 emape2/19.d0/19.d0
        write(6,999) tau,dsqrt(ermse/19.d0/19.d0),
     &                 dsqrt(ermse1/19.d0/19.d0),
     &                 dsqrt(ermse2/19.d0/19.d0),
     &                 emae/19.d0/19.d0,
     &                 emae1/19.d0/19.d0,
     &                 emae2/19.d0/19.d0,
     &                 emape/19.d0/19.d0,
     &                 emape1/19.d0/19.d0,
     &                 emape2/19.d0/19.d0
 999    format(f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,
     &         f10.6,1x,f10.6,1x,f10.6,1x,f10.6)
        write(6,*)
        write(1,*)
      enddo
c  END OF TIME t ITERATION
      call cpu_time ( t2 )
      write(1,*)'# For all time-steps:'
      write(6,*)'# For all time-steps:'
      write(1,*) '#N(1) CPU-time(2)             RMSE(3-5)
     &         MAE(6-8)                        MAPE (9-11)'
      write(6,*) '#N(1) CPU-time(2)             RMSE(3-5)
     &         MAE(6-8)                        MAPE (9-11)'
      write(1,9999) ns, t2-t1,
     &   dsqrt(trmse/19.d0/19.d0/ntime),
     &   dsqrt(trmse1/19.d0/19.d0/ntime),
     &   dsqrt(trmse2/19.d0/19.d0/ntime),
     &   tmae/19.d0/19.d0/ntime,
     &   tmae1/19.d0/19.d0/ntime,
     &   tmae2/19.d0/19.d0/ntime,
     &   tmape/19.d0/19.d0/ntime,
     &   tmape1/19.d0/19.d0/ntime,
     &   tmape2/19.d0/19.d0/ntime
      write(6,9999) ns, t2-t1,
     &   dsqrt(trmse/19.d0/19.d0/ntime),
     &   dsqrt(trmse1/19.d0/19.d0/ntime),
     &   dsqrt(trmse2/19.d0/19.d0/ntime),
     &   tmae/19.d0/19.d0/ntime,
     &   tmae1/19.d0/19.d0/ntime,
     &   tmae2/19.d0/19.d0/ntime,
     &   tmape/19.d0/19.d0/ntime,
     &   tmape1/19.d0/19.d0/ntime,
     &   tmape2/19.d0/19.d0/ntime
9999  format(2x,i2,3x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,
     &          f10.6,1x,f10.6,1x,f10.6,1x,f10.6,1x,f10.6)
      end


c  ------------------------------------------------------------------------
      subroutine findroot(rk11,rk12,rk22)
      implicit real*8 (a-h, o-z)
      complex*8 roots
      common /akar/ rtsdot, rtsddot
      if (rk11.le.0.d0 .or. rk22.le.0.d0) then
        write(6,*) '  k11 and k22 should be both positive. '
        write(1,*) '  k11 and k22 must be both positive. '
        stop
      else
      endif
      det=dble(rk11*rk22-rk12*rk12)
      if (det.lt.0.d0) then
        write(6,*) '  k11.k22-k12.k12 =',det
        write(1,*) '  k11.k22-k12.k12 =',det
        write(6,*) '  k11.k22-k12.k12 should be positive.'
        write(1,*) '  k11.k22-k12.k12 must be positive.'
        stop
      else
      endif
      roots=cmplx(-rk12/rk22,dsqrt(det)/rk22)
      rtdsot = real(roots)
      rtsddot = aimag(roots)
      return
      end

c  ------------------------------------------------------------------------
      subroutine fundpar(rk11,rk12,rk22,SmallK)
      implicit real*8 (a-h, o-z)
      complex*8 roots
      dimension VecV(2)
      common /vektorv/ v1, v2
      common /akar/ rtsdot, rtsddot
      common /constantmu/ Cmu
      common /constantD / CapD
      common /vectorV / VecV
      common /constantCapK/ CapK
      common /laplace_constant/ s
      common /alpha/ ralpha
      common /lambda/ rlambda
      CapD = .5d0*(rk11+2.d0*rtsdot*rk12+
     &            (rtsdot**2.d0+rtsddot**2.d0)*rk22)
      CapK = rtsddot/CapD
      VecV(1) = (v1 + rtsdot*v2)
      VecV(2) = (rtsddot*v2)
      Cmu = dsqrt( (VecV(1)**2.d0+VecV(2)**2.d0)/(4.d0*CapD**2.0d0)
     &      + SmallK/CapD )
      return
      end

c  ------------------------------------------------------------------------
      subroutine bode(ints,li,spt1,spt2,smp,smg,smpv,smp1,smg1,
     &                 smpv1,smp2,smg2,smpv2)
      implicit real*8 (a-h, o-z)
      parameter (nto=1000)
      dimension qlo1(nto),qlo2(nto),qhi1(nto),qhi2(nto),
     &          qn1(nto),qn2(nto)
      dimension w(5)
      complex*8 qkm1,qk,vectx,vecty
      common /vektorv/ v1, v2
      common /koord/ qlo1,qlo2,qhi1,qhi2
      common /normal/ qn1, qn2
      external Phi,Gamma,Phi1,Gamma1,Phi2,Gamma2
      data w / 2857.0d0, 15741.0d0, 1080.0d0, 19344.0d0, 5778.0d0 /
      smp=0.d0
      smg=0.d0
      smpv=0.d0
      smp1=0.d0
      smg1=0.d0
      smpv1=0.d0
      smp2=0.d0
      smg2=0.d0
      smpv2=0.d0
      qkm1=cmplx(qlo1(li),qlo2(li))
      qk=cmplx(qhi1(li),qhi2(li))
      do kkk=1,5
        vectx=qkm1+dfloat(kkk-1)*(qk-qkm1)/9.d0
        x1=real(vectx)
        x2=aimag(vectx)
        vecty=qkm1+dfloat(10-kkk)*(qk-qkm1)/9.d0
        y1=real(vecty)
        y2=aimag(vecty)
        call Pkig(li,x1,x2,pkix)
        call Pkig(li,y1,y2,pkiy)
        smp = smp + ( Phi(x1,x2,spt1,spt2)/h12(x1,x2) +
     &                Phi(y1,y2,spt1,spt2)/h12(y1,y2) ) * w(kkk)
        smg = smg + ( Gamma(li,x1,x2,spt1,spt2)*h12(x1,x2) +
     &                Gamma(li,y1,y2,spt1,spt2)*h12(y1,y2) ) * w(kkk)
        smpv = smpv + ( (pkix - (v1*qn1(li)+v2*qn2(li))*h12(x1,x2)) *
     &                  Phi(x1,x2,spt1,spt2)
     &              +   (pkiy - (v1*qn1(li)+v2*qn2(li))*h12(y1,y2)) *
     &                  Phi(y1,y2,spt1,spt2) ) * w(kkk)
        if (ints.eq.1) then
        smp1 = smp1 + ( Phi1(x1,x2,spt1,spt2)/h12(x1,x2) +
     &                  Phi1(y1,y2,spt1,spt2)/h12(y1,y2) ) * w(kkk)
        smg1 = smg1 + ( Gamma1(li,x1,x2,spt1,spt2)*h12(x1,x2) +
     &                  Gamma1(li,y1,y2,spt1,spt2)*h12(y1,y2) ) * w(kkk)
        smpv1 = smpv1 + ( (pkix - (v1*qn1(li)+v2*qn2(li))*h12(x1,x2)) *
     &                    Phi1(x1,x2,spt1,spt2)
     &                +   (pkiy - (v1*qn1(li)+v2*qn2(li))*h12(y1,y2)) *
     &                    Phi1(y1,y2,spt1,spt2) ) * w(kkk)
        smp2 = smp2 + ( Phi2(x1,x2,spt1,spt2)/h12(x1,x2) +
     &                  Phi2(y1,y2,spt1,spt2)/h12(y1,y2) ) * w(kkk)
        smg2 = smg2 + ( Gamma2(li,x1,x2,spt1,spt2)*h12(x1,x2) +
     &                  Gamma2(li,y1,y2,spt1,spt2)*h12(y1,y2) ) * w(kkk)
        smpv2 = smpv2 + ( (pkix - (v1*qn1(li)+v2*qn2(li))*h12(x1,x2)) *
     &                    Phi2(x1,x2,spt1,spt2)
     &                +   (pkiy - (v1*qn1(li)+v2*qn2(li))*h12(y1,y2)) *
     &                    Phi2(y1,y2,spt1,spt2) ) * w(kkk)
        else
        endif
      enddo
      h=abs(qkm1-qk)
      smp=smp*h/98600.0d0
      smg=smg*h/89600.0d0
      smpv=smpv*h/89600.0d0
      if (ints.eq.1) then
        smp1=smp1*h/89600.0d0
        smg1=smg1*h/89600.0d0
        smpv1=smpv1*h/89600.0d0
        smp2=smp2*h/89600.0d0
        smg2=smg2*h/89600.0d0
        smpv2=smpv2*h/89600.0d0
      else
      endif
      return
      end

c  ------------------------------------------------------------------------
      subroutine Pkig(kp,x1,x2,Pki)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension xn(nto),yn(nto)
      common /moduli/ rk11,rk12,rk22
      common /normal/ xn,yn
      der1=dh12dx1(x1,x2)
      der2=dh12dx2(x1,x2)
      Pki = rk11*der1*xn(kp) + rk12*(der1*yn(kp)+der2*xn(kp))
     &      + rk22*der2*yn(kp)
      return
      end

c  ------------------------------------------------------------------------
      subroutine laes(nn,aa,bb,xx)
      implicit real*8 (a-h, o-z)
      parameter (nto=1000)
      dimension aa(nto,nto),bb(nto),xx(nto),wk(nto),bp(nto),ap(nto),
     &       cp(nto)
      ig=nn
      jg=nn
      itt=0
  60  do i=1,ig
        ape=0.d0
        do j=1,ig
          ape=aa(i,j)*xx(j)+ape
        enddo
        wk(i)=bb(i)-ape
      enddo
      do i=1,ig
        ape=0.d0
        do j=1,ig
          ape=ape+aa(j,i)*wk(j)
        enddo
        ap(i)=ape
        cp(i)=ape
      enddo
      irr=0
  94  itt=itt+1
      irr=irr+1
      do i=1,ig
        bp(i)=0.d0
      enddo
      do j=1,ig
        do i=1,ig
          bp(i)=aa(i,j)*ap(j)+bp(i)
        enddo
      enddo
      daa=0.d0
      dab=0.d0
      do j=1,ig
        daa=daa+cp(j)*cp(j)
        dab=dab+bp(j)*bp(j)
      enddo
      if(dab.lt.1.0e-37.and.irr.eq.1) go to 115
      if(dab.lt.1.0e-37) go to 60
      dac=daa/dab
      do i=1,ig
        wk(i)=wk(i)-dac*bp(i)
        xx(i)=xx(i)+dac*ap(i)
      enddo
      err=0.d0
      do i=1,ig
        err=err+wk(i)*wk(i)
      enddo
      if(err.lt.1.0e-37.or.itt.eq.jg) go to 115
      if(irr.eq.ig) go to 60
      bong=0.d0
      do i=1,ig
        ape=0.d0
        do j=1,ig
          ape=ape+aa(j,i)*wk(j)
        enddo
        cp(i)=ape
        bong=bong+ape*ape
      enddo
      bong=bong/daa
      do i=1,ig
        ap(i)=cp(i)+bong*ap(i)
      enddo
      go to 94
 115  return
      end

c  ------------------------------------------------------------------------
      function Phi(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      dimension VecV(2), VecR(2)
      common /akar/ rtsdot, rtsddot
      common /constantmu/ Cmu
      common /constantD / CapD
      common /vectorV / VecV
      common /constantCapK/ CapK
      pi2 = 8.d0*datan(1.d0)
      r = Rcap(x1,x2,xi1,xi2)
      VecR(1) = (x1-xi1) + (x2-xi2)*rtsdot
      VecR(2) = (x2-xi2)*rtsddot
      ProductVR = VecV(1)*VecR(1) + VecV(2)*VecR(2)
      x = Cmu*r
      FFK0 = FK0(x)
      Phi = CapK/pi2*dexp(-ProductVR/(2.d0*CapD))*FFK0
      return
      end

      function Phi1(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      dimension VecV(2)
      common /constantmu/ Cmu
      common /constantD / CapD
      common /vectorV / VecV
      vx = VecV(1)
      PPhi = Phi(x1,x2,xi1,xi2)
      dR1 = dRcap1(x1,x2,xi1,xi2)
      r = Rcap(x1,x2,xi1,xi2)
      x = Cmu*r
      FFK0 = FK0(x)
      FFK1 = FK1(x)
      Phi1 = (vx/(2.d0*CapD) - Cmu*dR1*FFK1/FFK0) * PPhi
      return
      end

      function Phi2(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      dimension VecV(2)
      common /constantmu/ Cmu
      common /constantD / CapD
      common /vectorV / VecV
      common /constantCapK/ CapK
      common /akar/ rtsdot, rtsddot
      vx = VecV(1)
      vy = VecV(2)
      PPhi = Phi(x1,x2,xi1,xi2)
      dR2 = dRcap2(x1,x2,xi1,xi2)
      r = Rcap(x1,x2,xi1,xi2)
      x = Cmu*r
      FFK0 = FK0(x)
      FFK1 = FK1(x)
      Phi2 = ((vx*rtsdot+vy*rtsddot)/(2.d0*CapD) - Cmu*dR2*FFK1/
     &         FFK0) * PPhi
      return
      end

      function Phi11(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      dimension VecV(2)
      common /constantmu/ Cmu
      common /constantD / CapD
      common /vectorV / VecV
      common /constantCapK/ CapK
      vx = VecV(1)
      PPhi = Phi(x1,x2,xi1,xi2)
      PPhi1 = Phi1(x1,x2,xi1,xi2)
      dR1 = dRcap1(x1,x2,xi1,xi2)
      dR11 = dRcap11(x1,x2,xi1,xi2)
      r = Rcap(x1,x2,xi1,xi2)
      x = Cmu*r
      FFK0 = FK0(x)
      FFK1 = FK1(x)
      Phi11 = PPhi1*(vx/(2.d0*CapD)-Cmu*FFK1/FFK0*dR1)
     &   -PPhi*Cmu*(Cmu*dR1**2.d0*(-1.d0-FFK1/FFK0/Cmu/r+
     &   FFK1**2.d0/FFK0**2.d0)+FFK1/FFK0*dR11)
      return
      end

      function Phi12(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      dimension VecV(2)
      common /constantmu/ Cmu
      common /constantD / CapD
      common /vectorV / VecV
      common /constantCapK/ CapK
      vx = VecV(1)
      PPhi = Phi(x1,x2,xi1,xi2)
      PPhi2 = Phi2(x1,x2,xi1,xi2)
      dR1 = dRcap1(x1,x2,xi1,xi2)
      dR2 = dRcap2(x1,x2,xi1,xi2)
      dR12 = dRcap12(x1,x2,xi1,xi2)
      r = Rcap(x1,x2,xi1,xi2)
      x = Cmu*r
      FFK0 = FK0(x)
      FFK1 = FK1(x)
      Phi12 = PPhi2*(vx/(2.d0*CapD)-Cmu*FFK1/FFK0*dR1)
     &   -PPhi*Cmu*(Cmu*dR2*dR2*(-1.d0-FFK1/FFK0/Cmu/r+
     &   FFK1**2.d0/FFK0**2.d0)+FFK1/FFK0*dR12)
      return
      end

      function Phi22(x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      dimension VecV(2)
      common /akar/ rtsdot, rtsddot
      common /constantmu/ Cmu
      common /constantD / CapD
      common /vectorV / VecV
      common /constantCapK/ CapK
      vx = VecV(1)
      vy = VecV(2)
      PPhi = Phi(x1,x2,xi1,xi2)
      PPhi2 = Phi2(x1,x2,xi1,xi2)
      dR2 = dRcap2(x1,x2,xi1,xi2)
      dR22 = dRcap22(x1,x2,xi1,xi2)
      r = Rcap(x1,x2,xi1,xi2)
      x = Cmu*r
      FFK0 = FK0(x)
      FFK1 = FK1(x)
      Phi22 = PPhi2*((vx*rtsdot+vy*rtsddot)/(2.d0*CapD)-Cmu*
     &   FFK1/FFK0*dR2)-PPhi*Cmu*(Cmu*dR2**2.d0*(-1.d0-FFK1/
     &   FFK0/Cmu/r+FFK1**2.d0/FFK0**2.d0)+FFK1/FFK0*dR22)
      return
      end

c  ------------------------------------------------------------------------
      function Gamma(ki,x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension qn1(nto),qn2(nto)
      common /moduli/ rk11,rk12,rk22
      common /normal/ qn1,qn2
      dPhi1 = -Phi1(x1,x2,xi1,xi2)
      dPhi2 = -Phi2(x1,x2,xi1,xi2)
      Gamma = rk11*dPhi1*qn1(ki) + rk12*(dPhi1*qn2(ki)+dPhi2*qn1(ki))
     &        + rk22*dPhi2*qn2(ki)
      return
      end

      function Gamma1(ki,x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension qn1(nto),qn2(nto)
      common /moduli/ rk11,rk12,rk22
      common /normal/ qn1,qn2
c      external Phi11,Phi12
      dPhi11 = -Phi11(x1,x2,xi1,xi2)
      dPhi12 = -Phi12(x1,x2,xi1,xi2)
      Gamma1 = rk11*dPhi11*qn1(ki) + rk12*(dPhi11*qn2(ki)+
     &         dPhi12*qn1(ki)) + rk22*dPhi12*qn2(ki)
      return
      end

      function Gamma2(ki,x1,x2,xi1,xi2)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension qn1(nto),qn2(nto)
      common /moduli/ rk11,rk12,rk22
      common /normal/ qn1,qn2
      dPhi12 = -Phi12(x1,x2,xi1,xi2)
      dPhi22 = -Phi22(x1,x2,xi1,xi2)
      Gamma2 = rk11*dPhi12*qn1(ki) + rk12*(dPhi12*qn2(ki)+
     &         dPhi22*qn1(ki)) + rk22*dPhi22*qn2(ki)
      return
      end

c  ------------------------------------------------------------------------
      function FK0(x)
      implicit real*8 (a-h,o-z)
      DATA p1,p2,p3,p4,p5,p6,p7/-0.57721566d0,0.42278420d0,0.23069756d0,
     &  0.3488590d-1,0.262698d-2,0.10750d-3,0.74d-5/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,-0.7832358d-1,0.2189568d-1,
     &  -0.1062446d-1,0.587872d-2,-0.251540d-2,0.53208d-3/
      if (x.le.2.0d0) then
        y=x*x/4.0d0
        FK0=(-dlog(x/2.0d0)*FI0(x))+(p1+y*(p2+y*(p3+
     &  y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        y=(2.0d0/x)
        FK0=(dexp(-x)/dsqrt(x))*(q1+y*(q2+y*(q3+
     &  y*(q4+y*(q5+y*(q6+y*q7))))))
      endif
      return
      end

      function FK1(x)
      implicit real*8 (a-h,o-z)
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,0.15443144d0,-0.67278579d0,
     &  -0.18156897d0,-0.1919402d-1,-0.110404d-2,-0.4686d-4/
      DATA q1,q2,q3,q4,q5,q6,q7/1.25331414d0,0.23498619d0,-0.3655620d-1,
     &  0.1504268d-1,-0.780353d-2,0.325614d-2,-0.68245d-3/
      if (x.le.2.0d0) then
        y=x*x/4.0d0
        FK1=(dlog(x/2.0d0)*FI1(x))+(1.0d0/x)*(p1+y*(p2+
     &  y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        y=2.0d0/x
        FK1=(dexp(-x)/dsqrt(x))*(q1+y*(q2+y*(q3+
     &  y*(q4+y*(q5+y*(q6+y*q7))))))
      endif
      return
      end

c  ------------------------------------------------------------------------
      function FI0(x)
      implicit real*8 (a-h,o-z)
      DATA p1,p2,p3,p4,p5,p6,p7/1.0d0,3.5156229d0,3.0899424d0,
     &  1.2067492d0,0.2659732d0,0.360768d-1,0.45813d-2/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,0.1328592d-1,
     &  0.225319d-2,-0.157565d-2,0.916281d-2,-0.2057706d-1,
     &  0.2635537d-1,-0.1647633d-1,0.392377d-2/
      if (dabs(x).lt.3.75d0) then
        y=(x/3.75d0)**2.d0
        FI0=p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7)))))
      else
        ax=dabs(x)
        y=3.75d0/ax
        FI0=(dexp(ax)/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4
     &    +y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
      endif
      return
      end

      function FI1(x)
      implicit real*8 (a-h,o-z)
      DATA p1,p2,p3,p4,p5,p6,p7/0.5d0,0.87890594d0,0.51498869d0,
     &  0.15084934d0,0.2658733d-1,0.301532d-2,0.32411d-3/
      DATA q1,q2,q3,q4,q5,q6,q7,q8,q9/0.39894228d0,-0.3988024d-1,
     &  -0.362018d-2,0.163801d-2,-0.1031555d-1,0.2282967d-1,
     &  -0.2895312d-1,0.1787654d-1,-0.420059d-2/
      if (dabs(x).lt.3.75d0) then
        y=(x/3.75d0)**2.d0
        FI1=x*(p1+y*(p2+y*(p3+y*(p4+y*(p5+y*(p6+y*p7))))))
      else
        ax=dabs(x)
        y=3.75d0/ax
        FI1=(dexp(ax)/dsqrt(ax))*(q1+y*(q2+y*(q3+y*(q4+
     &  y*(q5+y*(q6+y*(q7+y*(q8+y*q9))))))))
        if (x.lt.0.0d0) FI1=-FI1
      endif
      return
      end

c  ------------------------------------------------------------------------
      function Rcap(x1,x2,a,b)
      implicit real*8 (a-h,o-z)
      common /akar/ rtsdot, rtsddot
      Rcap = dsqrt( (x1+rtsdot*x2-a-rtsdot*b)**2.d0 +
     &             (x2*rtsddot-b*rtsddot)**2.d0 )
      return
      end

      function dRcap1(x1,x2,a,b)
      implicit real*8 (a-h,o-z)
      common /akar/ rtsdot, rtsddot
      RR=Rcap(x1,x2,a,b)
      if (dabs(x1-a).le.1.e-10 .and. dabs(x2-b).le.1.e-10) then
        dRcap1 = 0.0d0
      else
        dRcap1 = -(x1+rtsdot*x2-a-rtsdot*b)/RR
      endif
      return
      end

      function dRcap2(x1,x2,a,b)
      implicit real*8 (a-h,o-z)
      common /akar/ rtsdot, rtsddot
      RR=Rcap(x1,x2,a,b)
      if (dabs(x1-a).le.1.e-10 .and. dabs(x2-b).le.1.e-10 .and.
     &    dabs(rtsdot).lt.1.0d0 .and. dabs(rtsddot).lt.1.0d0) then
        dRcap2 = 0.0d0
      else
        dRcap2 = -( rtsdot*(x1+rtsdot*x2-a-rtsdot*b)+rtsddot*
     &           (x2*rtsddot-b*rtsddot) ) / RR
      endif
      return
      end

      function dRcap11(x1,x2,a,b)
      implicit real*8 (a-h,o-z)
      RR=Rcap(x1,x2,a,b)
      RR1=dRcap1(x1,x2,a,b)
      dRcap11 = (1.d0-RR1**2.d0)/RR
      return
      end

      function dRcap12(x1,x2,a,b)
      implicit real*8 (a-h,o-z)
      common /akar/ rtsdot, rtsddot
      RR=Rcap(x1,x2,a,b)
      RR1=dRcap1(x1,x2,a,b)
      RR2=dRcap2(x1,x2,a,b)
      dRcap12 = (rtsdot-RR1*RR2)/RR
      return
      end

      function dRcap22(x1,x2,a,b)
      implicit real*8 (a-h,o-z)
      common /akar/ rtsdot, rtsddot
      RR=Rcap(x1,x2,a,b)
      RR2=dRcap2(x1,x2,a,b)
      dRcap22 = (rtsdot**2.d0+rtsddot**2.d0 - RR2**2.d0)/RR
      return
      end

c  ------------------------------------------------------------------------
      function fakt(n)
      implicit real*8 (a-h, o-z)
      fakt=1.d0
      do k=1,n
        fakt=fakt*real(k)
      enddo
      return
      end

c  ------------------------------------------------------------------------
      function h12(x,y)
      implicit real*8 (a-h, o-z)
      h12=dexp(.1d0*x-.1d0*y)
      return
      end

      function dh12dx1(x,y)
      implicit real*8 (a-h, o-z)
      dh12dx1=.1d0*dexp(.1d0*x-.1d0*y)
      return
      end

      function dh12dx2(x,y)
      implicit real*8 (a-h, o-z)
      dh12dx2=-.1d0*dexp(.1d0*x-.1d0*y)
      return
      end

c  ------------------------------------------------------------------------
      function ext(x,y,t)
      implicit real*8 (a-h, o-z)
      xxp = dcos(-.2d0*x+.1d0*y)+dsin(-.2d0*x+.1d0*y)
      ext=10.d0*xxp/h12(x,y)*t**2.d0
      return
      end

      function ext1(x,y,t)
      implicit real*8 (a-h, o-z)
      xxp = dcos(-.2d0*x+.1d0*y)+dsin(-.2d0*x+.1d0*y)
      xxm = dcos(-.2d0*x+.1d0*y)-dsin(-.2d0*x+.1d0*y)
      ext1=(-2.d0*xxm*h12(x,y)-10.d0*xxp*dh12dx1(x,y))
     &     /h12(x,y)**2.d0*t**2.d0
      return
      end

      function ext2(x,y,t)
      implicit real*8 (a-h, o-z)
      xxp = dcos(-.2d0*x+.1d0*y)+dsin(-.2d0*x+.1d0*y)
      xxm = dcos(-.2d0*x+.1d0*y)-dsin(-.2d0*x+.1d0*y)
      ext2=( xxm*h12(x,y)-10.d0*xxp*dh12dx2(x,y) )
     & /h12(x,y)**2.d0*t**2.d0
      return
      end

      function exs(x,y)
      implicit real*8 (a-h,o-z)
      common /laplace_constant/ rs
      xxp = dcos(-.2d0*x+.1d0*y)+dsin(-.2d0*x+.1d0*y)
      exs = 20.d0*xxp/h12(x,y)/rs**3.d0
      return
      end

      function exs1(x,y)
      implicit real*8 (a-h,o-z)
      common /laplace_constant/ rs
      xxp = dcos(-.2d0*x+.1d0*y)+dsin(-.2d0*x+.1d0*y)
      xxm = dcos(-.2d0*x+.1d0*y)-dsin(-.2d0*x+.1d0*y)
      exs1 = 20.d0*( -.2d0*xxm*h12(x,y) - xxp*dh12dx1(x,y) )
     &       / h12(x,y)**2.d0/rs**3.d0
      return
      end

      function exs2(x,y)
      implicit real*8 (a-h,o-z)
      common /laplace_constant/ rs
      xxp = dcos(-.2d0*x+.1d0*y)+dsin(-.2d0*x+.1d0*y)
      xxm = dcos(-.2d0*x+.1d0*y)-dsin(-.2d0*x+.1d0*y)
      exs2 = 20.d0*( .1d0*(xxm)*h12(x,y) - xxp*dh12dx2(x,y) )
     &       / h12(x,y)**2.d0/rs**3.d0
      return
      end

      function Pexs(ki,x,y)
      implicit real*8 (a-h,o-z)
      common /laplace_constant/ rs
      call Pkig(ki,x,y,Pki)
      call Ppsi(ki,x,y,Ppsiex)
      xxp = dcos(-.2d0*x+.1d0*y)+dsin(-.2d0*x+.1d0*y)
      psiex = 20.d0*xxp/rs**3.d0
      Pexs = -Pki*psiex + Ppsiex*h12(x,y)
      return
      end

      subroutine Ppsi(kp,x,y,Ppsiex)
      implicit real*8 (a-h,o-z)
      parameter (nto=1000)
      dimension xn(nto),yn(nto)
      common /laplace_constant/ rs
      common /moduli/ rk11,rk12,rk22
      common /normal/ xn,yn
      xxm = dcos(-.2d0*x+.1d0*y)-dsin(-.2d0*x+.1d0*y)
      der1 = -4.d0*xxm/rs**3.d0
      der2 = 2.d0*xxm/rs**3.d0
      Ppsiex = rk11*der1*xn(kp) + rk12*(der1*yn(kp)+der2*xn(kp))
     &         + rk22*der2*yn(kp)
      return
      end
c  ------------------------------------------------------------------------

