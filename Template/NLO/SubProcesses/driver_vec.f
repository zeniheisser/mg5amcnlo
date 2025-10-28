C**************************************************************************
C     Wrappers for vectorisation routines in NLO driver
C**************************************************************************

      subroutine generate_momenta_vec(iconfig,sum,proc_map,x,rwgt,vol1,vector_size)
        use driver
        ! use mint_module
        implicit none
C Included files for process information
        include 'nexternal.inc'
        include 'nFKSconfigs.inc'
        include 'run.inc'
        include 'orders.inc'
        include 'fks_info.inc'
        include 'genps.inc'
        include 'born_nhel.inc'
C Common blocks for sigint variables not passed as arguments
        logical       nbody
        common/cnbody/nbody
        integer         nndim
        common/tosigint/nndim
        character*4      abrv
        common /to_abrv/ abrv
        double precision p1_cnt(0:3,nexternal,-2:2),wgt_cnt(-2:2)
        double precision pswgt_cnt(-2:2),jac_cnt(-2:2)
        common/counterevnts/p1_cnt,wgt_cnt,pswgt_cnt,jac_cnt
        integer              nFKSprocess
        common/c_nFKSprocess/nFKSprocess
        integer             ini_fin_fks
        common/fks_channels/ini_fin_fks
        integer icolup_s(2,nexternal-1),icolup_h(2,nexternal)
        common /colour_connections/ icolup_s,icolup_h
        logical calculatedBorn
        common/ccalculatedBorn/calculatedBorn
C Momenta
        double precision p_born(0:3,nexternal-1), p_born_rot(0:3,nexternal-1)
        common /pborn/   p_born
        double precision p_born_coll(0:3,nexternal-1)
        common/pborn_coll/p_born_coll
        double precision p_born_ev(0:3,nexternal-1)
        common/pborn_ev/ p_born_ev
        double precision p_born_norad(0:3,nexternal-1)
        common/pborn_norad/p_born_norad
        double precision p(0:3,nexternal)
C Arguments
        integer proc_map(0:fks_configs,0:fks_configs)
        double precision x(99), rwgt, vol1
        integer vector_size, sum, iconfig
C Local variables
        integer i,ivec,iFKS,k,icoup,curr_ind,nFKS_picked_nbody,nFKS_in,nFKS_out
        double precision dummy_jac
        logical passcuts,passcuts_nbody,passcuts_n1body
C Local parameters
        integer izero,ione,itwo,mohdr
        parameter (izero=0,ione=1,itwo=2,mohdr=-100)

        ! if (ini_fin_fks.eq.0) then
        !     dummy_jac=1d0
        !  else
        !     dummy_jac=0.5d0
        !  endif
          
C ZW: Generate momenta and running couplings
         nFKS_picked_nbody=proc_map(proc_map(0,1),1)
         if (sum.eq.0) then
c For sum=0, determine nFKSprocess so that the soft limit gives a non-zero Born
            nFKS_in=nFKS_picked_nbody
            call get_born_nFKSprocess(nFKS_in,nFKS_out)
            nFKS_picked_nbody=nFKS_out
         endif
         nbody=.true.
         call update_fks_dir(nFKS_picked_nbody)
         dummy_jac=1d0
         call generate_momenta(nndim,iconfig,dummy_jac,x,p)
        !  if (p_born(0,1).lt.0d0) goto 12
         call set_alphaS(p1_cnt(0,1,0))
         call set_alphaS_vec(p1_cnt(0,1,0),4*FKS_configs +1)
         spb(:,:,0,1) = p_born(:,:)
         sp1_cnt(:,:,0,1) = p1_cnt(:,:,0)
         nbody=.false.

         do i=1,proc_map(proc_map(0,1),0)
            iFKS=proc_map(proc_map(0,1),i)
            call update_fks_dir(iFKS)
            dummy_jac=1d0
            ! icolup_s(1,1)=-1    ! set colour connection to -1: i.e., complete_xmcsubt has not been called
            call generate_momenta(nndim,iconfig,dummy_jac,x,p)
            if (p_born(0,1).lt.0d0) cycle
            ! call set_cms_stuff(izero)
            ! if (ickkw.eq.3) call set_FxFx_scale(-2,p1_cnt(0,1,0),nFKSprocess)
            passcuts_nbody=passcuts(p1_cnt(0,1,0),rwgt)
            ! call set_cms_stuff(mohdr)
            ! if (ickkw.eq.3) call set_FxFx_scale(-3,p,nFKSprocess)
            passcuts_n1body=passcuts(p,rwgt)
            icoup = 4*(FKS_configs - 1) + 1
            call set_alphaS(p_born_ev)
            call set_alphaS_vec(p_born_ev,icoup)
            spb_ev(:,:,iFKS,1) = p_born_ev(:,:)
            call set_alphaS(p_born_norad)
            call set_alphaS_vec(p_born_norad,icoup + 1)
            spb_norad(:,:,iFKS,1) = p_born_norad(:,:)
            call set_alphaS(p1_cnt(0,1,0))
            call set_alphaS_vec(p1_cnt(0,1,0),icoup + 2)
            spb(:,:,iFKS,1) = p_born(:,:)
            do k=1,nexternal-1
               p_born_rot(0,k)=p_born(0,k)
               p_born_rot(1,k)=-p_born(1,k)
               p_born_rot(2,k)=p_born(2,k)
               p_born_rot(3,k)=-p_born(3,k)
            enddo
            spb_rot(:,:,iFKS,1) = p_born_rot(:,:)
            sp1_cnt(:,:,iFKS,1) = p1_cnt(:,:,0)
            spb_coll(:,:,iFKS,1) = p_born_coll(:,:)
            call set_alphaS(p)
            call set_alphaS_vec(p,icoup + 3)
            sp1(:,:,iFKS,1) = p(:,:)
         enddo
        return
      end


      subroutine amplitudes_vec(proc_map,rwgt,vector_size,skip_iter)
        use driver
        implicit none
C Included files for process information
        include 'nexternal.inc'
        include 'nFKSconfigs.inc'
        include 'run.inc'
        include 'orders.inc'
        include 'fks_info.inc'
        include 'genps.inc'
        include 'born_nhel.inc'
      logical calculatedBorn
      common/ccalculatedBorn/calculatedBorn
      logical       nbody
      common/cnbody/nbody
C Momenta
        double precision p1_cnt(0:3,nexternal,-2:2),wgt_cnt(-2:2)
        double precision pswgt_cnt(-2:2),jac_cnt(-2:2)
        common/counterevnts/p1_cnt,wgt_cnt,pswgt_cnt,jac_cnt
        double precision p_born(0:3,nexternal-1), p_born_rot(0:3,nexternal-1)
        common /pborn/   p_born
        double precision p_born_coll(0:3,nexternal-1)
        common/pborn_coll/p_born_coll
        double precision p_born_ev(0:3,nexternal-1)
        common/pborn_ev/ p_born_ev
        double precision p_born_norad(0:3,nexternal-1)
        common/pborn_norad/p_born_norad
        double precision p(0:3,nexternal)
C Local variables
        integer i,ivec,iFKS,k,icoup,curr_ind,nFKS_picked_nbody,nFKS_in,nFKS_out
        double precision dummy_jac
        logical passcuts,passcuts_nbody,passcuts_n1body
C Born variables
      double precision born_amp2(ngraphs), born_jamp2(0:ncolor)
      complex*16 born_ans_cnt(2,nsplitorders)
      double precision born_amp_split(amp_split_size)
      double complex born_amp_split_cnt(amp_split_size,2,nsplitorders)
      double complex born_saveamp(ngraphs,max_bhel)
      double precision wgt_born
      double precision ev_amp2(ngraphs), ev_jamp2(0:ncolor)
      complex*16 ev_ans_cnt(2,nsplitorders)
      double precision ev_amp_split(amp_split_size)
      double complex ev_amp_split_cnt(amp_split_size,2,nsplitorders)
      double complex ev_saveamp(ngraphs,max_bhel)
      double precision wgt_ev
      double precision norad_amp2(ngraphs), norad_jamp2(0:ncolor)
      complex*16 norad_ans_cnt(2,nsplitorders)
      double precision norad_amp_split(amp_split_size)
      double complex norad_amp_split_cnt(amp_split_size,2,nsplitorders)
      double complex norad_saveamp(ngraphs,max_bhel)
      double precision wgt_norad
      double precision coll_amp2(ngraphs), coll_jamp2(0:ncolor)
      complex*16 coll_ans_cnt(2,nsplitorders)
      double precision coll_amp_split(amp_split_size)
      double complex coll_amp_split_cnt(amp_split_size,2,nsplitorders)
      double complex coll_saveamp(ngraphs,max_bhel)
      double precision wgt_coll
C Virtual variables
      double precision amp_split_virt(amp_split_size),
     &     amp_split_born_for_virt(amp_split_size),
     &     amp_split_avv(amp_split_size)
      double precision amp_split_wgtnstmp(amp_split_size),
     $                 amp_split_wgtwnstmpmuf(amp_split_size),
     $                 amp_split_wgtwnstmpmur(amp_split_size)
      double precision bsv_wgt,virt_wgt,born_wgt
C Real variables
      double precision zero, one
      parameter (zero=0d0,one=1d0)
      double precision real_amp_split(amp_split_size)
      double precision fx_ev
C n+1-kinematic borns
      double precision n1_amp2(ngraphs), n1_jamp2(0:ncolor)
      double precision rot_amp2(ngraphs), rot_jamp2(0:ncolor)
      complex*16 n1_ans_cnt(2,nsplitorders)
      complex*16 rot_ans_cnt(2,nsplitorders)
      double precision n1_amp_split(amp_split_size)
      double precision coll_n1_amp_split(amp_split_size)
      double precision rot_amp_split(amp_split_size)
      double complex n1_amp_split_cnt(amp_split_size,2,nsplitorders)
      double complex rot_amp_split_cnt(amp_split_size,2,nsplitorders)
      double complex n1_saveamp(ngraphs,max_bhel)
      double complex rot_saveamp(ngraphs,max_bhel)
      double precision coll_n1_amp2(ngraphs), coll_n1_jamp2(0:ncolor)
      complex*16 coll_n1_cnt(2,nsplitorders)
      double complex coll_n1_split_cnt(amp_split_size,2,nsplitorders)
      double complex coll_n1_saveamp(ngraphs,max_bhel)
      double precision wgt_n1, wgt_rot, wgt_coll_n1
C Born-like variables for FKS sector dependent contributions
      double precision nb_amp2(ngraphs), nb_jamp2(0:ncolor)
      complex*16 nb_ans_cnt(2,nsplitorders)
      double precision nb_amp_split(amp_split_size)
      double complex nb_amp_split_cnt(amp_split_size,2,nsplitorders)
      double complex nb_saveamp(ngraphs,max_bhel)
      double precision wgt_nb
C Arguments
        integer proc_map(0:fks_configs,0:fks_configs)
        double precision  rwgt
        integer vector_size, sum, iconfig
        logical skip_iter

        skip_iter = .false.

        do i=1,proc_map(proc_map(0,1),0)
            iFKS=proc_map(proc_map(0,1),i)
            call update_fks_dir(iFKS)
            p_born(:,:)=spb(:,:,iFKS,1)
            p_born_ev(:,:)=spb_ev(:,:,iFKS,1)
            p_born_norad(:,:)=spb_norad(:,:,iFKS,1)
            p_born_coll(:,:)=spb_coll(:,:,iFKS,1)
            p_born_rot(:,:)=spb_rot(:,:,iFKS,1)
            p1_cnt(:,:,0)=sp1_cnt(:,:,iFKS,1)
            p(:,:) =sp1(:,:,iFKS,1)
            if (p_born(0,1).lt.0d0) cycle
            passcuts_nbody=passcuts(p1_cnt(0,1,0),rwgt)
            passcuts_n1body=passcuts(p,rwgt)
            if (.not. (passcuts_nbody.or.passcuts_n1body)) cycle
            ! if (passcuts_nbody) then
               calculatedBorn=.false.
               call sborn_amp_vec(p_born_ev,ev_amp2,ev_jamp2,ev_amp_split
     $                    ,ev_amp_split_cnt,wgt_ev,ev_ans_cnt,ev_saveamp
     $                    ,4*(FKS_configs - 1)+1)
               sev_amp2(:,iFKS,1)=ev_amp2(:)
               calculatedBorn=.false.
               call sborn_amp_vec(p_born_norad,norad_amp2,norad_jamp2,norad_amp_split
     $                    ,norad_amp_split_cnt,wgt_norad,norad_ans_cnt,norad_saveamp
     $                    ,4*(FKS_configs - 1) + 2)
               snorad_amp2(:,iFKS,1)=norad_amp2(:)
               calculatedBorn=.false.
               call sborn_amp_vec(p_born_coll,coll_amp2,coll_jamp2,coll_amp_split
     $                    ,coll_amp_split_cnt,wgt_coll,coll_ans_cnt,coll_saveamp
     $                    ,4*(FKS_configs - 1) + 3)
               scb_amp_split(:,iFKS,1)=coll_amp_split(:)
               scb_amp_split_cnt(:,:,:,iFKS,1)=coll_amp_split_cnt(:,:,:)
               scb_ans_cnt(:,:,iFKS,1)=coll_ans_cnt(:,:)
               scb_saveamp(:,:,iFKS,1)=coll_saveamp(:,:)
               calculatedBorn=.false.
               call sborn_amp_vec(p_born,nb_amp2,nb_jamp2,nb_amp_split
     $                    ,nb_amp_split_cnt,wgt_nb,nb_ans_cnt,nb_saveamp
     $                    ,4*(FKS_configs - 1) + 3)
               snb_amp2(:,iFKS,1)=nb_amp2(:)
               snb_amp_split(:,iFKS,1)=nb_amp_split(:)
               snb_amp_split_cnt(:,:,:,iFKS,1)=nb_amp_split_cnt(:,:,:)
               snb_ans_cnt(:,:,iFKS,1)=nb_ans_cnt(:,:)
               snb_saveamp(:,:,iFKS,1)=nb_saveamp(:,:)
               calculatedBorn=.false.
               call sborn_amp_vec(p_born_rot,rot_amp2,rot_jamp2,rot_amp_split
     $                    ,rot_amp_split_cnt,wgt_rot,rot_ans_cnt,rot_saveamp
     $                    ,4*(FKS_configs - 1) + 4)
               srot_jamp2(:,iFKS,1)=rot_jamp2(:)
               ! srot_amp_split(:,iFKS,1)=rot_amp_split(:)
               srot_amp_split_cnt(:,:,:,iFKS,1)=rot_amp_split_cnt(:,:,:)
               srot_ans_cnt(:,:,iFKS,1)=rot_ans_cnt(:,:)
               ! srot_saveamp(:,:,iFKS,1)=rot_saveamp(:,:)
               calculatedBorn=.false.
               call sborn_amp_vec(p_born_coll,coll_n1_amp2,coll_n1_jamp2,coll_n1_amp_split
     $                    ,coll_n1_split_cnt,wgt_coll_n1,coll_n1_cnt,coll_n1_saveamp
     $                    ,4*(FKS_configs - 1) + 4)
               sc1_amp_split(:,iFKS,1)=coll_n1_amp_split(:)
               sc1_amp_split_cnt(:,:,:,iFKS,1)=coll_n1_split_cnt(:,:,:)
               sc1_ans_cnt(:,:,iFKS,1)=coll_n1_cnt(:,:)
               sc1_saveamp(:,:,iFKS,1)=coll_n1_saveamp(:,:)
               calculatedBorn=.false.
               call sborn_amp_vec(p_born,n1_amp2,n1_jamp2,n1_amp_split
     $                    ,n1_amp_split_cnt,wgt_n1,n1_ans_cnt,n1_saveamp
     $                    ,4*(FKS_configs - 1) + 4)
               sn1_amp2(:,iFKS,1)=n1_amp2(:)
               sn1_jamp2(:,iFKS,1)=n1_jamp2(:)
               sn1_amp_split(:,iFKS,1)=n1_amp_split(:)
               sn1_amp_split_cnt(:,:,:,iFKS,1)=n1_amp_split_cnt(:,:,:)
               sn1_ans_cnt(:,:,iFKS,1)=n1_ans_cnt(:,:)
               sn1_saveamp(:,:,iFKS,1)=n1_saveamp(:,:)
            ! endif
            ! if (passcuts_n1body) then
               call smatrix_real_vec(p,real_amp_split,fx_ev,4*(FKS_configs - 1) + 4)
               sreal_amp_split(:,iFKS,1)=real_amp_split(:)
               sfx_ev(iFKS,1)=fx_ev
            ! endif
         enddo

         
         
         nbody=.true.
!          calculatedBorn=.false.
! c Pick the first one because that's the one with the soft singularity
         nFKS_picked_nbody=proc_map(proc_map(0,1),1)
         if (sum.eq.0) then
! c For sum=0, determine nFKSprocess so that the soft limit gives a non-zero Born
            nFKS_in=nFKS_picked_nbody
            call get_born_nFKSprocess(nFKS_in,nFKS_out)
            nFKS_picked_nbody=nFKS_out
         endif
         call update_fks_dir(nFKS_picked_nbody)
c Pick the first one because that's the one with the soft singularity
         p_born(:,:) = spb(:,:,0,1)
         p1_cnt(:,:,0) = sp1_cnt(:,:,0,1)
         if (p_born(0,1).lt.0d0) skip_iter = .true.
            if (skip_iter) return
         ! p_born_nb(:,:) = p_born(:,:)
         passcuts_nbody=passcuts(p1_cnt(0,1,0),rwgt)
         calculatedBorn=.false.
         call sborn_amp_vec(p_born,born_amp2,born_jamp2,born_amp_split
     $                     ,born_amp_split_cnt,wgt_born,born_ans_cnt,born_saveamp
     $                     ,4*FKS_configs + 1)

         sborn_amp2(:,1)=born_amp2(:)
         sborn_jamp2(:,1)=born_jamp2(:)
         sborn_amp_split(:,1)=born_amp_split(:)
         sborn_amp_split_cnt(:,:,:,1)=born_amp_split_cnt(:,:,:)
         sborn_ans_cnt(:,:,1)=born_ans_cnt(:,:)
         sborn_saveamp(:,:,1)=born_saveamp(:,:)
         swgt_born(1)=wgt_born
12       continue

         return
         end


      
      subroutine update_fks_dir(nFKS)
      implicit none
      include 'run.inc'
      integer nFKS
      integer              nFKSprocess
      common/c_nFKSprocess/nFKSprocess
      nFKSprocess=nFKS
      call fks_inc_chooser()
      call leshouche_inc_chooser()
      call setcuts
      call setfksfactor(.true.)
      return
      end

      
      subroutine get_born_nFKSprocess(nFKS_in,nFKS_out)
      implicit none
      include 'nexternal.inc'
      include 'nFKSconfigs.inc'
      include 'fks_info.inc'
      integer nFKS_in,nFKS_out,iFKS,iiFKS,nFKSprocessBorn(fks_configs)
      logical firsttime
      data firsttime /.true./
      save nFKSprocessBorn
c
      if (firsttime) then
         firsttime=.false.
         do iFKS=1,fks_configs
            nFKSprocessBorn(iFKS)=0
            if ( need_color_links_D(iFKS) .or. 
     &           need_charge_links_D(iFKS) )then
               nFKSprocessBorn(iFKS)=iFKS
            endif
            if (nFKSprocessBorn(iFKS).eq.0) then
c     try to find the process that has the same j_fks but with i_fks a
c     gluon
               do iiFKS=1,fks_configs
                  if ( (need_color_links_D(iiFKS) .or.
     &                  need_charge_links_D(iiFKS)) .and.
     &                 fks_j_D(iFKS).eq.fks_j_D(iiFKS) ) then
                     nFKSprocessBorn(iFKS)=iiFKS
                     exit
                  endif
               enddo
            endif
c     try to find the process that has the j_fks initial state if
c     current j_fks is initial state (and similar for final state j_fks)
            if (nFKSprocessBorn(iFKS).eq.0) then
               do iiFKS=1,fks_configs
                  if ( need_color_links_D(iiFKS) .or.
     &                 need_charge_links_D(iiFKS) ) then
                     if ( fks_j_D(iiFKS).le.nincoming .and.
     &                    fks_j_D(iFKS).le.nincoming ) then
                        nFKSprocessBorn(iFKS)=iiFKS
                        exit
                     elseif ( fks_j_D(iiFKS).gt.nincoming .and.
     &                        fks_j_D(iFKS).gt.nincoming ) then
                        nFKSprocessBorn(iFKS)=iiFKS
                        exit
                     endif
                  endif
               enddo
            endif
c     If still not found, just pick any one that has a soft singularity
            if (nFKSprocessBorn(iFKS).eq.0) then
               do iiFKS=1,fks_configs
                  if ( need_color_links_D(iiFKS) .or.
     &                 need_charge_links_D(iiFKS) ) then
                     nFKSprocessBorn(iFKS)=iiFKS
                  endif
               enddo
            endif
c     if there are no soft singularities at all, just do something trivial
            if (nFKSprocessBorn(iFKS).eq.0) then
               nFKSprocessBorn(iFKS)=iFKS
            endif
         enddo
         write (*,*) 'Total number of FKS directories is', fks_configs
         write (*,*) 'For the Born we use nFKSprocesses:'
         write (*,*)  nFKSprocessBorn
      endif
      if (nFKSprocessBorn(nFKS_in).eq.0) then
         write(*,*) 'Could not find the correct map to Born '/
     &        /'FKS configuration for the NLO FKS '/
     &        /'configuration', nFKS_in
         stop 1
      else
         nFKS_out=nFKSprocessBorn(nFKS_in)
      endif
      return
      end