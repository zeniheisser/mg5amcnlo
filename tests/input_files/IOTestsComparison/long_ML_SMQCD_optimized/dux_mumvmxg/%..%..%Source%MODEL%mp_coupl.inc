ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      written by the UFO converter
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      REAL*16 MP__G
      COMMON/MP_STRONG/ MP__G

      COMPLEX*32 MP__GAL(2)
      COMMON/MP_WEAK/ MP__GAL

      COMPLEX*32 MP__MU_R
      COMMON/MP_RSCALE/ MP__MU_R


      REAL*16 MP__MDL_MB,MP__MDL_MH,MP__MDL_MT,MP__MDL_MTA,MP__MDL_MW
     $ ,MP__MDL_MZ

      COMMON/MP_MASSES/ MP__MDL_MB,MP__MDL_MH,MP__MDL_MT,MP__MDL_MTA
     $ ,MP__MDL_MW,MP__MDL_MZ


      REAL*16 MP__MDL_WH,MP__MDL_WT,MP__MDL_WW,MP__MDL_WZ

      COMMON/MP_WIDTHS/ MP__MDL_WH,MP__MDL_WT,MP__MDL_WW,MP__MDL_WZ


      COMPLEX*32 MP__GC_11

      COMPLEX*32 MP__GC_4,MP__GC_5,MP__R2_GQQ,MP__R2_QQQ
     $ ,MP__UV_GQQG_1EPS,MP__UV_GQQB_1EPS,MP__UVWFCT_G_1_1EPS
     $ ,MP__R2_BXTW,MP__UV_GQQB,MP__UV_GQQT,MP__UVWFCT_G_1
     $ ,MP__UVWFCT_G_2

      COMMON/MP_COUPLINGS/ MP__GC_4,MP__GC_5,MP__R2_GQQ,MP__R2_QQQ
     $ ,MP__UV_GQQG_1EPS,MP__UV_GQQB_1EPS,MP__UVWFCT_G_1_1EPS
     $ ,MP__R2_BXTW,MP__UV_GQQB,MP__UV_GQQT,MP__UVWFCT_G_1
     $ ,MP__UVWFCT_G_2,MP__GC_11


