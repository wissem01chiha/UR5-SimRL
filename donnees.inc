      real*8 long,lrg,epai,Iqua,Aire,YOUNG,NU,Gcis,
     & LIAIV(1000),LOADV(1000),CRIT,nlgm
     
      integer nmaix,nmaiy,nbnoe,LIAI(2,1000),NLIAI,LOAD(2,1000),
     & NLOAD,nbincr ,NORDRE,nctr,xctr,dincsto,pilot
     
      character*16 fdat,typ_elt(2),TYPLOAD(10)

      common /donnees/ long,lrg,epai,Iqua,Aire,YOUNG,NU,Gcis,
     & LIAIV,LOADV,CRIT,nlgm,nmaix,nmaiy,
     & nbnoe,LIAI,NLIAI,LOAD,NLOAD,nbincr,NORDRE,nctr,xctr,dincsto,
     & pilot,fdat,typ_elt,TYPLOAD
     
