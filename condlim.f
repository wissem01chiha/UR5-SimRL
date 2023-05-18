      SUBROUTINE CALIMPNNL(MTG,LIAI,LIAIV,F,NLIAI,NDPN,NIV,NDL)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %        PRISE EN COMPTE DES DEPLACEMENTS IMPOSES NON NULS       %
COMM %        MODIFICATION DU SECOND MEMBRE F                         %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION F(*),LIAIV(*),MTG(NDL,*)
      DIMENSION LIAI(2,*)
    

      
      DO I=1,NLIAI
       N =LIAI(1,I) !noeud
       ND=LIAI(2,I) !ddl du noeud
       IF(ND.LE.NDPN) THEN
        NN=NDPN*(N-1)+ND+NIV
        DO J=1,NDL
         F(J)=F(J)-MTG(J,NN)*LIAIV(I) 
        ENDDO
       ELSE
        STOP ' Pb BLOCAGE : CALIMP  - Verifier les blocages - '
       ENDIF
      ENDDO


      
      DO I=1,NLIAI
       N =LIAI(1,I) !noeud
       ND=LIAI(2,I) !ddl du noeud
       NN=NDPN*(N-1)+ND+NIV
       F(NN)=LIAIV(I) ! insertion des déplacement imposés dans le SM
      ENDDO
      END     
      
      SUBROUTINE CALIMP0(LIAI,F,NLIAI,NDPN,NIV,NDL)
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %        PRISE EN COMPTE DES DEPLACEMENTS IMPOSES NON NULS       %
COMM %        MODIFICATION DU SECOND MEMBRE F                         %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION F(*)
      DIMENSION LIAI(2,*)
      
      DO I=1,NLIAI
       N =LIAI(1,I) !noeud
       ND=LIAI(2,I) !ddl du noeud
       NN=NDPN*(N-1)+ND+NIV
       F(NN)=0. ! insertion des déplacement imposés nuls dans le SM
      ENDDO
      
      END           
            
      SUBROUTINE CALIMP(LIAI,LIAIV,F,REAC,NLIAI,NDPN,NDL)
                       
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
COMM %                                                                %
COMM %        PRISE EN COMPTE DES DEPLACEMENTS IMPOSES NULS           %
COMM %        MODIFICATION DU SECOND MEMBRE F                         %
COMM %                                                                %
COMM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      DOUBLE PRECISION F(1000),LIAIV(1000),REAC(1000)
      DIMENSION LIAI(2,1000)
      
c       DO I=1,NDL
c       REAC(I)=REAC(I)+F(I)
c       ENDDO
      
      
      DO 1 I=1,NLIAI
      N =LIAI(1,I)
      ND=LIAI(2,I)
      IF(ND.LE.NDPN) THEN
        NN=NDPN*(N-1)+ND
        REAC(NN)=REAC(NN)+2*F(NN) !2*: prise en compte de l'integ dans largeur de la force
        F(NN)= 0 !LIAIV(I)

	
      ELSE
        STOP ' Pb BLOCAGE : CALIMP  - Verifier les blocages - '
      ENDIF
 1    CONTINUE
      END        
      
      
      
      SUBROUTINE CLTERMUNIT(MTG,LIAI,NLIAI,NDPN,NIV,NDL)
      DOUBLE PRECISION MTG(ndl,*)
      DIMENSION LIAI(2,*)      
      
      
      DO I=1,NLIAI
       N =LIAI(1,I) !noeud
       ND=LIAI(2,I) !ddl du noeud  
       IF(ND.LE.NDPN) THEN
        NN=NDPN*(N-1)+ND+NIV
        
        DO J=1,NDL
         IF(J.NE.NN)THEN
          MTG(NN,J)=0.
          MTG(J,NN)=0.         
         ELSE       
          MTG(NN,J)=1.
         ENDIF
         
        ENDDO
       ELSE
        STOP ' Pb BLOCAGE : CALIMP  - Verifier les blocages - '
       ENDIF
      ENDDO      
      
      END
      
      
      SUBROUTINE CLMUV_U(MTG,LIAI,NLIAI,NDPN,NIV,NDL)
      DOUBLE PRECISION LIAIV(1000),MTG(1000,1000)
      DIMENSION LIAI(2,1000)      
      
      
      DO I=1,NLIAI
       N =LIAI(1,I) !noeud
       ND=LIAI(2,I) !ddl du noeud  
       IF(ND.LE.NDPN) THEN
        NN=NDPN*(N-1)+ND+NIV
        
        DO J=1,NDL
         MTG(NN,J)=0.        
        ENDDO
        
       ELSE
        STOP ' Pb BLOCAGE : CALIMP  - Verifier les blocages - '
       ENDIF
      ENDDO      
      
      END     
      
      
      SUBROUTINE CLMUV_V(MTG,LIAI,NLIAI,NDPN,NIV,NDL)
      DOUBLE PRECISION LIAIV(1000),MTG(1000,1000)
      DIMENSION LIAI(2,1000)      
      
      
      DO I=1,NLIAI
       N =LIAI(1,I) !noeud
       ND=LIAI(2,I) !ddl du noeud  
       IF(ND.LE.NDPN) THEN
        NN=NDPN*(N-1)+ND+NIV
        
        DO J=1,NDL
         MTG(J,NN)=0.        
        ENDDO
        
       ELSE
        STOP ' Pb BLOCAGE : CALIMP  - Verifier les blocages - '
       ENDIF
      ENDDO      
      
      END         
