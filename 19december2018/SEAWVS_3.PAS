                        PROGRAM SEAWVS_3;
{        ëóÖí ÑÇìïäêÄíçéÉé àçíÖÉêÄãÄ
                  WIND WAVES + SLICK  !! §Æ°†¢®´® -kzxy
                  Æ‚ 3/10/2009 - Ø•‡•™‡•·‚≠Î© Á´•≠ · ¨®≠„·Æ¨
                  ØÆ xy ® · Ø´Ó·Æ¨ ØÆ xt
}

  USES
       GREAD_V,CRT,USEFUL_V;

  VAR
     INT,INTEG,A,B,C,D,E,E_2 : EXTENDED;
     NG1,N1X,N2X,N1Y,N2Y,NUM,M : LONGINT;

     { Ñãü îìçäñàâ }
     KM,WM,ALPHA,GAMMA,LIMIT_1K,LIMIT_2K,LIMIT_3K : EXTENDED;
     EH,LIMIT_H,LIMIT_1,LIMIT_2,LIMIT_3,LIMIT_4K, LL,U10T : EXTENDED;
     LIMIT_4,G : EXTENDED;
     INTEG_SMALL,SIGMAX,SIGMAY,U10,U10K,HW : EXTENDED;
     SIGMAH,SIGMAXX,SIGMAYY,SIGMATTX,FIW0,F0,FK,KL_3,U19 : EXTENDED;
     SIGMATT,KZXT,KZYT,KZXY, BETA0,SIGMA_W,SIGMA_M,EMASLO,UFTR: EXTENDED;
     UFRICTION,XRAZGON, XRG,OMEGA_M : EXTENDED;
     X_TILD,SLICK_KONTRAST : EXTENDED;

     N_POINT,I : INTEGER;
     F : TEXT;
     LOG_LARGE,LOG_SLICK : BOOLEAN;

{
 ---------------------------------------------------------------------------
}
         FUNCTION   FFF(X,Z : EXTENDED) : EXTENDED;
{
         Ñãü   èéàëäÄ   çìãÖâ, X - ëäéêéëíú íêÖçàü, Z - ÇõëéíÄ ,
	     èéëíéüççÄü äÄêåÄçÄ 0.4
}
     VAR
          Z0 : EXTENDED;
     BEGIN
         Z0:= 0.684/X + 0.0000428*X*X - 0.0443;
         FFF:= X/0.4 * LN(Z/Z0)
     END;

{
---------------------------------------------------------------------------
}
	     FUNCTION  WIND(UBEG,ZBEG,ZEND : EXTENDED;VAR UFRICTION:EXTENDED)
															     : EXTENDED;
{
           îìçäñàü çÄïéÑàí ëäéêéëíú ÇÖíêÄ çÄ çìÜçéâ ÇõëéíÖ ZEND
           èé ÇÇÖÑÖççéâ ÇõëéíÖ à ëäéêéëíà ÇÖíêÄ. äêéåÖ íéÉé ÇõÑÄÖí
           àçîéêåÄñàû é ëäéêéëíà íêÖçàü
}

     VAR
		 YWINDOW,XWINDOW : BYTE;

{
***** çÄóÄãé ÇçìíêÖççÖâ îìçäñàà *****
*************************************
}
                FUNCTION  FIND_FRICTION(X0,HX,EPSX: EXTENDED) : EXTENDED;
                {
                        çÄïéÑàí  çìãà  îìçäñàà. Ç  ÑÄççéå  ëãìóÄÖ
                        çÄïéÑàí ëäéêéëíú íêÖçàü
                }
		    VAR XA,XB,X : EXTENDED;
			    YA,YB,Y : BOOLEAN;
		    BEGIN
				XA:= X0;
				XB:= X0 + HX;

				YA:= (UBEG - FFF(XA,ZBEG)) > 0.0;
				YB:= (UBEG - FFF(XB,ZBEG)) > 0.0;

				WHILE  ( (YA AND YB) OR (NOT YA AND NOT YB) )
				DO   BEGIN
						XA:= XB;
						YA:= YB;
						XB:= XA + HX;
						YB:= (UBEG - FFF(XB,ZBEG)) > 0.0;
					 END;

				WHILE  (XB-XA) > EPSX
				DO   BEGIN
						 X:= (XA + XB) / 2.0;
						 Y:= (UBEG - FFF(X,ZBEG)) > 0.0;

						 IF  (YA AND Y) OR (NOT YA AND NOT Y)
						 THEN   BEGIN
									  XA:= X;
									  YA:= Y
								END
						 ELSE   BEGIN
									  XB:= X;
									  YB:= Y;
								END;
					 END;

				FIND_FRICTION:= (XA+XB) / 2.0;

		    END;

{**** äéçÖñ ÇëèéåéÉÄíÖãúçéâ îìçäñàà ****
****************************************
}

	BEGIN     {**** ÉãÄÇçÄü îìçäñàü ****}
        YWINDOW:= WHEREY;
        WRITELN;
        WRITE_NOVISIBLE('         ëäéêéëíú ÇÖíêÄ à ÇõëéíÄ ÇÇéÑüíëü Ç ëå/C ',
                                                                       1200);
        WRITELN;

        IF  UBEG > 8000.0
        THEN BEGIN
                 WRITELN(' ëãàòäéå ëàãúçõâ ÇÖíÖê ');
                 HALT;
             END;
	    IF  UBEG < 40.0
        THEN BEGIN
                 WRITELN(' ÇÖíÖê éóÖçú ëãÄÅõâ ');
                 HALT;
             END;

        UFRICTION:= FIND_FRICTION(0.0001,80.0,0.001);
        WIND:= FFF(UFRICTION,ZEND);
    END;

{
 ------------------------------------------------------------------------
}
         FUNCTION  ANGLE_OUR(K,KM: EXTENDED) : EXTENDED;
{ èêéÉêÄååÄ ÑÄÖí ìÉãéÇéÖ êÄëèêÖÑÖãÖçàÖ Ñãü çÄòÖâ îìçäñàà
  ÉÑÖ ÇÇÖÑÖçÄ ÑÇéâäÄ (ÇõóàëãÖç àçíÖÉêÄã ÄçÄãàíàóÖëäà)
}
     VAR
         Y,X,A,B : EXTENDED;

     BEGIN
          X:= K/KM;

          IF (X < 0.05) OR (X > 400000.0)
          THEN BEGIN
                 WRITELN(' FUNCTION ANGLE_OUR = [0.05, 400000] ');
                 HALT;
               END;

          IF (K/KM) < 0.4
          THEN K:= KM*0.4;
          B:= -0.28 + 0.65*EXPO(-0.75*LN(K/KM))
                    + 0.01*EXPO(-0.2 + 0.7*LOG10(K/KM));
          B:= POWER(10.0,B);

          A:= 2.0*B;

          Y:= 2.0 * ARCTAN(SH(PI*A)) / A;
          ANGLE_OUR:= 1.0/Y
     END;

{
 -----------------------------------------------------------------------
}
         FUNCTION  DIRECTION5(FIW,KM,K: EXTENDED) : EXTENDED;
{     åéÑàîàñàêéÇÄççéÖ êÄëèêÖÑÖãÖçàÖ  ìÉãéÇéÖ BANNER, 1990, J. PH. OC.
      ÑéÅÄÇãÖç äéùîîàñàÖçí 2 , àáåÖçÖçé áçÄóÖçàÖ -0.4
      à èÖêÖïéÑ éí éÑçéâ åéÑÖãà ä ÑêìÉéâ
}
     VAR
         Y,B : EXTENDED;

     BEGIN
          IF FIW > PI
          THEN FIW:= 2.0*PI - FIW
          ELSE IF FIW < -PI
               THEN FIW:= 2.0*PI + FIW;

          IF (K/KM) < 0.4
          THEN K:= KM*0.4;

          B:= -0.28 + 0.65*EXPO(-0.75*LN(K/KM))
                    + 0.01*EXPO(-0.2 + 0.7*LOG10(K/KM));
          B:= POWER(10.0,B);
          Y:= 2.0 /( EXPO(2.0*B*FIW) + EXPO(-2.0*B*FIW) );

         DIRECTION5:= Y;
	 END;

{
-----------------------------------------------------------------------------
                         é è à ë Ä ç à Ö   î ì ç ä ñ à â
}
         FUNCTION  FIND_DECISION(WX: EXTENDED) : EXTENDED;
{         çÄïéÑàí  êÖòÖçàÖ ìêÄÇçÖçàü íêÖíÖâ ëíÖèÖçà
          ÅÖá äéåèãÖäëçõï äéêçÖâ
          WX - óÄëíéíÄ
          P - äéùîîàñàÖçí èêà èÖêÇéâ ëíÖèÖçà
          Q - ëÇéÅéÑçõâ óãÖç
          ÇõóàëãüÖíëü ä àá ÑàëèÖêëàéçäà
}
     VAR
         KX,X1,X2,P,Q : EXTENDED;

	 BEGIN
        P:= 9.8*1000.0/0.074;
        Q:= -1000.0*WX*WX/0.074;

        X1:= -Q/2.0 + SQRT( SQR(Q/2.0) + POWER(P/3.0,3.0) );
        X2:= -Q/2.0 - SQRT( SQR(Q/2.0) + POWER(P/3.0,3.0) );

        KX:= POWER(X1,1.0/3.0) - POWER(-X2,1.0/3.0);
        FIND_DECISION:= KX;
     END;

{
 --------------------------------------------------------------------------
}
         FUNCTION  OMEGA_K(X : EXTENDED) : EXTENDED;
{             îìçäñàü ÇõèéãçüÖí èÖêÖëóÖí ÇéãçéÇéÉé óàëãÄ
         Ç OMEGA èé èéãçéâ ÑàëèÖêëàéçäÖ
}
    VAR
        W : EXTENDED;

    BEGIN
        W:= SQRT(9.8*X + X*X*X*0.074/1000.0);
        OMEGA_K:= W
    END;

{
 --------------------------------------------------------------------------
}
         FUNCTION  KM_LIMIT(U10 : EXTENDED) : EXTENDED;
{        áÄÇàëàåéëíú  èêÖÑÖãÄ éí ëäéêéëíà ÇÖíêÄ
}

    BEGIN
         IF (U10 < 2.5) OR (U10 > 22.0)
         THEN BEGIN
                 WRITELN(' WIND SPEED IS NOT CORRECT ');
                 WRITELN(' FUNCTION KM_LIMIT, éí [2.5 - 22.0 M/S] ');
                 HALT
              END;

         KM_LIMIT:=    0.371347584096022408
                     + 0.290241610467870486*U10
                     + 0.290178032985796564/U10;
    END;


{
 -------------------------------------------------------------------------
}
         FUNCTION  KL_LIMIT(U10,LL : EXTENDED) : EXTENDED;
{        áÄÇàëàåéëíú  §•´•≠®Ô èé ä ÇéãçÖçàü
}
    VAR
         Y : EXTENDED;
         NLS : INTEGER;

    BEGIN
         IF (U10 < 2.0) OR (U10 > 23.0)
         THEN BEGIN
                 WRITELN(' WIND SPEED IS NOT CORRECT ');
                 WRITELN(' FUNCTION KM_LIMIT [2, 23 M/S]');
                 HALT
              END;

         NLS:= ROUND(1000*LL);

         CASE NLS OF

             8 :  BEGIN
                         Y:=   5.95638868922385001
                             + 0.0720347041917690456*U10
                             + 584.172993830735945/SQRT(U10)
                             - 495.241636423855795*EXPO(-U10);
                   END;
             21 :  BEGIN
                         Y:=   171.844249137808972
                             - 29.9715115049259685*U10
                             + 11.0462317475020238*U10*LN(U10)
                             - 0.872695543932866229*U10*U10
                             + 0.139324811425696058*U10*U10*LN(U10)
                   END;
             30 :  BEGIN
                         Y:=   111.815140123056191
                             - 17.9371868277202707*U10
                             + 3.09504183073677383*U10*U10
                             - 1.2655962097427494*U10*U10*LN(U10)
                             + 0.314107775309787268*POWER(U10,2.5)
                   END;
             55:   BEGIN
                          Y:=   7.49318039534745137
                              - 0.1074769345744707*U10
                              + 0.494588232199318583*SQRT(U10)
                              + 81.5368708904397176*LN(U10)/U10
                              + 43.8821633213996845/SQR(U10)
                   END;
             100:  BEGIN
                          Y:=   3.50331310640510917
                              - 0.13811269308311842*U10
                              + 0.00839835244584838482*SQR(U10)
                              - 0.000303039793256592249*U10*SQR(U10)
                              + 4.47267177226626458e-06*SQR(U10)*SQR(U10);
                          Y:= EXPO(Y);
                   END;
             230:  BEGIN
                          Y:=   2.90160061283935079
                              - 0.187283073237593592*SQRT(U10)
                              + 17.8301004517745961*LN(U10)/U10
                              + 10.337585052957506/SQR(U10)
                   END;
         ELSE  BEGIN
                  WRITELN(' SUCH WAVELENGHT IS ABSENT IN PROGRAM ');
                  WRITELN(' FUNCTION KL_LIMIT ');
                  HALT;
               END;
         END;

         KL_LIMIT:= Y;

    END;

{
 --------------------------------------------------------------------------
}
           FUNCTION  DET(X : EXTENDED) : EXTENDED;
{             îìçäñàü êÄÇçÄ ÑÖíÖêåàçÄçíì èêà
              èÖêÖïéÑÖ éí óÄëíéíõ ä ÇéãçéÇõå óàëãÄå
              èé èéãçéâ ÑàëèÖêëàéçäÖ
              X - ÇéãçéÇéÖ óàëãé
}
    VAR
        D : EXTENDED;

    BEGIN
        D:=    ( 9.8 + 3.0*SQR(X)*0.074/1000.0 )
             / ( 2.0*SQRT(9.8*X + X*X*X*0.074/1000.0) );
        DET:= D;
    END;

{
 -------------------------------------------------------------------------
                  PARAMETERS FOR DEVELOPING WIND WAVE
}
              FUNCTION OMEGA_DEV(X : EXTENDED) : EXTENDED;
{         ÇõóàëãüÖíëü OMEGA~  áÄÇàëüôÖÖ éí êÄáÉéçÄ X
          ùíé ÅÖáêÄáåÖêçõâ êÄáÉéç
}
     VAR
         Y: EXTENDED;
     BEGIN
          IF X < 1430.0
          THEN BEGIN
                  WRITELN(' WIND FETCH IS LESS 1430 , FUNCTION OMEGA_DEV ');
                  HALT;
               END;

          IF X >= 20170.0
          THEN OMEGA_DEV:= 0.835
          ELSE BEGIN
                 Y:= 0.61826357843576103 + 3.52883010586243843e-06*X
                          - 0.00197508032233982112*SQRT(X)
                          + 62.5540113059129759/SQRT(X)
                           - 290.214120684236224/X;

                 OMEGA_DEV:= Y
               END;
     END;

              FUNCTION GAMMA_DEV(X : EXTENDED) : EXTENDED;
{         ÇõóàëãüÖíëü GAMMA~  áÄÇàëüôÖÖ éí êÄáÉéçÄ X
          ùíé ÅÖáêÄáåÖêçõâ êÄáÉéç
}
     VAR
         Y: EXTENDED;
     BEGIN
          IF X < 1430.0
          THEN BEGIN
                  WRITELN(' WIND FETCH IS LESS 1430 , FUNCTION GAMMA_DEV ');
                  HALT;
               END;

          IF X >= 20170.0
          THEN GAMMA_DEV:= 1.0
          ELSE BEGIN
                 Y:=  5.25366092888284494 + 0.000107621806575821253*X
                            - 0.0377677617998704552*SQRT(X)
                            - 162.983465281012658/SQRT(X)
                            + 253251.456471558638/(X*SQRT(X));
                 GAMMA_DEV:= Y;
               END;
     END;

              FUNCTION ALPHA_DEV(X : EXTENDED) : EXTENDED;
{         ÇõóàëãüÖíëü  ALPHA~ áÄÇàëüôÄü éí êÄáÉéçÄ X
          ùíé ÅÖáêÄáåÖêçõâ êÄáÉéç
}
     VAR
         Y: EXTENDED;
     BEGIN
          IF X < 1430.0
          THEN BEGIN
                  WRITELN(' WIND FETCH IS LESS 1430 , FUNCTION ALPHA_DEV ');
                  HALT;
               END;

          IF X >= 20170.0
          THEN ALPHA_DEV:= 0.0081
          ELSE BEGIN
                 Y:= 0.0311936666714071662 - 0.0023277357121814146*LN(X)
                              - 8367.86787861323644/SQR(X)
                              + 4.51145991018772117e+617*EXPO(-X);
                 ALPHA_DEV:= Y;
               END;
     END;

{
----------------------------------------------------------------------------
				   ******* ëèÖäíê ÇéãçÖçàü ********
}
		   FUNCTION   SPECTR_WAVE(X : EXTENDED) : EXTENDED;
{                ëèÖäíê JONSWAP
                X - ÇéãçéÇéÖ óàëãé
                èêà àçíÖÉêàêéÇÄçàà çÖ çìÜçé ÑéåçéÜÄíú çÄ ÇéãçéÇéÖ
                óàëãé, í.Ö. KdK
}
	VAR
		G,SW,SW1,ARG,B: EXTENDED;

	BEGIN
			G:= 9.8;

			IF  X < KM
			THEN  B:= 0.07
			ELSE  B:= 0.09;

			SW:= (ALPHA/2.0) * POWER(X,-3.0);
		    ARG:= -1.25 * POWER(KM/X,2.0);
			SW:= SW * EXPO(ARG);
			SW1:= EXPO( - SQR( SQRT(X/KM) - 1.0 ) / (2.0*SQR(B)) );
			SW1:= POWER(GAMMA,SW1);
			SW:= SW*SW1;

            SPECTR_WAVE:= SW;
	END;


{
----------------------------------------------------------------------------
}
          FUNCTION   NEW_SPECTR1(KT : EXTENDED) : EXTENDED;
{                 ëèÖäíê ÇéãçÖçàü èêÖÑãéÜÖççõâ äàíÄâÉéêéÑëäàå
        çÄóàçÄÖíëü  éí 1.5 Ñé 3.5 WM
        W - óÄëíéíÄ
        ! çÖí äéùîîàñàÖçíÄ BETA, éç àôÖíëü ëíõäéÇäéâ
        ëèÖäíêéÇ Ç íéóäÖ LIMIT_1K
}
   VAR
        Y,BETA0,W0 : EXTENDED;

   BEGIN
        W0:=  OMEGA_K(LIMIT_1K);
        BETA0:=  SPECTR_WAVE(LIMIT_1K)*POWER(W0,4.0) / DET(LIMIT_1K);

        W0:=  OMEGA_K(KT);
        Y:=   BETA0 / POWER(W0,4.0);

        NEW_SPECTR1:= Y * DET(KT);
   END;


{
----------------------------------------------------------------------------
}
          FUNCTION   NEW_SPECTR2(KT : EXTENDED) : EXTENDED;
{                 ëèÖäíê ÇéãçÖçàü çÄòÄ åéÑÖãú
        çÄóàçÄÖíëü  éí 3.5 Ñé 10.5 WM
        W - óÄëíéíÄ
        ! çÖí äéùîîàñàÖçíÄ BETA, éç àôÖíëü ëíõäéÇäéâ
        ëèÖäíêéÇ Ç íéóäÖ LIMIT_2K
}
   VAR
        Y,BETA0,W0 : EXTENDED;

   BEGIN
        W0:=  OMEGA_K(LIMIT_2K);
        BETA0:=  NEW_SPECTR1(LIMIT_2K)*POWER(W0,5.0) / DET(LIMIT_2K);

        W0:=  OMEGA_K(KT);
        Y:=   BETA0 / POWER(W0,5.0);

        NEW_SPECTR2:= Y * DET(KT);
   END;


{
----------------------------------------------------------------------------
}
          FUNCTION   NEW_SPECTR3(KT : EXTENDED) : EXTENDED;
{                 ëèÖäíê ÇéãçÖçàü
        çÄóàçÄÖíëü  éí 30  WM  Ñé 100 çÄòÄ åéÑÖãú
        W - óÄëíéíÄ
        ! çÖí äéùîîàñàÖçíÄ BETA, éç àôÖíëü ëíõäéÇäéâ
        ëèÖäíêéÇ Ç íéóäÖ LIMIT_3K
}
   VAR
        Y,BETA0,W0 : EXTENDED;

   BEGIN
        W0:=  OMEGA_K(LIMIT_3K);
        BETA0:=  NEW_SPECTR2(LIMIT_3K)*POWER(W0,2.7) / DET(LIMIT_3K);

        W0:=  OMEGA_K(KT);
        Y:=   BETA0 / POWER(W0,2.7);

        NEW_SPECTR3:= Y * DET(KT);
   END;

{
----------------------------------------------------------------------------
}
          FUNCTION   NEW_SPECTR4(KT : EXTENDED) : EXTENDED;
{                 ëèÖäíê ÇéãçÖçàü
        çÄóàçÄÖíëü  éí 30  WM  Ñé 100 çÄòÄ åéÑÖãú
        W - óÄëíéíÄ
        ! çÖí äéùîîàñàÖçíÄ BETA, éç àôÖíëü ëíõäéÇäéâ
        ëèÖäíêéÇ Ç íéóäÖ LIMIT_3K
}
   VAR
        Y,BETA0,W0 : EXTENDED;

   BEGIN
        W0:=  OMEGA_K(LIMIT_4K);
        BETA0:=  NEW_SPECTR3(LIMIT_4K)*POWER(W0,5.0) / DET(LIMIT_4K);

        W0:=  OMEGA_K(KT);
        Y:=   BETA0 / POWER(W0,5.0);

        NEW_SPECTR4:= Y * DET(KT);
   END;

{
 -----------------------------------------------------------------------
 }
       FUNCTION GAMMA_KONTRAST(X,E,OMEGA : EXTENDED) : EXTENDED;
{  X - VOLNOVOE CHISLO, BETA - INKREMENT RAZGONA
   SIGMA - VYZKOST' - 72 - VODA, 40 MASLO  E - UPRUGOST' PLENKI = 15-20 }
  VAR
     G,K,RO, NU,X1,X2,X3,X4,X5, GAMMA : EXTENDED;

  BEGIN
     G:= 981;
     K:= X;
     RO:= 1;
     NU:= 0.01;

     X1:= 2*NU*K*K/OMEGA;
     X2:= E*K*K*K*SQRT(X1)/(RO*OMEGA*OMEGA);
     X3:= E*E*POWER(K,6)*SQRT(1/X1)/(2*RO*RO*POWER(OMEGA,4));
     X4:= 2*E*K*K*K/(RO*OMEGA*OMEGA)*SQRT(X1);
     X5:= 2.0*E*E*POWER(K,6)/(RO*RO*POWER(OMEGA,4));

     GAMMA:= 2*NU*K*K*(X1 - X2 + X3)/(X1-X4+X5);

     GAMMA_KONTRAST:= GAMMA;

  END;


{
 -----------------------------------------------------------------------
 }
       FUNCTION KONTRAST(X,BETA0,SIGMA_W,SIGMA_M,E,UFTR:EXTENDED): EXTENDED;
{  X - VOLNOVOE CHISLO, SIGMA - VYZKOST' - 72 - VODA, 40 MASLO
  E - UPRUGOST' PLENKI UFTR - SKOROST' TRENIYA }
  VAR
     G,K,RO, KON, OMEGA_W, OMEGA_M, NU,G1,G0,BETA : EXTENDED;
     LOG_KONT : BOOLEAN;

  BEGIN
     G:= 981;
     K:= X;
     RO:= 1;
     NU:= 0.01;

     OMEGA_W:= SQRT(G*K+SIGMA_W/RO*K*K*K);
     OMEGA_M:= SQRT(G*K+SIGMA_M/RO*K*K*K);

     G0:= GAMMA_KONTRAST(K,0,OMEGA_W);
     G1:= GAMMA_KONTRAST(K,E,OMEGA_M);

     BETA:= BETA0*SQR(K*UFTR)/OMEGA_W;
     LOG_KONT:= TRUE;

     IF (BETA > 2*G0) AND (BETA > 2*G1)
     THEN BEGIN
            KON:= (BETA - 2.0*G1) / (BETA - 2.0*G0);
          END
     ELSE IF (BETA < 2.0*G0) AND (BETA < 2.0*G1)
          THEN BEGIN
                 KON:= (2.0*G0 - BETA) / (2.0*G1 - BETA);
               END
          ELSE BEGIN
{                 WRITELN(' LARGE CONTRAST ');
}
                 KON:= 0.1;
                 LOG_KONT:= FALSE;
               END;

     IF LOG_KONT
     THEN SLICK_KONTRAST:= KON
     ELSE KON:= SLICK_KONTRAST * KON;

     KONTRAST:= KON;
  END;




{
----------------------------------------------------------------------------
}
          FUNCTION   FULL_SPECTR2(KT : EXTENDED) : EXTENDED;
{         èéãçõâ ëèÖäíê - JONSWAP + ÑêìÉàÖ ëèÖäíêõ
        ëíõäéÇäÄ ëèÖäíêéÇ - W/WM
        Ñé LIMIT_1K - 1/W**5
     LIMIT_2K    WM - 1/(W**4)
     LIMIT_3K    1/(W**5)
     LIMIT_4K    1/(W**2.7)
                 1/W**5
}
      VAR
         FULL_SPECTRX, COEFFICIENT: EXTENDED;

   BEGIN
      IF KT <= LIMIT_1K
      THEN FULL_SPECTRX:=  SPECTR_WAVE(KT)
      ELSE IF KT <= LIMIT_2K
           THEN FULL_SPECTRX:=  NEW_SPECTR1(KT)
           ELSE IF KT <= LIMIT_3K
                THEN  FULL_SPECTRX:= NEW_SPECTR2(KT)
                ELSE  IF KT <= LIMIT_4K
                      THEN FULL_SPECTRX:= NEW_SPECTR3(KT)
                      ELSE FULL_SPECTRX:= NEW_SPECTR4(KT);

      IF LOG_SLICK
      THEN BEGIN
             COEFFICIENT:= KONTRAST(KT/100,BETA0,SIGMA_W,SIGMA_M,
                                             EMASLO,UFTR);

             FULL_SPECTR2:= FULL_SPECTRX*COEFFICIENT;
           END
      ELSE FULL_SPECTR2:= FULL_SPECTRX

   END;


{
 ----------------------------------------------------------------------
}
         FUNCTION  F1(X,Y : EXTENDED) : EXTENDED;
{
                   ÑÇìåÖêçÄü îìçäñàü ëèÖäíêÄ ÇéãçÖçàü
                   X - ÇéãçéÇéÖ óàëãé
                   Y - ìÉéã
                   ÑàëèÖêëàü Çõëéí
}
  VAR
      Z : EXTENDED;

  BEGIN
      Z:=  FULL_SPECTR2(X);
      Z:=  Z * DIRECTION5(Y,KM,X);

      F1:= Z * ANGLE_OUR(X,KM);

  END;


         FUNCTION  F2(X,Y : EXTENDED) : EXTENDED;
{
                   ÑÇìåÖêçÄü îìçäñàü ëèÖäíêÄ ÇéãçÖçàü
                   X - ÇéãçéÇéÖ óàëãé
                   Y - ìÉéã
                   ÑàëèÖêëàü çÄäãéçéÇ SIGMAXX
                   FIW0 - íÖäìôàâ ìÉéã,í.Ö. èéÇéêéí ÇéãçÖçàü
                   éíçéëàíÖãúçé éëà X
}
  VAR
      Z : EXTENDED;

  BEGIN
      Z:=  FULL_SPECTR2(X) * X*X * SQR(COS(Y));
      Z:=  Z * DIRECTION5(FIW0-Y,KM,X);

      F2:= Z * ANGLE_OUR(X,KM);
  END;


         FUNCTION  F3(X,Y : EXTENDED) : EXTENDED;
{
                   ÑÇìåÖêçÄü îìçäñàü ëèÖäíêÄ ÇéãçÖçàü
                   X - ÇéãçéÇéÖ óàëãé
                   Y - ìÉéã
                   ÑàëèÖêëàü çÄäãéçéÇ, SIGMAYY
}
  VAR
      Z : EXTENDED;

  BEGIN
      Z:=  FULL_SPECTR2(X) * X*X * SQR(SIN(Y));
      Z:=  Z * DIRECTION5(FIW0-Y,KM,X);

      F3:= Z * ANGLE_OUR(X,KM);
  END;


         FUNCTION  F4(X,Y : EXTENDED) : EXTENDED;
{
                   ÑÇìåÖêçÄü îìçäñàü ëèÖäíêÄ ÇéãçÖçàü
                   X - ÇéãçéÇéÖ óàëãé
                   Y - ìÉéã
                   éêÅàíÄãúçÄü ëäéêéëíú
                   SIGMATT èêéÖäñàü çÄ éëú X
}
  VAR
      Z,WT : EXTENDED;

  BEGIN
      WT:= OMEGA_K(X);

      Z:=  FULL_SPECTR2(X) * WT*WT * SQR(COS(Y));
      Z:=  Z * DIRECTION5(FIW0-Y,KM,X);

      F4:= Z * ANGLE_OUR(X,KM);
  END;


         FUNCTION  F5(X,Y : EXTENDED) : EXTENDED;
{
                   ÑÇìåÖêçÄü îìçäñàü ëèÖäíêÄ ÇéãçÖçàü
                   X - ÇéãçéÇéÖ óàëãé
                   Y - ìÉéã
                   éêÅàíÄãúçÄü ëäéêéëíú
                   SIGMATT
}
  VAR
      Z,WT : EXTENDED;

  BEGIN
      WT:= OMEGA_K(X);

      Z:=  FULL_SPECTR2(X) * WT*WT;
      Z:=  Z * DIRECTION5(FIW0-Y,KM,X);

      F5:= Z * ANGLE_OUR(X,KM);
  END;


         FUNCTION  F6(X,Y : EXTENDED) : EXTENDED;
{
                   ÑÇìåÖêçÄü îìçäñàü ëèÖäíêÄ ÇéãçÖçàü
                   X - ÇéãçéÇéÖ óàëãé
                   Y - ìÉéã
                   éêÅàíÄãúçÄü ëäéêéëíú
                   KZXT ≠•‚ ¨®≠„·†
}
  VAR
      Z,WT : EXTENDED;

  BEGIN
      WT:= OMEGA_K(X);

      Z:=  FULL_SPECTR2(X) * WT*X*COS(Y);
      Z:=  Z * DIRECTION5(FIW0-Y,KM,X);

      F6:= Z * ANGLE_OUR(X,KM);
  END;


         FUNCTION  F7(X,Y : EXTENDED) : EXTENDED;
{
                   ÑÇìåÖêçÄü îìçäñàü ëèÖäíêÄ ÇéãçÖçàü
                   X - ÇéãçéÇéÖ óàëãé
                   Y - ìÉéã
                   éêÅàíÄãúçÄü ëäéêéëíú
                   KZYT - ≠•‚ ¨®≠„·†
}
  VAR
      Z,WT : EXTENDED;

  BEGIN
      WT:= OMEGA_K(X);

      Z:=  - FULL_SPECTR2(X) * WT*X*SIN(Y);
      Z:=  Z * DIRECTION5(FIW0-Y,KM,X);

      F7:= Z * ANGLE_OUR(X,KM);
  END;


         FUNCTION  F8(X,Y : EXTENDED) : EXTENDED;
{
                   ÑÇìåÖêçÄü îìçäñàü ëèÖäíêÄ ÇéãçÖçàü
                   X - ÇéãçéÇéÖ óàëãé
                   Y - ìÉéã

                   KZXY - •·‚Ï ¨®≠„· - ¢·‚†¢®´® ¢ ™Æ≠•Á≠„Ó ‰„≠™Ê®Ó
}
  VAR
      Z,WT : EXTENDED;

  BEGIN
      WT:= OMEGA_K(X);

      Z:=  FULL_SPECTR2(X) * X*X*SIN(Y)*COS(Y);
      Z:=  Z * DIRECTION5(FIW0-Y,KM,X);

      F8:= Z * ANGLE_OUR(X,KM);
  END;


{
---------------------------------------------------------------------------
}
	 FUNCTION FGAUSS(X,Y:EXTENDED;M:INTEGER):EXTENDED;
	   BEGIN

			 CASE M OF
					   1: FGAUSS:= F1(X,Y);        { SIGMAH }
					   2: FGAUSS:= F2(X,Y);        { SIGMAXX }
					   3: FGAUSS:= F3(X,Y);        { SIGMAYY }
                       4: FGAUSS:= F4(X,Y);        { SIGMATT XX }
                       5: FGAUSS:= F5(X,Y);        { SIGMATT }
                       6: FGAUSS:= F6(X,Y);        { KZXT }
                       7: FGAUSS:= F7(X,Y);        { KZYT }
                       8: FGAUSS:= - F8(X,Y);        { KZXY }
             ELSE
                    WRITELN('  NUMBER OF FUNCTION IS NOT CORRECT ');
                    DELAY(3000);
                    HALT;
			 END;
	   END;

{
---------------------------------------------------------------------------
}
	 FUNCTION    GAUSS_Y(XT, A,B :EXTENDED ; M,NG:INTEGER) : EXTENDED ;
		 { M = çéåÖP îìHKñàà , NG - HOMEP BAPàAHTA ÉAìCCA }
	 VAR
		BMA ,BPA,GAUSSB : EXTENDED;
		GAUSSA: ARRAY[1..100] OF EXTENDED;
		I : INTEGER;
 BEGIN
		BMA := (B-A) / 2. ;
		BPA := (B+A) / 2. ;
		GAUSSB:=0.0;
   CASE NG OF
   8: BEGIN          { èOPüÑOK GAUSSA 8 }
	   GAUSSA[1]:=  0.362683783378362
                            * FGAUSS (XT, BMA*(+0.183434642495650) +BPA ,M);
	   GAUSSA[2]:=	0.313706645877887
                            * FGAUSS (XT, BMA*(+0.525532409916329) +BPA ,M);
	   GAUSSA[3]:=	0.222381034453374 * FGAUSS (XT, BMA*(+0.796666477413627) +BPA ,M);
	   GAUSSA[4]:=  0.101228536290376 * FGAUSS (XT, BMA*(+0.960289856497536) +BPA ,M);

	   GAUSSA[5]:=  0.362683783378362 * FGAUSS (XT, BMA*(-0.183434642495650) +BPA ,M);
	   GAUSSA[6]:=  0.313706645877887 * FGAUSS (XT, BMA*(-0.525532409916329) +BPA ,M);
	   GAUSSA[7]:=	0.222381034453374 * FGAUSS (XT, BMA*(-0.796666477413627) +BPA ,M);
	   GAUSSA[8]:=  0.101228536290376 * FGAUSS (XT, BMA*(-0.960289856497536) +BPA ,M);

	   FOR I:= 1 TO 8
	   DO  GAUSSB:=GAUSSB+GAUSSA[I];
	   GAUSS_Y:= GAUSSB*BMA
	  END;
   END;  { äéçÖñ ELSE }
 END;


{
---------------------------------------------------------------------------
}
	     FUNCTION   INTEGR_Y ( XT, C,D: EXTENDED; N1_Y, N2_Y:LONGINT;
                                   E:EXTENDED ; M,NG1: INTEGER) : EXTENDED;
{
      ÇõóàëãüÖå àçíÖÉêÄã èé Y,  XT - íÖäìôÖÖ áçÄóÖçàÖ îìçäñàà èé X
      KK - ëóÖíóàä éòàÅéä,  E  - éíçéëàíÖãúçÄü èOÉPEòHOCTú,
      M  - çéåÖê îìçäñàà,  C,D - èêÖÑÖãõ àçíÖÉêàêéÇÄçàü èé Y
      HMIN -  MàHàMAãúHõâ òAÉ èé Y,      NG1 - èéêüÑäà ÉÄìëëéÇõï èéãàçéåéÇ,
      HMAX -  MAKëàMAãúHõâ òAÉ èé Y
      N1_Y, N2_Y - Ñãü éèêÖÑÖãÖçàü ÇÖãàóàçõ òÄÉÄ - ãàçÖâçéÖ áçÄóÖçàÖ
      EH - êÖÉìãàêìÖí ÑêéÅãÖçàÖ òÄÉÄ èé  X
}

		 VAR
			 I,KK : LONGINT;
			 HY,INTNG1,INTNG2, FUN_MIN, FUN_PLUS:EXTENDED;
			 INTEGRAL1,Y0,F_ERRORY,HMAX,HMIN:EXTENDED;
             SIGN: BOOLEAN;

		 BEGIN   { çÄóÄãé }

                HMAX:= (D - C) / (N1_Y + 1);
                HMIN:= (D - C) / (N2_Y + 1);

{ ***** èéÑÉéíéÇäÄ ä ëóÖíì *****
********************************
}
			    HY:= HMIN;
				INTEGRAL1:= 0.0;
				Y0 := C;
				KK:= 0;
                F_ERRORY:= 0.0;
                SIGN:= FALSE;

{ ***** ëóÖí *****
******************
}
                REPEAT    { ÉãÄÇçõâ ñàäã }
			       IF Y0+HY >= D
			       THEN BEGIN
					      SIGN:=TRUE;
					      HY:= D-Y0        { ìMEHúòEHàE ÑãàHõ ùãEMEHTA }
					    END;

  { BTOPOâ  KPàTEPàâ }

	               INTNG1:=   GAUSS_Y( XT , Y0,Y0+HY, M,NG1);
	               INTNG2:=   GAUSS_Y( XT , Y0,Y0+HY/2.0, M,NG1)
                            + GAUSS_Y( XT , Y0+HY/2.0,Y0+HY, M,NG1);

                   FUN_MIN:= ABS( INTNG2 - INTNG1 );
                   FUN_PLUS:= ABS( INTNG2 + INTNG2 );

	               IF  ( FUN_MIN <= E*FUN_PLUS )
	               THEN BEGIN    {  ùãÖåÖçí ëéëóàíÄç  }
			              INTEGRAL1:= INTEGRAL1 + INTNG2;
                          Y0:= Y0 + HY;
                      { Ñãü ìåÖçúòÖçàü ëäéêéëíà êÄáÉéçÄ }
                          IF FUN_MIN  < E_2*E*FUN_PLUS
                          THEN  HY:=HY*2.0;
		                END
	               ELSE BEGIN    { ìåÖçúòÖçàÖ Ñãàçõ ùãÖåÖçíÄ }
                          IF FUN_PLUS/HY <= E*INTEGRAL1/(Y0-C)
                          THEN BEGIN
                                INTEGRAL1:= INTEGRAL1 + INTNG2;
                                Y0:= Y0 + HY;
                               END
                          ELSE BEGIN
				                 IF HY > HMIN/EH
				                 THEN BEGIN
					                    HY:= HY / 3.0;
					                    SIGN:= FALSE
					                  END
				                 ELSE BEGIN
					                    KK:= KK+1;
					                    Y0:= Y0+HY;
					                    INTEGRAL1:=INTEGRAL1 + INTNG2;
					                    HY:= HMIN;
                                        F_ERRORY:= F_ERRORY + (INTNG2-INTNG1);
				                      END
                               END;
		                END;

                UNTIL  SIGN ;     { äéçÖñ ÉãÄÇçéÉé ñàäãÄ }

	            INTEGR_Y:= INTEGRAL1
		 END;      { äéçÖñ îìçäñàà àçíÖÉêÄã }

{
----------------------------------------------------------------------------
}
	 FUNCTION    GAUSS_X(A,B,  C,D: EXTENDED;
                                N1_Y,N2_Y : LONGINT;
                                E:EXTENDED;
                                            M,NG,NG1:INTEGER): EXTENDED ;
{            M = çéåÖP îìHKñàà , NG - HOMEP BAPàAHTA ÉAìCCA
             N1_Y, N2_Y - áçÄóÖçàü Ñãü òÄÉÄ èé X
}
	 VAR
		BMA ,BPA,GAUSSI ,GAUSSJ: EXTENDED;
		GAUSSA: ARRAY[1..100] OF EXTENDED;
		I : INTEGER;
		W : ARRAY[1..100] OF EXTENDED;
		X : ARRAY[1..100] OF EXTENDED;

 BEGIN
		BMA := (B-A) / 2. ;
		BPA := (B+A) / 2. ;

   CASE NG OF

   8: BEGIN          { èOPüÑOK GAUSSA 8 }

	   W[1]:=  0.362683783378362 ;     X[1]:= +0.183434642495650;
	   W[2]:=  0.313706645877887 ;     X[2]:= +0.525532409916329;
	   W[3]:=  0.222381034453374 ;     X[3]:= +0.796666477413627;
	   W[4]:=  0.101228536290376 ;     X[4]:= +0.960289856497536;

	   W[5]:=  0.362683783378362 ;     X[5]:= -0.183434642495650;
	   W[6]:=  0.313706645877887 ;     X[6]:= -0.525532409916329;
	   W[7]:=  0.222381034453374 ;     X[7]:= -0.796666477413627;
	   W[8]:=  0.101228536290376 ;     X[8]:= -0.960289856497536;

	  END;

   END;
		   GAUSSI:=0.0;

        	 FOR I:=1 TO NG
			 DO	GAUSSI:= GAUSSI + W[I]*INTEGR_Y(BMA*X[I]+BPA,
                                           C,D, N1_Y,N2_Y,E, M,NG1);

                GAUSS_X:= GAUSSI * BMA;
 END;  { END OF GAUSS_X }



{
---------------------------------------------------------------------------
}
	     FUNCTION   TWO_INTEGRALN ( A,B,C,D, E:EXTENDED ;
                                          M,NG1: INTEGER;
                                          N1_X,N2_X,N1_Y,N2_Y:LONGINT)
                                                                : EXTENDED;

{
N1_X - óàëãé åàçàåÄãúçéÖ êÄáÅàÖçàâ èé ï, N2_X - åÄäëàMAãúHOÖ óàCãé PAáÅàEHàâ,
N1_Y - óàëãé åàçàåÄãúçéÖ êÄáÅàÖçàâ èé Y, N2_Y - åÄäëàMAãúHOÖ óàCãé PAáÅàEHàâ,
  KK - ëóÖíóàä éòàÅéä , ÄçÄãéÉàóçé èé Y
  E  - éíçéëàíÖãúçÄü èOÉPEòHOCTú ,     A,B - èêÖÑÖãõ àçíÖÉêàêéÇÄçàü  èé X
  M - çéåÖê îìçäñàà,                   C,D - èêÖÑÖãõ àçíÖÉêàêéÇÄçàü èé Y
  HMIN_X -  MàHàMAãúHõâ òAÉ èé X,      NG1 - èéêüÑäà ÉÄìëëéÇõï èéãàçéåéÇ,
  HMAX_X -  MAKëàMAãúHõâ òAÉ èé X,  HMIN_Y - åàçàåÄãúçõâ òÄÉ èé Y
  HMAX_Y - åÄäëàåÄãúçõâ òÄÉ èé Y
}

		 VAR
		   I,KK : LONGINT;
		   HMIN_X,HMAX_X,HX,D_INTNG1,D_INTNG2:EXTENDED;
		   TWO_INTEGRAL1,X0,F_ERRORX:EXTENDED;
           SIGN : BOOLEAN;

	BEGIN   { çÄóÄãé }

		  IF ( A > B )
		  THEN   BEGIN
		  			WRITELN(' A > B Ç îìçäñàà INTEGRAL ');
					HALT;
				 END;
		   IF  NOT (NG1 = 8)
		   THEN   BEGIN
					  WRITELN(' éòàÅäÄ Ç ÇõÅéêÖ èéêüÑäÄ ÉÄìëÄ ');
					  WRITELN(' îìçäñàü ÑÇìïäêÄíçõâ àçíÖÉêÄã ');
					  HALT;
				  END;

{ *** èéÑÉéíéÇäÄ ä ëóÖíì *** }
		   HMAX_X:= (B - A) / (N1_X + 1);
		   HMIN_X:= (B - A) / (N2_X + 1);

		   HX:= HMIN_X;
		   TWO_INTEGRAL1:= 0.0;
           F_ERRORX:= 0.0;
		   X0 := A;
		   KK:= 0;
           SIGN:= FALSE;

          REPEAT    { ÉãÄÇçõâ ñàäã }
			 IF X0+HX >= B
			 THEN   BEGIN
					   SIGN:=TRUE;
					   HX:=B-X0        { ìMEHúòEHàE ÑãàHõ ùãEMEHTA }
					END;

	  { BTOPOâ  KPàTEPàâ }
		     D_INTNG1:= GAUSS_X( X0,X0+HX , C,D, N1_Y,N2_Y, E,M,NG1,NG1);
	         D_INTNG2:= GAUSS_X( X0,X0+HX/2.0 , C,D, N1_Y,N2_Y, E,
                                                            M,NG1,NG1)+
                        GAUSS_X( X0+HX/2.0,X0+HX , C,D, N1_Y,N2_Y, E,
                                                            M,NG1,NG1);

           IF ( ABS(D_INTNG2 - D_INTNG1) <= E*ABS(D_INTNG2 + D_INTNG1) )
           THEN BEGIN
		          TWO_INTEGRAL1:= TWO_INTEGRAL1 + D_INTNG2;
                  X0:= X0 + HX;
                               { Ñãü ìåÖçúòÖçàü ëäéêéëíà êÄáÉéçÄ }
                  IF ABS(D_INTNG2-D_INTNG1)
                              < E_2*E*(ABS(D_INTNG2)+ABS(D_INTNG1))
                  THEN  HX:=HX*1.8;
		        END
	       ELSE BEGIN    { ìÇÖãàóÖçàÖ Ñãàçõ ùãÖåÖçíÄ }
                  IF (X0-A)*ABS(D_INTNG2 + D_INTNG1)/HX <= E*TWO_INTEGRAL1
                  THEN BEGIN
		                 TWO_INTEGRAL1:= TWO_INTEGRAL1 + D_INTNG2;
                         X0:= X0 + HX;
                       END
                  ELSE BEGIN
			             IF HX > HMIN_X/EH
			             THEN BEGIN
				                HX:=HX / 3.0;
				                SIGN:= FALSE
				              END
			             ELSE BEGIN
				                KK:= KK+1;
                                F_ERRORX:= F_ERRORX + (D_INTNG2-D_INTNG1);
			                    X0:= X0+HX;
					            TWO_INTEGRAL1:=TWO_INTEGRAL1+D_INTNG2;
					            HX:= HMIN_X;
				              END
                       END;
		        END;

          UNTIL  SIGN ;     { äéçÖñ ÉãÄÇçéÉé ñàäãÄ }

	      TWO_INTEGRALN:= TWO_INTEGRAL1
	END;      { äéçÖñ îìçäñàà àçíÖÉêÄã }



{***************************************************************************}
{                             ÉãÄÇçÄü óÄëíú                                 }

 BEGIN
   TEXTBACKGROUND(BLACK);
   CLRSCR;
   TEXTCOLOR(10);

WRITELN;
WRITELN('                ⁄ƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒø');
WRITELN('                ≥       èêéÉêÄååÄ SEAWVS_3.PAS    ≥');
WRITELN('                ¿ƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒƒŸ');
WRITELN;

   WINDOW(1, 5, 80, 25);
   CLRSCR;
   TEXTCOLOR(12);

   NG1:= 8; { IREAD(' ÇÇÖÑàíÖ NG1 '); }
   N1X:=  30;{IREAD(' óàëãé êÄáÅàÖçàâ èé X (¨®≠®¨†´Ï≠Æ• Á®·´Æ,´®≠•©≠Æ) ');}
   N2X:= N1X*100;
   N1Y:= 1; {IREAD(' óàëãé êÄáÅàÖçàâ èé Y (¨®≠®¨†´Ï≠Æ• Á®·´Æ,´®≠•©≠Æ) ');}
   N2Y:= N1Y*100;

   E:= 0.0001;      { éíçéëàíÖãúçÄü èéÉêÖòçéëíú }
   E_2:= 0.2;       { èéÉêÖòçéëíú ìÇÖãàóÖçàü òÄÉÄ }
   EH := 100000.0;  { ìåÖçúòÖçàÖ åàçàåÄãúçéÉé òÄÉÄ }

   IF ANT(' REWRITTEN FILE ? ')
   THEN BEGIN
           ASSIGN(F,'SEAWVS_3.DAT');
           REWRITE(F);
        END
   ELSE HALT;


{***** INITIAL SET *****
 *******************************}
   G:= 9.8;
   X_TILD:= RREAD(' INPUT WIND FETCH (1430-20170) ');
   GAMMA:= GAMMA_DEV(X_TILD);
   ALPHA:= ALPHA_DEV(X_TILD);
   OMEGA_M:= OMEGA_DEV(X_TILD);
   WRITELN(' ALPHA= ',ALPHA:6:5,' GAMMA= ',GAMMA:5:4);

   LOG_SLICK:= ANT(' WORK WITH SLICK ? ');
   IF LOG_SLICK
   THEN BEGIN
          WRITELN(' PARAMETERS OF SLICK ');
          BETA0:= RREAD(' INPUT  BETA0 (0.04) ');
          SIGMA_W:= RREAD(' INPUT SIGMA OF WATER (72) ');
          SIGMA_M:= RREAD(' INPUT SIGMA OF OIL (40) ');
          EMASLO:= RREAD(' INPUT E OF OIL (5-30) ');
        END;

{ =======================================================================}
         M:= 1;                {ÇÇéÑàíëü çéåÖê îìçäñàà}

         LIMIT_1:= 1.2;   { çÄóàçÄÖíëü  - 1/(W**4) }

         LIMIT_3:= 270.0; { çÄóàçÄÖíëü 1/W**2.7 }
         LIMIT_4:= 1020.0; { çÄóàçÄÖíëü 1/W*5 }

         LIMIT_H:= 6000.0; { ÇÖêïçüü ÉêÄçàñÄ àçíÖÉêàêéÇÄçàü }
                                { §Æ·‚†‚ÆÁ≠†Ô ‚ÆÁ≠Æ·‚Ï, °Æ´ÏË• ≠• ‡†·‚•‚}
         WRITELN('          L = 0.008, 0.021, 0.03, 0.055, 0.1, 0.23 M ');
         LL:= RREAD('Input radar wavelength L ');

         IF ANT(' DEPENDENCE ON WIND SPEED ? ')
	     THEN BEGIN
                WRITELN(F,' U10T  ','  FI0      SIGMAH          SIGMAXX   ',
                     '  SIGMAYY   ',' SIGMATTX    ','    SIGMATT   ',
                     '   KZXT   ','   KZYT   ','   KZXY ');


		        WRITELN(' CYCLE ON WIND SPEED ');
			    U10:=RREAD(' Input initial wind speed U10, m/s ');
			    U10K:= RREAD(' Input final wind speed U10, m/s ');
			    N_POINT:= IREAD(' Input number of points ');
			    HW:= (U10K - U10) / (N_POINT-1.0);

                XRG:= X_TILD;

                FIW0:= RREAD(' Input DIRECTION OF WIND ');
                FIW0:= FIW0*PI/180.0;

                LOG_LARGE:= ANT(' ONLY Parameters of Large-scale Waves ? ');
                F0:= FIW0;

			    FOR I:= 1 TO N_POINT
			    DO BEGIN
			 	     U10T:= U10 + (I-1)*HW;

                     XRAZGON:= XRG*SQR(U10T)/9.8;
                     WRITELN(' wind fetch, m = ',XRAZGON:7:1);
                     WRITELN(' ALFA ',ALPHA:5:4,' GAMMA ',GAMMA:5:4);

                     U19:= WIND(U10T*100.0,1000.0,1950.0, UFRICTION)/100.0;
                     UFRICTION:= UFRICTION / 100.0;
                     UFTR:= UFRICTION*100;
                     WM:= OMEGA_M*G/U10T;
                     KM:= SQR(WM)/G;
                     WRITELN;

                     WRITELN('U10= ',U10T:5:2,
                             'M/S   U* friction= ',UFRICTION:4:2,' M/S');
                     WRITELN('KM = ',KM:7:5,'  WM= ',WM:7:5,
                             '     LM= ',2.0*PI/KM:5:2);

                     LIMIT_2:= KM_LIMIT(U10T);

                     A:= KM/4.0;
                     B:= KL_LIMIT(U10T,LL);

                     WRITELN(' A = ',A:10:8,' B = ',B:10:8);
                     C:= -PI + FIW0;
                     D:=  PI + FIW0;

                     LIMIT_1K:=  FIND_DECISION(LIMIT_1*WM);
                                               { çÄóàçÄÖíëü  - 1/(W**4) }
                     LIMIT_2K:=  FIND_DECISION(LIMIT_2*WM);
                                               { -//-  1/(W**5) }
                     LIMIT_3K:=  LIMIT_3; { 1/W**2.7 }
                     LIMIT_4K:=  LIMIT_4; { 1/W**5.0 }

                     WRITELN(' KM LIMIT_1 = ',LIMIT_1K:6:3,
                                 ' 2 = ',LIMIT_2K:6:3,
                                 ' 3 = ',LIMIT_3K:6:3,
                                 ' 4 = ',LIMIT_4K:6:3);
{***** Çõëéíõ *****
 ******************
}
                     INTEG:= TWO_INTEGRALN(A,B,C,D, E,
                                             M,NG1,N1X,N2X,N1Y,N2Y);
                     INTEG_SMALL:= TWO_INTEGRALN(B,LIMIT_H,C,D, E,
                                             M,NG1,N1X,N2X,N1Y,N2Y);
                     SIGMAH:= INTEG+INTEG_SMALL;
                     WRITELN(' VARIANCE OF HEIGHT:  (LARGE-SCALE) ',
                                 INTEG:9:7,
                                 ' SMALL-SCALE ',INTEG_SMALL:9:8);
                     WRITELN('              TOTAL VARIANCE OF HEIGHT = ',
                                                               SIGMAH:9:8);
                     WRITELN(' ------------------------------------- ');


{***** çÄäãéçõ *****
 *******************
}
            { XX }   INTEG:= TWO_INTEGRALN(A,B,C,D, E,
                                            M+1,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' VARIANCE OF SLOPES XX: LARGE-SCALE= ',
                                                               INTEG:9:8);

                     INTEG_SMALL:= TWO_INTEGRALN(B,LIMIT_H,C,D, E,
                                            M+1,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' VARIANCE OF SLOPES YY: SMALL-SCALE= ',
                                                          INTEG_SMALL:9:8);
                     IF LOG_LARGE
                     THEN  SIGMAXX:= INTEG
                     ELSE  SIGMAXX:= INTEG + INTEG_SMALL;

                     WRITELN('               TOTAL VARIANCE OF SLOPES XX= ',
                                                    (INTEG+INTEG_SMALL):9:8);

            { YY }   INTEG:= TWO_INTEGRALN(A,B,C,D, E,
                                                 M+2,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' VARIANCE OF SLOPES YY: LARGE-SCALE= ',
                                                                INTEG:9:8);

                     INTEG_SMALL:= TWO_INTEGRALN(B,LIMIT_H,C,D, E,
                                                   M+2,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' VARIANCE OF SLOPES YY: SMALL-SCALE= ',
                                                            INTEG_SMALL:9:8);
                     IF LOG_LARGE
                     THEN  SIGMAYY:= INTEG
                     ELSE  SIGMAYY:= INTEG + INTEG_SMALL;
                     WRITELN('              TOTAL VARIANCE OF SLOPE YY = ',
                                                    (INTEG+INTEG_SMALL):9:8);

{***** ëäéêéëíà *****
 ********************
}
            {TT XX}  INTEG:= TWO_INTEGRALN(A,B,C,D, E,
                                                 M+3,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN('VARIANCE OF ORBITAL VELOCITY XX: LARGE-SCALE= ',
                                                                  INTEG:9:8);

                     INTEG_SMALL:= TWO_INTEGRALN(B,LIMIT_H,C,D, E,
                                                 M+3,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN('VARIANCE OF ORB. VEL. XX: SMALL-SCALE= ',
                                                            INTEG_SMALL:9:8);
                     IF LOG_LARGE
                     THEN  SIGMATTX:= INTEG
                     ELSE  SIGMATTX:= INTEG + INTEG_SMALL;
                     WRITELN('              TOTAL VARIANCE OF ',
                             'ORBITAL VELOCITY XX= ',(INTEG+INTEG_SMALL):9:8);
{***** ëäéêéëíà *****
 ********************
}
              {TT}   INTEG:= TWO_INTEGRALN(A,B,C,D, E,
                                            M+4,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN('VARIANCE OF ORBITAL VELOCITY : LARGE-SCALE= ',
                                                                  INTEG:9:8);
                     INTEG_SMALL:= TWO_INTEGRALN(B,LIMIT_H,C,D, E,
                                            M+4,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' VARIANCE OF ORB. VEL.: SMALL-SCALE = ',
                                                            INTEG_SMALL:9:8);
                     IF LOG_LARGE
                     THEN  SIGMATT:= INTEG
                     ELSE  SIGMATT:= INTEG + INTEG_SMALL;
                     WRITELN('              TOTAL VARIANCE OF ORB. VEL. = ',
                                                   (INTEG+INTEG_SMALL):9:8);

              {XT}   INTEG:= TWO_INTEGRALN(A,B,C,D, E,
                                                M+5,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' KZXT CORRELATION : LARGE-SCALE = ',INTEG:9:8);

                     INTEG_SMALL:= TWO_INTEGRALN(B,LIMIT_H,C,D, E,
                                                 M+5,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' KZXT CORRELATION: SMALL-SCALE= ',
                                                          INTEG_SMALL:9:8);

                     IF LOG_LARGE
                     THEN  KZXT:= INTEG
                     ELSE  KZXT:= INTEG + INTEG_SMALL;
                     WRITELN('              TOTAL KZXT CORRELATION= ',
                                                   (INTEG+INTEG_SMALL):9:8);

             {YT}    INTEG:= TWO_INTEGRALN(A,B,C,D, E,
                                                 M+6,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' KZYT CORRELATION: LARGE-SCALE = ',INTEG:9:8);

                     INTEG_SMALL:= TWO_INTEGRALN(B,LIMIT_H,C,D, E,
                                            M+6,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' KZYT CORRELATION: SMALL-SCALE= ',
                                                           INTEG_SMALL:9:8);
                     IF LOG_LARGE
                     THEN  KZYT:= INTEG
                     ELSE  KZYT:= INTEG + INTEG_SMALL;
                     WRITELN('              TOTAL KZYT CORRELATION= ',
                                                    (INTEG+INTEG_SMALL):9:8);

             {XY}    INTEG:= TWO_INTEGRALN(A,B,C,D, E,
                                            M+7,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' KZXY CORRELATION: LARGE-SCALE= ',INTEG:9:8);

                     INTEG_SMALL:= TWO_INTEGRALN(B,LIMIT_H,C,D, E,
                                            M+7,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' KZYT CORRELATION: SMALL-SCALE= ',
                                                           INTEG_SMALL:9:8);

                     IF LOG_LARGE
                     THEN  KZXY:= INTEG
                     ELSE  KZXY:= INTEG + INTEG_SMALL;
                     WRITELN('              TOTAL KZXY CORRELATION= ',
                                                  (INTEG+INTEG_SMALL):9:8);

                     WRITELN(' ----------------------------------------- ');
                     WRITELN(' SIGMAXX =  ',SIGMAXX:9:8);
                     WRITELN(' SIGMAYY =  ',SIGMAYY:9:8);
                     WRITELN(' SIGMATTX=  ',SIGMATTX:9:8);
                     WRITELN(' SIGMATT =  ',SIGMATT:9:8);
                     WRITELN(' KZXT    =  ',KZXT:9:8);
                     WRITELN(' KZYT    =  ',KZYT:9:8);
                     WRITELN(' KZXY    =  ',KZXY:9:8);
                     WRITELN(' U10     =  ',U10T:4:2,'  FIW ',F0*180/PI:5:2);
                     WRITELN(' ------------------------------------------ ');

                     WRITELN(F,'  ',U10T:4:2,'   ',F0*180/PI:4:2,'        ',
                                       SIGMAH:9:8,'    ',
                                       SIGMAXX:7:6,'    ',SIGMAYY:7:6,'    ',
                                       SIGMATTX:7:6,'     ',
                                       SIGMATT:7:6,'    ',KZXT:7:6,'   ',
                                       KZYT:7:6,'    ',KZXY:7:6);

		           END;
				WRITELN(' CALCULATIONS WERE CARRIED OUT ')
		      END;

         IF ANT(' DEPENDENCE ON AZIMUTHAL ANGLE ? ')
	     THEN BEGIN
                WRITELN(F,' U10T  ','  FI0         SIGMAXX   ',
                     '  SIGMAYY   ',' SIGMATTX    ','    SIGMATT   ',
                     '   KZXT   ','   KZYT   ','   KZXY ');

	            WRITELN(' CYCLE ON AZIMUTHAL ANGLE ');
			    F0:= RREAD(' Input initial angle F0 ');
			    FK:= RREAD(' Input final angle FK  ');
                F0:= F0*PI/180.0;
                FK:= FK*PI/180.0;

			    N_POINT:= IREAD(' Input number of points ');
			    HW:= (FK - F0) / (N_POINT-1.0);
                LOG_LARGE:= ANT(' ONLY Parameters of Large-scale Waves ? ');

                U10T:= RREAD(' Input wind speed, M/S ');
                XRG:= X_TILD;

                XRAZGON:= XRG*SQR(U10T)/9.8;
                WRITELN(' WIND FETCH, M = ',XRAZGON:7:1);
                WRITELN(' ALFA ',ALPHA:5:4,' GAMMA ',GAMMA:5:4);

                U19:= WIND(U10T*100.0,1000.0,1950.0, UFRICTION)/100.0;
                UFRICTION:= UFRICTION / 100.0;
                UFTR:= UFRICTION*100;
                WM:= OMEGA_M*G/U10T;
                KM:= SQR(WM)/G;
                WRITELN;

                LIMIT_2:= KM_LIMIT(U10T);
                LIMIT_1K:=  FIND_DECISION(LIMIT_1*WM);
                                               { çÄóàçÄÖíëü  - 1/(W**4) }
                LIMIT_2K:=  FIND_DECISION(LIMIT_2*WM);
                                               { -//-  1/(W**5) }

                LIMIT_3K:=  LIMIT_3; { 1/W**2.7 }
                LIMIT_4K:=  LIMIT_4; { 1/W**5.0 }

                WRITELN(' KM LIMIT_1 = ',LIMIT_1K:6:3,
                                 ' 2 = ',LIMIT_2K:6:3,
                                 ' 3 = ',LIMIT_3K:6:3,
                                 ' 4 = ',LIMIT_4K:6:3);

                A:= KM/4.0;
                B:= KL_LIMIT(U10T,LL);
                WRITELN(' A = ',A,' B = ',B);

			    FOR I:= 1 TO N_POINT
			    DO BEGIN
                     C:= - PI + F0;
                     D:=   PI + F0;
                     FIW0:= F0;
                     WRITELN(' C = ',C*180/PI:5:2,' D = ',D*180/PI:5:2);

                     SIGMAH:= 0.0;

{***** çÄäãéçõ *****
 *******************
}
            { XX }   INTEG:= TWO_INTEGRALN(A,B,C,D, E,
                                            M+1,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' VARIANCE OF SLOPES XX: LARGE-SCALE= ',
                                                                  INTEG:9:8);

                     INTEG_SMALL:= TWO_INTEGRALN(B,LIMIT_H,C,D, E,
                                            M+1,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' VARIANCE OF SLOPES YY: SMALL-SCALE= ',
                                                            INTEG_SMALL:9:8);
                     IF LOG_LARGE
                     THEN  SIGMAXX:= INTEG
                     ELSE  SIGMAXX:= INTEG + INTEG_SMALL;

                     WRITELN('              TOTAL VARIANCE OF SLOPES XX= ',
                                                    (INTEG+INTEG_SMALL):9:8);

            { YY }   INTEG:= TWO_INTEGRALN(A,B,C,D, E,
                                                 M+2,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' VARIANCE OF SLOPES YY: LARGE-SCALE= ',
                                                                  INTEG:9:8);

                     INTEG_SMALL:= TWO_INTEGRALN(B,LIMIT_H,C,D, E,
                                                   M+2,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' VARIANCE OF SLOPES YY: SMALL-SCALE= ',
                                                             INTEG_SMALL:9:8);
                     IF LOG_LARGE
                     THEN  SIGMAYY:= INTEG
                     ELSE  SIGMAYY:= INTEG + INTEG_SMALL;
                     WRITELN('              TOTAL VARIANCE OF SLOPE YY = ',
                                                   (INTEG+INTEG_SMALL):9:8);
{***** ëäéêéëíà *****
 ********************
}
            {TT XX}  INTEG:= TWO_INTEGRALN(A,B,C,D, E,
                                                 M+3,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN('VARIANCE OF ORBITAL VELOCITY XX: LARGE-SCALE= ',
                                                                  INTEG:9:8);

                     INTEG_SMALL:= TWO_INTEGRALN(B,LIMIT_H,C,D, E,
                                                 M+3,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN('VARIANCE OF ORB. VEL. XX: SMALL-SCALE= ',
                                                            INTEG_SMALL:9:8);
                     IF LOG_LARGE
                     THEN  SIGMATTX:= INTEG
                     ELSE  SIGMATTX:= INTEG + INTEG_SMALL;
                     WRITELN('              TOTAL VARIANCE OF ORBITAL ',
                             'VELOCITY XX= ',(INTEG+INTEG_SMALL):9:8);
{***** ëäéêéëíà *****
 ********************
}
              {TT}   INTEG:= TWO_INTEGRALN(A,B,C,D, E,
                                            M+4,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN('VARIANCE OF ORBITAL VELOCITY : LARGE-SCALE= ',
                                                                  INTEG:9:8);
                     INTEG_SMALL:= TWO_INTEGRALN(B,LIMIT_H,C,D, E,
                                            M+4,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' VARIANCE OF ORB. VEL.: SMALL-SCALE = ',
                                                            INTEG_SMALL:9:8);
                     IF LOG_LARGE
                     THEN  SIGMATT:= INTEG
                     ELSE  SIGMATT:= INTEG + INTEG_SMALL;
                     WRITELN('              TOTAL VARIANCE OF ORB. VEL. = ',
                                                   (INTEG+INTEG_SMALL):9:8);

              {XT}   INTEG:= TWO_INTEGRALN(A,B,C,D, E,
                                                M+5,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' KZXT CORRELATION : LARGE-SCALE = ',INTEG:9:8);

                     INTEG_SMALL:= TWO_INTEGRALN(B,LIMIT_H,C,D, E,
                                                 M+5,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' KZXT CORRELATION: SMALL-SCALE= ',
                                                          INTEG_SMALL:9:8);

                     IF LOG_LARGE
                     THEN  KZXT:= INTEG
                     ELSE  KZXT:= INTEG + INTEG_SMALL;
                     WRITELN('              TOTAL KZXT CORRELATION= ',
                                                   (INTEG+INTEG_SMALL):9:8);

             {YT}    INTEG:= TWO_INTEGRALN(A,B,C,D, E,
                                                 M+6,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' KZYT CORRELATION: LARGE-SCALE = ',INTEG:9:8);

                     INTEG_SMALL:= TWO_INTEGRALN(B,LIMIT_H,C,D, E,
                                            M+6,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' KZYT CORRELATION: SMALL-SCALE= ',
                                                            INTEG_SMALL:9:8);

                     IF LOG_LARGE
                     THEN  KZYT:= INTEG
                     ELSE  KZYT:= INTEG + INTEG_SMALL;
                     WRITELN('              TOTAL KZYT CORRELATION= ',
                                                    (INTEG+INTEG_SMALL):9:8);

             {XY}    INTEG:= TWO_INTEGRALN(A,B,C,D, E,
                                            M+7,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' KZXY CORRELATION: LARGE-SCALE= ',INTEG:9:8);

                     INTEG_SMALL:= TWO_INTEGRALN(B,LIMIT_H,C,D, E,
                                            M+6,NG1,N1X,N2X,N1Y,N2Y);
                     WRITELN(' KZYT CORRELATION: SMALL-SCALE= ',
                                                            INTEG_SMALL:9:8);

                     IF LOG_LARGE
                     THEN  KZXY:= INTEG
                     ELSE  KZXY:= INTEG + INTEG_SMALL;
                     WRITELN('              TOTAL KZXY CORRELATION= ',
                                                    (INTEG+INTEG_SMALL):9:8);

                     WRITELN(' ----------------------------------------- ');
                     WRITELN(' SIGMAXX =  ',SIGMAXX:9:8);
                     WRITELN(' SIGMAYY =  ',SIGMAYY:9:8);
                     WRITELN(' SIGMATTX=  ',SIGMATTX:9:8);
                     WRITELN(' SIGMATT =  ',SIGMATT:9:8);
                     WRITELN(' KZXT    =  ',KZXT:9:8);
                     WRITELN(' KZYT    =  ',KZYT:9:8);
                     WRITELN(' KZXY    =  ',KZXY:9:8);
                     WRITELN(' U10     =  ',U10T:4:2,'  FT ',F0*180/PI:5:2);
                     WRITELN(' ------------------------------------------ ');

                     WRITELN(F,'  ',U10T:4:2,'   ',F0*180/PI:4:2,'   ',
                                       SIGMAXX:7:6,'    ',SIGMAYY:7:6,'    ',
                                       SIGMATTX:7:6,'    ',
                                       SIGMATT:7:6,'    ',KZXT:7:6,'   ',
                                       KZYT:7:6,'    ',KZXY:7:6);

			 	     F0:= F0 + HW;
		           END;
				WRITELN(' END OF CALCULATIONS ')
		      END;


      CLOSE(F);

 END.
