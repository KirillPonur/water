			                PROGRAM  SEAWVS_2;
{      PROGRAM CALCULATES WAVENUMVER SPECTRUM OF Swell

                   ŒŽ„…ˆ‡€–ˆŸ
       1) HEIGHT AND SLOPES - LONG CALCULATION
       2) SWELL
       3) WIND WAVES
       4) MIXED SEA
}
   USES
		GREAD_V,USEFUL_V,CRT,UCHAIN_V;

   VAR
      ALPHA,KM,WM,U10,K0,KK,HK,KT,G,AU,BIG,GAMMA : EXTENDED;
      UFRICT, LIMIT_1,LIMIT_2, LIMIT_3,LIMIT_4,WT,D : EXTENDED;
      LIMIT_1K,LIMIT_2K,LIMIT_3K,LIMIT_4K : EXTENDED;
      MAX_SWH,X_TILD,OMEGA_M,L_SWELL,KM_SWELL,WM_SWELL,ALPHA_SWELL: EXTENDED;
      U10_SWELL,U19_SWELL,UFRICTION,U19,SP_SWELL,SP_WAVE : EXTENDED;
      M1: INTEGER;
      FNAME : STRING;
      F : TEXT;
      LOG_LONG,LOG_SWELL,LOG_TILTS :  BOOLEAN;

{
----------------------------------------------------------------------------
}
           FUNCTION   FFF(X,Z : EXTENDED) : EXTENDED;
{
         „‹Ÿ   Žˆ‘Š€   “‹…‰, X - ‘ŠŽŽ‘’œ ’…ˆŸ, Z - ‚›‘Ž’€ ,
	     Ž‘’ŽŸ€Ÿ Š€Œ€€ 0.4
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
{          ”“Š–ˆŸ €•Ž„ˆ’ ‘ŠŽŽ‘’œ ‚…’€ € “†Ž‰ ‚›‘Ž’… ZEND
           Ž ‚‚…„…Ž‰ ‚›‘Ž’… ˆ ‘ŠŽŽ‘’ˆ ‚…’€. ŠŽŒ… ’ŽƒŽ ‚›„€…’
           ˆ”ŽŒ€–ˆž Ž ‘ŠŽŽ‘’ˆ ’…ˆŸ
}
    VAR
		 YWINDOW,XWINDOW : BYTE;
{
***** €—€‹Ž ‚“’……‰ ”“Š–ˆˆ *****
*************************************
}
                 FUNCTION  FIND_FRICTION(X0,HX,EPSX: EXTENDED) : EXTENDED;
                    { €•Ž„ˆ’  “‹ˆ  ”“Š–ˆˆ. ‚  „€ŽŒ  ‘‹“—€…
                      €•Ž„ˆ’ ‘ŠŽŽ‘’œ ’…ˆŸ  }
		         VAR XA,XB,X : EXTENDED;
			     YA,YB,Y : BOOLEAN;
		         BEGIN
			 	      XA:= X0;
				      XB:= X0 + HX;

				      YA:= (UBEG - FFF(XA,ZBEG)) > 0.0;
				      YB:= (UBEG - FFF(XB,ZBEG)) > 0.0;

				      WHILE  ( (YA AND YB) OR (NOT YA AND NOT YB) )
				      DO BEGIN
						   XA:= XB;
						   YA:= YB;
						   XB:= XA + HX;
						   YB:= (UBEG - FFF(XB,ZBEG)) > 0.0;
					     END;

					  WHILE  (XB-XA) > EPSX
					  DO BEGIN
						   X:= (XA + XB) / 2.0;
						   Y:= (UBEG - FFF(X,ZBEG)) > 0.0;

						   IF  (YA AND Y) OR (NOT YA AND NOT Y)
						   THEN BEGIN
								  XA:= X;
								  YA:= Y
								END
						   ELSE BEGIN
								  XB:= X;
								  YB:= Y;
								END;
						 END;
					  FIND_FRICTION:= (XA+XB) / 2.0;

		         END;

{**** ŠŽ…– ‚‘ŽŒŽƒ€’…‹œŽ‰ ”“Š–ˆˆ ****
****************************************
}
	BEGIN     {**** ƒ‹€‚€Ÿ ”“Š–ˆŸ ****}
        YWINDOW:= WHEREY;
        WRITELN;
        WRITE_NOVISIBLE('  WIND SPEED (HEIGHT OF ANEMOMETER) - IN CM/S (CM) ',
                                                                       1200);
        GOTOXY(30,YWINDOW-1);

        IF  UBEG > 8000.0
        THEN BEGIN
                 WRITELN(' VERY STRONG WIND SPEED ');
                 READLN;
                 HALT;
             END;
	    IF  UBEG < 40.0
        THEN BEGIN
                 WRITELN(' VERY LOW WIND SPEED ');
                 READLN;
                 HALT;
             END;

         UFRICTION:= FIND_FRICTION(0.0001,80.0,0.001);
         WIND:= FFF(UFRICTION,ZEND);
    END;

{
 -------------------------------------------------------------------------
}
       FUNCTION  FIND_DECISION(WX: EXTENDED) : EXTENDED;
{         €•Ž„ˆ’  …˜…ˆ… “€‚…ˆŸ ’…’…‰ ‘’……ˆ
          …‡ ŠŽŒ‹…Š‘›• ŠŽ…‰
          WX - —€‘’Ž’€
          P - ŠŽ””ˆ–ˆ…’ ˆ …‚Ž‰ ‘’……ˆ
          Q - ‘‚ŽŽ„›‰ —‹…
          ‚›—ˆ‘‹Ÿ…’‘Ÿ Š ˆ‡ „ˆ‘…‘ˆŽŠˆ
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
{             ”“Š–ˆŸ ‚›Ž‹Ÿ…’ ……‘—…’ ‚Ž‹Ž‚ŽƒŽ —ˆ‘‹€
         ‚ OMEGA Ž Ž‹Ž‰ „ˆ‘…‘ˆŽŠ…
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
{        ‡€‚ˆ‘ˆŒŽ‘’œ  …„…‹€ Ž’ ‘ŠŽŽ‘’ˆ ‚…’€
}
    BEGIN
         IF (U10 < 2.5) OR (U10 > 22.0)
         THEN BEGIN
                 WRITELN(' WIND SPEED IN NOT CORRECT ');
                 WRITELN(' FUNCTION KM_LIMIT, U10 = [2.5 - 22M/S] ');
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
{        ‡€‚ˆ‘ˆŒŽ‘’œ  ¤¥«¥­¨ï Ž Š ‚Ž‹…ˆŸ
          € Š“ŽŒ€‘˜’€Ž… ˆ Ÿœ Ž‚›‰ ‚€ˆ€’
}
    VAR
         Y : EXTENDED;
         NLS : INTEGER;

    BEGIN
         IF (U10 < 2.0) OR (U10 > 23.0)
         THEN BEGIN
                 WRITELN(' WIND SPEED IS NOT CORRECT ');
                 WRITELN(' FUNCTION KM_LIMIT ');
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
                  WRITELN('THIS RADAR WAVELENGTH IS ABSENT ');
                  WRITELN('L = 0.008,0.021,0.03,0.055,0.1,0.23 M ');
                  HALT;
               END;
         END;

         KL_LIMIT:= Y;
    END;


{
 --------------------------------------------------------------------------
}
           FUNCTION  DET(X : EXTENDED) : EXTENDED;
{             ”“Š–ˆŸ €‚€ „…’…Œˆ€’“ ˆ
              ……•Ž„… Ž’ —€‘’Ž’› Š ‚Ž‹Ž‚›Œ —ˆ‘‹€Œ
              Ž Ž‹Ž‰ „ˆ‘…‘ˆŽŠ…
              X - ‚Ž‹Ž‚Ž… —ˆ‘‹Ž
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
{         ‚›—ˆ‘‹Ÿ…’‘Ÿ OMEGA~  ‡€‚ˆ‘Ÿ™…… Ž’ €‡ƒŽ€ X
          ’Ž …‡€‡Œ…›‰ €‡ƒŽ
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
{         ‚›—ˆ‘‹Ÿ…’‘Ÿ GAMMA~  ‡€‚ˆ‘Ÿ™…… Ž’ €‡ƒŽ€ X
          ’Ž …‡€‡Œ…›‰ €‡ƒŽ
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
{         ‚›—ˆ‘‹Ÿ…’‘Ÿ  ALPHA~ ‡€‚ˆ‘Ÿ™€Ÿ Ž’ €‡ƒŽ€ X
          ’Ž …‡€‡Œ…›‰ €‡ƒŽ
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
				   ******* ‘…Š’ ‚Ž‹…ˆŸ ********
}
           FUNCTION  SPECTR_SWELL(X : EXTENDED) : EXTENDED;
{                ‘…Š’  „€‚ˆ„€€, X - ‚Ž‹Ž‚Ž… —ˆ‘‹Ž
                Ž —€‘’Ž’›‰ ‚€ˆ€’
}
    VAR
        W,M0,ALPHA_W : EXTENDED;

    BEGIN
          W:= OMEGA_K(X);  {’…Š“™€Ÿ —€‘’Ž’€}

          M0:= SQR(MAX_SWH)/16.0;
          SPECTR_SWELL:= 6.0 * M0 * POWER(WM_SWELL/W,5.0) *
                          expo(-1.2*POWER(WM_SWELL/W,5.0))/ W
                           *SQRT(9.8/X)/2.0; {ŠŽ…Š–ˆŸ }
    END;


{
 -------------------------------------------------------------------------
}
		   FUNCTION   SPECTR_WAVE(X : EXTENDED) : EXTENDED;
{               ‘…Š’ JONSWAP
                X - ‚Ž‹Ž‚Ž… —ˆ‘‹Ž
                ˆ ˆ’…ƒˆŽ‚€ˆˆ … “†Ž „ŽŒŽ†€’œ € ‚Ž‹Ž‚Ž…
                —ˆ‘‹Ž, ’.…. KdK
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
{                 ‘…Š’ ‚Ž‹…ˆŸ …„‹Ž†…›‰ Šˆ’€‰ƒŽŽ„‘ŠˆŒ
        €—ˆ€…’‘Ÿ  Ž’ 1.5 „Ž 3.5 WM
        W - —€‘’Ž’€
        ! …’ ŠŽ””ˆ–ˆ…’€ BETA, Ž ˆ™…’‘Ÿ ‘’›ŠŽ‚ŠŽ‰
        ‘…Š’Ž‚ ‚ ’Ž—Š… LIMIT_1K, 1/W**4
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
{                 ‘…Š’ ‚Ž‹…ˆŸ €˜€ ŒŽ„…‹œ
        €—ˆ€…’‘Ÿ  Ž’ 3.5 „Ž 10.5 WM
        W - —€‘’Ž’€
        ! …’ ŠŽ””ˆ–ˆ…’€ BETA, Ž ˆ™…’‘Ÿ ‘’›ŠŽ‚ŠŽ‰
        ‘…Š’Ž‚ ‚ ’Ž—Š… LIMIT_2K         1/W**5
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
{                 ‘…Š’ ‚Ž‹…ˆŸ
        €—ˆ€…’‘Ÿ  Ž’ 30  WM  „Ž 100 €˜€ ŒŽ„…‹œ
        W - —€‘’Ž’€
        ! …’ ŠŽ””ˆ–ˆ…’€ BETA, Ž ˆ™…’‘Ÿ ‘’›ŠŽ‚ŠŽ‰
        ‘…Š’Ž‚ ‚ ’Ž—Š… LIMIT_3K
        W**3
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
{                 ‘…Š’ ‚Ž‹…ˆŸ
        €—ˆ€…’‘Ÿ  Ž’ 30  WM  „Ž 100 €˜€ ŒŽ„…‹œ
        W - —€‘’Ž’€
        ! …’ ŠŽ””ˆ–ˆ…’€ BETA, Ž ˆ™…’‘Ÿ ‘’›ŠŽ‚ŠŽ‰
        ‘…Š’Ž‚ ‚ ’Ž—Š… LIMIT_3K
        W**5
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
----------------------------------------------------------------------------
}
          FUNCTION   FULL_SPECTR2(KT : EXTENDED) : EXTENDED;
{         Ž‹›‰ ‘…Š’ - JONSWAP + „“ƒˆ… ‘…Š’›
        ‘’›ŠŽ‚Š€ ‘…Š’Ž‚ - W/WM
        „Ž LIMIT_1K - 1/W**5
     LIMIT_2K    WM - 1/(W**4)
     LIMIT_3K    1/(W**5)
     LIMIT_4K    1/(W**3.5)
}
   BEGIN
      IF KT <= LIMIT_1K
      THEN FULL_SPECTR2:=  SPECTR_WAVE(KT)
      ELSE IF KT <= LIMIT_2K
           THEN FULL_SPECTR2:=  NEW_SPECTR1(KT)
           ELSE IF KT <= LIMIT_3K
                THEN  FULL_SPECTR2:= NEW_SPECTR2(KT)
                ELSE  IF KT <= LIMIT_4K
                      THEN FULL_SPECTR2:= NEW_SPECTR3(KT)
                      ELSE FULL_SPECTR2:= NEW_SPECTR4(KT);
   END;


{
----------------------------------------------------------------------------
}
          FUNCTION   WAVE_SPECTRUM(M : INTEGER; K : EXTENDED) : EXTENDED;
{             “€‚‹Ÿž™€Ÿ ”“Š–ˆŸ „‹Ÿ Ž‘’Ž…ˆŸ ‘…Š’€
              K - ‚Ž‹Ž‚Ž… —ˆ‘‹Ž
}
   BEGIN

         CASE  M OF
           1:   WAVE_SPECTRUM:= SPECTR_WAVE(K);
           2:   WAVE_SPECTRUM:= FULL_SPECTR2(K); { ‘“ŒŒ€›‰ ‘…Š’ }
                                            { Ž‚›‰ ‚€ˆ€’ 14 ˆž‹Ÿ 99}
           3:   WAVE_SPECTRUM:= SPECTR_SWELL(K); { SWELL }
         ELSE  BEGIN
                  WRITELN('     THIS NUMBER IS ABSENT ');
                  WRITELN('     FUNCTION     Wave_Spectrum ');
                  HALT;
               END;
         END;
   END;

{
  ========================================================================
                                 MAIN PROGRAM
}
   BEGIN

     WRITELN(' ----------------------------------------------------------');
     WRITELN('                     NEW WAVE MODEL (SWELL)                ');
     WRITELN(' ==========================================================');

{***** INITIAL SET *****
 *******************************}
     LOG_SWELL:= ANTWRT(' CALCULATION OF SWELL (Y) OR MIXED SEA (N) ? ');

     G:= 9.8;

     IF NOT LOG_SWELL
     THEN BEGIN
            X_TILD:= RREAD(' INPUT WIND FETCH (1430-20170) ');
            GAMMA:= GAMMA_DEV(X_TILD);
            ALPHA:= ALPHA_DEV(X_TILD);
            OMEGA_M:= OMEGA_DEV(X_TILD);
          END;

{***** ƒ€ˆ–› ‘…Š’€ *****
 ***************************}
     LIMIT_1:= 1.2;   { €—ˆ€…’‘Ÿ  - 1/(W**4) }
     LIMIT_3:= 270.0; { €—ˆ€…’‘Ÿ 1/W**2.7 }
     LIMIT_4:= 1020.0; { €—ˆ€…’‘Ÿ 1/W*5 }

     LOG_LONG:= ANT(' LONG CALCULATIONS ? ');
     LOG_TILTS:= ANT(' SLOPES (Y) OR HEIGHT (N) ? ');
     WRITELN(' THE NEXT FUNCTIONS ');
     WRITELN('1- JONSWAP(K), 2 - NEW MODEL, 3 - SWELL ');
     WRITELN(' -------------------------------------------------------- ');

     REPEAT
          IF LOG_SWELL
          THEN BEGIN
                 L_SWELL:= RREAD(' Input wavelenght of swell ');
                 KM_SWELL:= 2*PI/L_SWELL;
                 WM_SWELL:= SQRT(G*KM_SWELL);
                 ALPHA_SWELL:= 0.0081;

                 MAX_SWH:= ALPHA_SWELL*SQR(G)*EXPO(-0.05)
                            /( 6.0*POWER(WM_SWELL,4.0) );
                 WRITELN(' MAX SWH = ',4.0*SQRT(MAX_SWH):4:3);
                 WRITELN(' M0 = ',MAX_SWH:5:4);
                 MAX_SWH:= RREAD(' Input SWH OF SWELL ');

                 U10_SWELL:= SQRT(0.697*G/KM_SWELL);
                 U19_SWELL:= WIND(U10_SWELL*100.0,1000.0,1950.0,UFRICTION)
                                                                     /100.0;
                 UFRICTION:= UFRICTION / 100.0;
                 WRITELN(' U10 (SWELL) = ',U10_SWELL:4:2,' U19 (SWELL) ',
                                                           U19_SWELL:4:2);
               END
          ELSE BEGIN
                 WRITELN(' INPUT OF SWELL PARAMETERS ');
                 L_SWELL:= RREAD(' Input wavelenght of swell ');
                 KM_SWELL:= 2*PI/L_SWELL;
                 WM_SWELL:= SQRT(G*KM_SWELL);
                 ALPHA_SWELL:= 0.0081;

                 MAX_SWH:= ALPHA_SWELL*SQR(G)*EXPO(-0.05)
                            /( 6.0*POWER(WM_SWELL,4.0) );
                 WRITELN(' MAX SWH = ',4.0*SQRT(MAX_SWH):4:3);
                 WRITELN(' M0 = ',MAX_SWH:5:4);
                 MAX_SWH:= RREAD(' Input SWH OF SWELL ');

                 U10_SWELL:= SQRT(0.697*G/KM_SWELL);
                 U19_SWELL:= WIND(U10_SWELL*100.0,1000.0,1950.0,UFRICTION)
                              /100.0;
                 UFRICTION:= UFRICTION / 100.0;
                 WRITELN(' U10 (SWELL) = ',U10_SWELL:4:2,' U19 (SWELL) ',
                                                           U19_SWELL:4:2);

                 WRITELN(' INPUT WIND WAVES PARAMETERS ');
   		         U10:= RREAD(' Input U10,  Œ/‘  ');
                 U19:= WIND(U10*100.0,1000.0,1950.0, UFRICTION)/100.0;
                 UFRICTION:= UFRICTION / 100.0;
                 WM:= OMEGA_M*9.8/U10;
                 KM:= SQR(WM)/9.8;

                 WRITELN(' U10 = ',U10:4:2,' UFRICTION ',UFRICTION:4:3);
                 WRITELN(' KM = ',KM:8:6,'  WM = ',WM:8:6,
                         ' LM = ',2.0*PI/KM:6:3);
                 WRITELN;
               END;

          IF NOT LOG_SWELL
          THEN BEGIN
                 LIMIT_1K:= FIND_DECISION(LIMIT_1*WM);{€—ˆ€…’‘Ÿ - 1/(W**4) }

                 LIMIT_2:= KM_LIMIT(U10);

                 LIMIT_2K:=  FIND_DECISION(LIMIT_2*WM);  { -//-  1/(W**5) }
                 LIMIT_3K:=  LIMIT_3; { 1/W**2.7 }
                 LIMIT_4K:=  LIMIT_4; { 1/W**5.0 }

                 WRITELN(' KM LIMIT_1 = ',LIMIT_1K:6:3,
                         ' 2 = ',LIMIT_2K:6:3,
                         ' 3 = ',LIMIT_3K:6:3,
                         ' 4 = ',LIMIT_4K:6:3);

                 WRITELN;
               END;

          IF LOG_SWELL
          THEN WRITELN(' KM_SWELL = ',KM_SWELL:5:2)
          ELSE WRITELN(' KM = ',KM:5:2,' KM_SWELL= ',KM_SWELL:5:2);

          K0:= RREAD(' Input K0 ');
          KK:= RREAD(' Input KK ');
          HK:= RREAD(' Input HK ');

          FNAME:= 'SEAWVS_2.DAT';
          ASSIGN(F,FNAME);
          IF NOT ANT(' REWRITTEN FILE ? ')
          THEN HALT;
          REWRITE(F);

          WRITELN(F,'      KT  ','       SWELL   ','        WIND   ',
                    '      MIXED ');

          KT:= K0;
          WHILE  KT <= KK
          DO BEGIN
               WT:= OMEGA_K(KT);

               IF LOG_SWELL
               THEN BEGIN
                      SP_SWELL:= WAVE_SPECTRUM(3,KT);
                      SP_WAVE:= 0;
                    END
               ELSE BEGIN
                      SP_SWELL:= WAVE_SPECTRUM(3,KT);
                      SP_WAVE:= WAVE_SPECTRUM(2,KT);
                    END;

               IF LOG_TILTS
               THEN BEGIN
                       WRITELN(' ',KT:7:4,'   ',
                                    KT*KT*KT*SP_SWELL:15:12,
                               '   ',KT*KT*KT*SP_WAVE:15:12,
                               '   ',KT*KT*KT*(SP_SWELL+SP_WAVE):15:12);

                       WRITELN(F,'  ',KT:7:4,
                                 '    ',KT*KT*KT*SP_SWELL:12,
                                 '    ',KT*KT*KT*SP_WAVE:12,
                                 '    ',KT*KT*KT*(SP_SWELL+SP_WAVE):12);
                    END
               ELSE BEGIN
                       WRITELN(' ',KT:7:4,'   ',
                                    SP_SWELL:15:12,
                               '   ',SP_WAVE:15:12,
                               '   ',(SP_SWELL+SP_WAVE):15:12);

                       WRITELN(F,'  ',KT:7:4,
                                 '    ',SP_SWELL:12,
                                 '    ',SP_WAVE:12,
                                 '    ',(SP_SWELL+SP_WAVE):12);
                    END;

               IF LOG_LONG
               THEN KT:= KT + KT*HK
               ELSE KT:= KT + HK;

             END;

          CLOSE(F);
          WRITELN(' SHW SWELL = ',MAX_SWH:4:2);
	 UNTIL NOT ANT(' CONTINUE ? ');

   END.