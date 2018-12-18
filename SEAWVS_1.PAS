			                PROGRAM  SEAWVS_1;
{      PROGRAM CALCULATES WAVENUMBER SPECTRUM OF SEA WAVES

                   МОДЕРНИЗАЦИЯ
       1) HEIGHT AND SLOPES - LONG/SHORT CALCULATION
       2) JONSWAP + Elfouhaily
       3) DEVELOPING WIND WAVES
       4) SLICK
}
   USES
		GREAD_V,USEFUL_V,CRT,UCHAIN_V;

   VAR
      ALPHA,KM,WM,U10,K0,KK,HK,KT,G,AU,BIG,GAMMA : EXTENDED;
      UFRICT, LIMIT_1,LIMIT_2, LIMIT_3,LIMIT_4,WT,D : EXTENDED;
      LIMIT_1K,LIMIT_2K,LIMIT_3K,LIMIT_4K : EXTENDED;
      MAX_SWH,X_TILD, SLICK_KONTRAST : EXTENDED;
      BETA0,SIGMA_W,SIGMA_M,EMASLO, UFTR : EXTENDED;
      M1: INTEGER;
      FNAME : STRING;
      F : TEXT;
      LOG_LONG,LOG_SLICK :  BOOLEAN;

{  ***** НОВЫЕ ПЕРЕМЕННЫЕ ДЛЯ СПЕКТРА CHAPRON *****
}
       OMEGA_C,OMEGA_P,K_P,C_P,GAMMA_CH,U19 : EXTENDED;
       ALPHA_P,ALPHA_M,C_M,K_M,OMEGA_M,UFRICTION : EXTENDED;
       LOG_TILTS : BOOLEAN;

{
----------------------------------------------------------------------------
}
           FUNCTION   FFF(X,Z : EXTENDED) : EXTENDED;
{
         ДЛЯ   ПОИСКА   НУЛЕЙ, X - СКОРОСТЬ ТРЕНИЯ, Z - ВЫСОТА ,
	     ПОСТОЯННАЯ КАРМАНА 0.4
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
{          ФУНКЦИЯ НАХОДИТ СКОРОСТЬ ВЕТРА НА НУЖНОЙ ВЫСОТЕ ZEND
           ПО ВВЕДЕННОЙ ВЫСОТЕ И СКОРОСТИ ВЕТРА. КРОМЕ ТОГО ВЫДАЕТ
           ИНФОРМАЦИЮ О СКОРОСТИ ТРЕНИЯ
}
    VAR
		 YWINDOW,XWINDOW : BYTE;
{
***** НАЧАЛО ВНУТРЕННЕЙ ФУНКЦИИ *****
*************************************
}
                 FUNCTION  FIND_FRICTION(X0,HX,EPSX: EXTENDED) : EXTENDED;
                    { НАХОДИТ  НУЛИ  ФУНКЦИИ. В  ДАННОМ  СЛУЧАЕ
                      НАХОДИТ СКОРОСТЬ ТРЕНИЯ  }
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

{**** КОНЕЦ ВСПОМОГАТЕЛЬНОЙ ФУНКЦИИ ****
****************************************
}
	BEGIN     {**** ГЛАВНАЯ ФУНКЦИЯ ****}
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
{         НАХОДИТ  РЕШЕНИЕ УРАВНЕНИЯ ТРЕТЕЙ СТЕПЕНИ
          БЕЗ КОМПЛЕКСНЫХ КОРНЕЙ
          WX - ЧАСТОТА
          P - КОЭФФИЦИЕНТ ПРИ ПЕРВОЙ СТЕПЕНИ
          Q - СВОБОДНЫЙ ЧЛЕН
          ВЫЧИСЛЯЕТСЯ К ИЗ ДИСПЕРСИОНКИ
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
{             ФУНКЦИЯ ВЫПОЛНЯЕТ ПЕРЕСЧЕТ ВОЛНОВОГО ЧИСЛА
         В OMEGA ПО ПОЛНОЙ ДИСПЕРСИОНКЕ
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
{        ЗАВИСИМОСТЬ  ПРЕДЕЛА ОТ СКОРОСТИ ВЕТРА
                   НОВЫЙ ВАРИАНТ  ОТ 14/7-99
}
    BEGIN
         IF (U10 < 2.5) OR (U10 > 22.0)
         THEN BEGIN
                 WRITELN(' Wind speed is not correct ');
                 WRITELN(' U10 = [2.5 - 22M/S] ');
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
{        ЗАВИСИМОСТЬ  деления ПО К ВОЛНЕНИЯ
          НА КРУПНОМАСШТАБНОЕ И РЯБЬ НОВЫЙ ВАРИАНТ
}
    VAR
         Y : EXTENDED;
         NLS : INTEGER;

    BEGIN
         IF (U10 < 2.0) OR (U10 > 23.0)
         THEN BEGIN
                 WRITELN(' wind speed is not correct ');
                 WRITELN(' ФУНКЦИЯ KM_LIMIT ');
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
                  WRITE(' program is working for 0.08 m, 0.021m, 0.03 m ');
                  WRITELN(' 0.055 m, 0.1 m, 0.23 m ');
                  WRITELN(' function KL_LIMIT ');
                  HALT;
               END;
         END;

         KL_LIMIT:= Y;
    END;

{
 --------------------------------------------------------------------------
}
           FUNCTION  DET(X : EXTENDED) : EXTENDED;
{             ФУНКЦИЯ РАВНА ДЕТЕРМИНАНТУ ПРИ
              ПЕРЕХОДЕ ОТ ЧАСТОТЫ К ВОЛНОВЫМ ЧИСЛАМ
              ПО ПОЛНОЙ ДИСПЕРСИОНКЕ
              X - ВОЛНОВОЕ ЧИСЛО
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
{         ВЫЧИСЛЯЕТСЯ OMEGA~  ЗАВИСЯЩЕЕ ОТ РАЗГОНА X
          ЭТО БЕЗРАЗМЕРНЫЙ РАЗГОН
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
{         ВЫЧИСЛЯЕТСЯ GAMMA~  ЗАВИСЯЩЕЕ ОТ РАЗГОНА X
          ЭТО БЕЗРАЗМЕРНЫЙ РАЗГОН
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
{         ВЫЧИСЛЯЕТСЯ  ALPHA~ ЗАВИСЯЩАЯ ОТ РАЗГОНА X
          ЭТО БЕЗРАЗМЕРНЫЙ РАЗГОН
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
 ------------------------------------------------------------------------
}
           FUNCTION  SPECTR_CHAPRON(X : EXTENDED) : EXTENDED;
{                    СПЕКТР  CHAPRON
        X - ВОЛНОВОЕ ЧИСЛО
        ПРИ ИНТЕГРИРОВАНИИ НЕ НУЖНО ДОМНОЖАТЬ НА ВОЛНОВОЕ
        ЧИСЛО, Т.К. КdK
}
    VAR
        SPECTR, BL,BH,C,W,FP,LPM,JP,GAMMA,SIGMA,FM : EXTENDED;

    BEGIN
          W:= OMEGA_K(X);  {ТЕКУЩАЯ ЧАСТОТА}
          C:= W / X;  {ТЕКУЩАЯ ФАЗОВАЯ СКОРОСТЬ}
          SIGMA:= 0.08 * ( 1 + 4.0*POWER(OMEGA_C,-3.0) );
          GAMMA:= EXPO(- SQR( SQRT(X/K_P)-1.0 )/( 2.0*SQR(SIGMA) ));

          LPM:= EXPO(-1.25*SQR(K_P/X));
          JP:= POWER(GAMMA_CH,GAMMA);
          FP:= LPM * JP * EXPO(- OMEGA_P/SQRT(10.0) * ( SQRT(X/K_P) - 1.0 ));
          BL:= 0.5 * ALPHA_P * C_P/C * FP;

          FM:= EXPO(-0.25*SQR(X/K_M - 1.0));
          BH:= 0.5 * ALPHA_M * C_M/C * FM   * LPM*JP;
                                            { ДОБАВКА }
          SPECTR_CHAPRON:= (BL + BH) / (X*X*X);
    END;

{
 -------------------------------------------------------------------------
}
		   FUNCTION   SPECTR_WAVE(X : EXTENDED) : EXTENDED;
{               СПЕКТР JONSWAP
                X - ВОЛНОВОЕ ЧИСЛО
                ПРИ ИНТЕГРИРОВАНИИ НЕ НУЖНО ДОМНОЖАТЬ НА ВОЛНОВОЕ
                ЧИСЛО, Т.Е. KdK
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
{                 СПЕКТР ВОЛНЕНИЯ ПРЕДЛОЖЕННЫЙ КИТАЙГОРОДСКИМ
        НАЧИНАЕТСЯ  ОТ 1.2WM ДО AM WM
        W - ЧАСТОТА
        ! НЕТ КОЭФФИЦИЕНТА BETA, ОН ИЩЕТСЯ СТЫКОВКОЙ
        СПЕКТРОВ В ТОЧКЕ LIMIT_1K, 1/W**4
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
{                 СПЕКТР ВОЛНЕНИЯ НАША МОДЕЛЬ
        НАЧИНАЕТСЯ  ОТ AM WM ДО WGC
        W - ЧАСТОТА
        ! НЕТ КОЭФФИЦИЕНТА BETA, ОН ИЩЕТСЯ СТЫКОВКОЙ
        СПЕКТРОВ В ТОЧКЕ LIMIT_2K         1/W**5
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
{                 СПЕКТР ВОЛНЕНИЯ
        НАЧИНАЕТСЯ  ОТ WGC  ДО WC НАША МОДЕЛЬ
        W - ЧАСТОТА
        ! НЕТ КОЭФФИЦИЕНТА BETA, ОН ИЩЕТСЯ СТЫКОВКОЙ
        СПЕКТРОВ В ТОЧКЕ LIMIT_3K
        W**2.7
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
{                 СПЕКТР ВОЛНЕНИЯ
        НАЧИНАЕТСЯ  ОТ WC   НАША МОДЕЛЬ
        W - ЧАСТОТА
        ! НЕТ КОЭФФИЦИЕНТА BETA, ОН ИЩЕТСЯ СТЫКОВКОЙ
        СПЕКТРОВ В ТОЧКЕ LIMIT_3K
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
 -----------------------------------------------------------------------
 }
       FUNCTION GAMMA_KONTRAST(X,E,OMEGA : EXTENDED) : EXTENDED;
{  X - ВОЛНОВОЕ ЧИСЛО, BETA - ИНКРЕМЕНТ РАЗГОНА
   SIGMA - ВЯЗКОСТЬ - 72 - ВОДА, 40 МАСЛО  E - УПРУГОСТЬ ПЛЕНКИ = 15-20 }
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
{         ПОЛНЫЙ СПЕКТР - JONSWAP + ДРУГИЕ СПЕКТРЫ
        СТЫКОВКА СПЕКТРОВ - W/WM
        ДО LIMIT_1K - 1/W**5
     LIMIT_2K    WM - 1/(W**5)
     LIMIT_3K    1/(W**2.7)
     LIMIT_4K    1/(W**5)
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
{             УПРАВЛЯЮЩАЯ ФУНКЦИЯ ДЛЯ ПОСТРОЕНИЯ СПЕКТРА
              K - ВОЛНОВОЕ ЧИСЛО
}
   VAR
      KONTRASTX : EXTENDED;
   BEGIN

         CASE  M OF
           1:   WAVE_SPECTRUM:= SPECTR_WAVE(K);
           2:   BEGIN
                  IF LOG_SLICK
                  THEN BEGIN
                        KONTRASTX:= KONTRAST(K/100,BETA0,SIGMA_W,SIGMA_M,
                                             EMASLO,UFTR);
                        WAVE_SPECTRUM:= FULL_SPECTR2(K)*KONTRASTX;
                                          { СУММАРНЫЙ СПЕКТР }
                       END
                  ELSE BEGIN
                        WAVE_SPECTRUM:= FULL_SPECTR2(K); { СУММАРНЫЙ СПЕКТР }
                       END;
                END;
           3:   WAVE_SPECTRUM:= SPECTR_CHAPRON(K); { СПЕКТР CHAPRON }
         ELSE  BEGIN
                  WRITELN('       THIS NUMBER IS ABSENT  ');
                  WRITELN('       FUNCTION     Wave_Spectrum ');
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
     WRITELN('                     NEW WAVE MODEL                        ');
     WRITELN(' ==========================================================');

{***** INITIAL SET *****
 *******************************}
     G:= 9.8;
     SLICK_KONTRAST:= 0;
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


{***** ГРАНИЦЫ СПЕКТРА *****
 ***************************}
     LIMIT_1:= 1.2;   { НАЧИНАЕТСЯ  - 1/(W**4) }
     LIMIT_3:= 270.0; { НАЧИНАЕТСЯ 1/W**2.7 }
     LIMIT_4:= 1020.0; { НАЧИНАЕТСЯ 1/W*5 }

     LOG_LONG:= ANT(' LONG CALCULATIONS ? ');
     LOG_TILTS:= ANT(' CURVES (Y) OR HEIGHT (N) ? ');
     WRITELN(' THE NEXT FUNCTIONS : ');
     WRITELN('1 - JONSWAP (K), 2 - NEW MODEL (K), 3 - CHAPRON (K) ');
     WRITELN(' -------------------------------------------------------- ');

     REPEAT
   		  U10:= RREAD(' Input U10,  M/S  ');
          U19:= WIND(U10*100.0,1000.0,1950.0, UFRICTION)/100.0;
          UFRICTION:= UFRICTION / 100.0;
          UFTR:= UFRICTION*100;
          WM:= OMEGA_M*G/U10;
          KM:= SQR(WM)/G;
          WRITELN;

          WRITELN('U10= ',U10:5:2,'M/S   U* friction= ',UFRICTION:4:2,' M/S');
          WRITELN('KM = ',KM:7:5,'  WM= ',WM:7:5,'     LM= ',2.0*PI/KM:5:2);
          WRITELN;

          LIMIT_1K:=  FIND_DECISION(LIMIT_1*WM);

          LIMIT_2:= KM_LIMIT(U10);
          LIMIT_2K:=  FIND_DECISION(LIMIT_2*WM);
          LIMIT_3K:=  LIMIT_3;
          LIMIT_4K:=  LIMIT_4;

{***** ПАРАМЕТРЫ ДЛЯ СПЕКТРА CHAPRON *****
 *****************************************}
          C_P:= WM / KM;  { ФАЗОВАЯ СКОРОСЬ ПИКА }
          OMEGA_P:= U10 / C_P;  { РАБОТАЕМ С РАЗВИТЫМ ВОЛНЕНИЕМ }
          OMEGA_C:= U10 / C_P;  { ПОКАЗАТЕЛЬ СТЕПЕНИ РАЗВИТИЯ }
          C_M:= 0.23;
          K_M:= 2.0*G/SQR(0.23);
          OMEGA_M:= UFRICTION / C_M;
          IF UFRICTION < C_M
          THEN ALPHA_M:= 0.01 * ( 1.0 + LN(UFRICTION/C_M) )
          ELSE ALPHA_M:= 0.01 * ( 1.0 + 3.0*LN(UFRICTION/C_M) );

          IF (0.83 <= OMEGA_C) OR (OMEGA_C <= 5.0)
          THEN  BEGIN
                     ALPHA_P:= 0.006 * POWER(OMEGA_C,0.55);
                     K_P:= KM;

                     IF (OMEGA_C <= 1.0)
                     THEN  GAMMA_CH:= 1.7
                     ELSE  GAMMA_CH:= 1.7 + 6.0*LOG10(OMEGA_C)
                END
          ELSE  BEGIN
                     WRITELN(' ВОЛНЕНИЕ ВЫХОДИТ ЗА РАМКИ 0.83 - 5.0 ');
                     DELAY(2000);
                     HALT
                END;

          WRITELN;
          K0:= RREAD(' Input  K0 ');
          KK:= RREAD(' Input  KK ');
          HK:= RREAD(' Input  HK ');

          FNAME:= 'SEAWVS_1.DAT';
          ASSIGN(F,FNAME);
          IF NOT ANT(' FILE IS REWRITTEN ? ')
          THEN HALT;
          REWRITE(F);

          M1:= 1;

          WRITELN(F,'   KT  ','         JONSWAP  ','        OUR ',
                     '          ELFOUHAILY             WT  ');

          KT:= K0;
          WRITELN(' FRICTION VELOCITY, M/S ',UFTR:4:2);

          WHILE  KT <= KK
          DO BEGIN
               WT:= OMEGA_K(KT);

               IF LOG_TILTS
               THEN BEGIN
                       WRITELN('  ',KT:7:4,'   ',
                                    KT*KT*KT*WAVE_SPECTRUM(M1,KT):15:12,
                               '   ',KT*KT*KT*WAVE_SPECTRUM(2,KT):15:12,
                               '   ',KT*KT*KT*WAVE_SPECTRUM(3,KT):15:12,
                               '   ',WT:7:4);

                       WRITELN(F,'  ',KT:7:4,
                                 '    ',KT*KT*KT*WAVE_SPECTRUM(M1,KT):12,
                                 '    ',KT*KT*KT*WAVE_SPECTRUM(2,KT):12,
                                 '    ',KT*KT*KT*WAVE_SPECTRUM(3,KT):12,
                                 '    ',WT:7:4);
                    END
               ELSE BEGIN
                       WRITELN(KT:6:3,'   ',WAVE_SPECTRUM(M1,KT):15:12,
                                      '   ',WAVE_SPECTRUM(2,KT):15:12,
                                      '   ',WAVE_SPECTRUM(3,KT):15:12,
                                      '   ',WT:7:4);

                       WRITELN(F,'  ',KT:7:4,
                                 '    ',WAVE_SPECTRUM(M1,KT):12,
                                 '    ',WAVE_SPECTRUM(2,KT):12,
                                 '    ',WAVE_SPECTRUM(3,KT):12,
                                 '    ',WT:7:4);
                    END;

               IF LOG_LONG
               THEN KT:= KT + KT*HK
               ELSE KT:= KT + HK;
             END;

          CLOSE(F);
	 UNTIL NOT ANT(' CONTINUE ? ');
   END.