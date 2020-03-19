

      BLOCK DATA BDPROB
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C             INITIALIZE OPTIMAL CONTROL OF ODE TEST PROBLEM
C             DATA.
C
      INCLUDE '../commons/odeprb.cmn'
C
      COMMON /OSMRYR/ FOPTNW(NTC),FOPTRF(NTC),TIM1NW(NTC),TIM1RF(NTC), 
     &                TIM2NW(NTC),TIM2RF(NTC),TIM3NW(NTC),TIM3RF(NTC),
     &                RODENW(NTC),RODERF(NTC)
      COMMON /OSMRYI/ IPRVEC(NTC),NVARNW(NTC),MCONNW(NTC),MZERNW(NTC),
     &                NFCLNW(NTC),NGCLNW(NTC),NHCLNW(NTC),NFCLRF(NTC),
     &                NGCLRF(NTC),NHCLRF(NTC),NFEVNW(NTC),IERNNW(NTC),
     &                NGRDNW(NTC),NGRDRF(NTC),NPTSNW(NTC),NPTSRF(NTC),
     &                NHNZNW(NTC),NHNZRF(NTC),NJNZNW(NTC),NJNZRF(NTC),
     &                NRHSNW(NTC),NRHSRF(NTC),IERNRF(NTC),
     &                NFEVRF(NTC),NVARRF(NTC),MCONRF(NTC),
     &                MZERRF(NTC),NPROB
      COMMON /OSMRYC/ REFDAT,REFTIM
      CHARACTER*8  REFDAT,REFTIM
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk =   1, 15) /
     a    1,  1, 1, 2, 1,   10,  0.00    , 'qM  ',  'qlin01',       
     b    2,  1, 2, 2, 1,   10,  0.00    , 'qM  ',  'qlin02',       
     c    3,  2, 1, 2, 1,   10,  0.00    , 'qF  ',  'lnts01',       
     d    4,  2, 2, 1, 4,   10,  0.00    , 'qFM ',  'lnts02',       
     e    5,  2, 2, 1, 4,   20,  0.00    , 'qFM ',  'lnts03',       
     f    6,  2, 2, 1, 4,  200,  0.00    , 'qFM ',  'lnts04',       
     g    7,  2, 2, 2, 1,   10,  0.00    , 'qFM ',  'lnts05',       
     h    8,  2, 2, 2, 1,   20,  0.00    , 'bFM ',  'lnts06',       
     i    9,  2, 2, 2, 1,  200,  0.00    , 'bFM ',  'lnts07',       
     j   10,  2, 2, 2, 1,  400,  0.00    , 'bFM ',  'lnts08',       
     k   11,  2, 2, 2, 1, 2000,  0.00    , 'bFM ',  'lnts09',       
     l   12,  2, 2, 3, 1,   10,  0.00    , 'qFM ',  'lnts10',       
     m   13,  2, 2, 3, 1,   20,  0.00    , 'qFM ',  'lnts11',       
     n   14,  2, 2, 3, 1,  200,  0.00    , 'bFM ',  'lnts12',       
     o   15,  3, 1, 1, 4,    5,  0.00    , 'bM  ',  'traj01'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk =  16, 30) /
     a   16,  3, 1, 1, 4,   20,  0.00    , 'qFM ',  'traj02',       
     b   17,  3, 1, 2, 1,    5,  0.00    , 'qFM ',  'traj03',       
     c   18,  3, 1, 2, 1,   20,  0.00    , 'qFME',  'traj04',       
     d   19,  3, 1, 3, 1,    5,  0.00    , 'qFME',  'traj05',       
     e   20,  3, 1, 3, 1,   20,  0.00    , 'bFM ',  'traj06',       
     f   21,  3, 2, 1, 4,   51,  0.00    , 'qFM ',  'traj07',       
     g   22,  3, 2, 1, 4,  101,  0.00    , 'qM  ',  'traj08',       
     h   23,  3, 2, 2, 1,   51,  0.00    , 'qFME',  'traj09',       
     i   24,  3, 2, 2, 1,  101,  0.00    , 'qFM ',  'traj10',       
     j   25,  3, 2, 2, 1,  201,  0.00    , 'bFM ',  'traj11',       
     k   26,  3, 2, 2, 1,  401,  0.00    , 'qM  ',  'traj12',       
     l   27,  3, 2, 3, 1,   51,  0.00    , 'qFM ',  'traj13',       
     m   28,  3, 2, 3, 1,  101,  0.00    , 'qFME',  'traj14',       
     n   29,  3, 2, 3, 1,  201,  0.00    , 'qM  ',  'traj15',       
     o   30,  3, 2, 3, 1,  401,  0.00    , 'qFME',  'traj16'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk =  31, 45) /
     a   31,  3, 3, 1, 4,   51,  90.0    , 'qFM ',  'traj17',       
     b   32,  3, 3, 1, 4,   51,  70.0    , 'qFM ',  'traj18',       
     c   33,  3, 3, 1, 4,  101,  90.0    , 'qFM ',  'traj19',       
     d   34,  3, 3, 1, 4,  101,  70.0    , 'qFM ',  'traj20',       
     e   35,  3, 3, 2, 1,   51,  90.0    , 'qFM ',  'traj21',       
     f   36,  3, 3, 2, 1,   51,  70.0    , 'qFM ',  'traj22',       
     g   37,  3, 3, 2, 1,   51, -70.0    , 'qM  ',  'traj23',       
     h   38,  3, 3, 2, 1,  101,  90.0    , 'qFM ',  'traj24',       
     i   39,  3, 3, 2, 1,  101,  70.0    , 'qFM ',  'traj25',       
     j   40,  3, 3, 2, 1,  201,  90.0    , 'qFM ',  'traj26',       
     k   41,  3, 3, 2, 1, 1430,  90.0    , 'qFM ',  'traj27',       
     l   42,  3, 3, 3, 1,   51,  90.0    , 'qFM ',  'traj28',       
     m   43,  3, 3, 3, 1,   51,  70.0    , 'qM  ',  'traj29',       
     n   44,  3, 3, 3, 1,  101,  90.0    , 'qFM ',  'traj30',       
     o   45,  3, 3, 3, 1,  101,  70.0    , 'qFM ',  'traj31'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk =  46, 60) /
     a   46,  3, 3, 3, 1,  201,  90.0    , 'qFM ',  'traj32',       
     b   47,  3, 3, 3, 1, 1430,  90.0    , 'qFM ',  'traj33',       
     c   48,  3, 4, 1, 4,   51,  70.0    , 'qFM ',  'traj34',       
     d   49,  3, 4, 1, 4,  101,  70.0    , 'qFM ',  'traj35',       
     e   50,  3, 4, 2, 1,   51,  70.0    , 'qFM ',  'traj36',       
     f   51,  3, 4, 2, 1,  101,  70.0    , 'qFM ',  'traj37',       
     g   52,  3, 4, 3, 1,   51,  70.0    , 'qFM ',  'traj38',       
     h   53,  3, 4, 3, 1,  101,  70.0    , 'qFM ',  'traj39',       
     i   54,  4, 1, 1, 4,   10,  0.00    , 'qM  ',  'gdrd01',       
     j   55,  4, 1, 2, 1,   10,  0.00    , 'qM  ',  'gdrd02',       
     k   56,  4, 1, 3, 1,   10,  0.00    , 'qM  ',  'gdrd03',       
     l   57,  4, 2, 1, 4,   50,  0.00    , 'qFM ',  'gdrd04',       
     m   58,  4, 2, 2, 1,   50,  0.00    , 'qFM ',  'gdrd05',       
     n   59,  4, 2, 3, 1,   50,  0.00    , 'qM  ',  'gdrd06',       
     o   60,  4, 3, 2, 1,   50,  0.00    , 'qF  ',  'gdrd07'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk =  61, 75) /
     a   61,  5, 1, 2, 1,   50, 0.100E-02, 'qFM ',  'brgr01',       
     b   62,  6, 1, 2, 1,    5,  0.00    , 'qFME',  'nlqr01',       
     c   63,  6, 2, 2, 1,    5,  0.00    , 'qFME',  'nlqr02',       
     d   64,  6, 2, 2, 1,   20,  0.00    , 'qFME',  'nlqr03',       
     e   65,  6, 2, 2, 1,   50,  0.00    , 'qM  ',  'nlqr04',       
     f   66,  7, 1, 1, 4,   20,  1.00    , 'qFM ',  'clym01',       
     g   67,  7, 1, 1, 4,   50,  1.00    , 'qFM ',  'clym02',       
     h   68,  7, 1, 1, 4,  100,  1.00    , 'qFM ',  'clym03',       
     i   69,  7, 1, 2, 1,   10,  1.00    , 'qFM ',  'clym04',       
     j   70,  7, 1, 2, 1,   50,  1.00    , 'qM  ',  'clym05',       
     k   71,  7, 1, 2, 1,  100,  1.00    , 'qM  ',  'clym06',       
     l   72,  7, 1, 3, 1,   10,  1.00    , 'qFM ',  'clym07',       
     m   73,  7, 1, 3, 1,   50,  1.00    , 'qFM ',  'clym08',       
     n   74,  7, 1, 3, 1,  100,  1.00    , 'qFM ',  'clym09',       
     o   75,  7, 2, 1, 4,   20,  1.00    , 'qFM ',  'clym10'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk =  76, 90) /
     a   76,  7, 2, 1, 4,   50,  1.00    , 'qFM ',  'clym11',       
     b   77,  7, 2, 1, 4,  100,  1.00    , 'qFM ',  'clym12',       
     c   78,  7, 2, 2, 1,   10,  1.00    , 'qFM ',  'clym13',       
     d   79,  7, 2, 2, 1,   50,  1.00    , 'qM  ',  'clym14',       
     e   80,  7, 2, 2, 1,  100,  1.00    , 'qM  ',  'clym15',       
     f   81,  7, 2, 3, 1,   10,  1.00    , 'qFM ',  'clym16',       
     g   82,  7, 2, 3, 1,   50,  1.00    , 'qFM ',  'clym17',       
     h   83,  7, 2, 3, 1,  100,  1.00    , 'qFM ',  'clym18',       
     i   84,  7, 3, 2, 1,   10,  1.00    , 'qFM ',  'clym19',       
     j   85,  7, 4, 2, 1,   10,  1.00    , 'qFM ',  'clym20',       
     k   86,  7, 5, 2, 1,   10,  1.00    , 'qF  ',  'clym21',       
     l   87,  7, 6, 2, 1,   10,  1.00    , 'qF  ',  'clym22',       
     m   88,  7, 7, 2, 1,   10,  1.00    , 'qF  ',  'clym23',       
     n   89,  8, 1, 2, 1,   10,  0.00    , 'qFM ',  'capt01',       
     o   90,  8, 1, 2, 1,   10,  1.00    , 'qFM ',  'capt02'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk =  91,105) /
     a   91,  8, 2, 2, 1,   10,  0.00    , 'qFM ',  'capt03',       
     b   92,  8, 2, 2, 1,   10,  1.00    , 'qFME',  'capt04',       
     c   93,  8, 3, 2, 1,   10,  0.00    , 'qFM ',  'capt05',       
     d   94,  8, 3, 2, 1,   10,  1.00    , 'qFM ',  'capt06',       
     e   95,  8, 4, 2, 1,   10,  0.00    , 'qFM ',  'capt07',       
     f   96,  8, 4, 2, 1,   10,  1.00    , 'qFME',  'capt08',       
     g   97,  8, 5, 2, 1,   10,  0.00    , 'qM  ',  'capt09',       
     h   98,  8, 5, 2, 1,   51,  1.00    , 'qFM ',  'capt10',       
     i   99,  8, 6, 2, 1,   10,  0.00    , 'qFM ',  'capt11',       
     j  100,  8, 6, 2, 1,   51,  1.00    , 'qM  ',  'capt12',       
     k  101,  8, 7, 2, 1,   10,  1.00    , 'qM  ',  'capt13',       
     l  102,  8, 8, 2, 1,   10,  1.00    , 'qFME',  'capt14',       
     m  103,  9, 1, 2, 1,   20,  0.00    , 'qM  ',  'gydn01',       
     n  104,  9, 2, 2, 1,   20,  0.00    , 'qM  ',  'gydn02',       
     o  105, 10, 1, 2, 1,   10,  1.00    , 'qM  ',  'brac01'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk = 106,120) /
     a  106, 11, 1, 2, 1,   10,  0.00    , 'bFM ',  'wind01',       
     b  107, 12, 1, 2, 1,   10,  0.00    , 'qFM ',  'orbt01',       
     c  108, 12, 2, 2, 1,   10,  0.00    , 'qFM ',  'orbt02',       
     d  109, 12, 5, 2, 1,   10, -1.00    , 'qFME',  'orbt03',       
     e  110, 13,-1, 2, 1,   50,  2.00    , 'qF  ',  'orbe01',       
     f  111, 13, 1, 2, 1,   50,  1.00    , 'qM  ',  'orbe02',       
     g  112, 13, 1, 2, 1,   50,  2.00    , 'qM  ',  'orbe03',       
     h  113, 13, 1, 2, 1,   50, -1.00    , 'qFM ',  'orbe04',       
     i  114, 13, 2, 2, 1,   50,  1.00    , 'qM  ',  'orbe05',       
     j  115, 13, 2, 2, 1,   50,  2.00    , 'qM  ',  'orbe06',       
     k  116, 13, 2, 2, 1,   50, -1.00    , 'qFM ',  'orbe07',       
     l  117, 13, 3, 2, 1,  100,  1.00    , 'qFM ',  'orbe08',       
     m  118, 13, 3, 2, 1,  100,  2.00    , 'qFM ',  'orbe09',       
     n  119, 13, 3, 2, 1,  100, -1.00    , 'qFM ',  'orbe10',       
     o  120, 13, 4, 2, 1,  150,  1.00    , 'qFM ',  'orbe11'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk = 121,135) /
     a  121, 13, 4, 2, 1,  150,  2.00    , 'qFM ',  'orbe12',       
     b  122, 13, 4, 2, 1,  150, -1.00    , 'qFM ',  'orbe13',       
     c  123, 14, 1, 2, 1,   25,  0.00    , 'bFM ',  'hang01',       
     d  124, 14, 2, 2, 1,   25,  0.00    , 'bFM ',  'hang02',       
     e  125, 15, 1, 2, 1,   25,  800.    , 'qFM ',  'aotv01',       
     f  126, 16, 2, 2, 1,    1,  0.00    , 'qFM ',  'plnt01',       
     g  127, 17, 1, 2, 1,   10, 0.100E-03, 'qFM ',  'ashr01',       
     h  128, 18, 1, 2, 1,   10,  0.00    , 'qFM ',  'kplr01',       
     i  129, 10, 2, 2, 1,   10, 0.100    , 'qFM ',  'brac02',       
     j  130, 19, 1, 2, 1,   10,  0.00    , 'qFM ',  'lowt01',       
     k  131, 20, 1, 2, 1,   10,  0.00    , 'qFM ',  'chmr01',       
     l  132, 20, 2, 2, 1,   10,  0.00    , 'qFM ',  'chmr02',       
     m  133, 20, 3, 2, 1,   10,  0.00    , 'qFM ',  'chmr03',       
     n  134, 20, 4, 2, 1,   10,  0.00    , 'qFM ',  'chmr04',       
     o  135, 20, 5, 2, 1,   10,  0.00    , 'qFM ',  'chmr05'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk = 136,150) /
     a  136, 20, 6, 2, 1,   10,  0.00    , 'qFM ',  'chmr06',       
     b  137, 20, 7, 2, 1,   20,  0.00    , 'qFM ',  'chmr07',       
     c  138, 20, 8, 2, 1,   40,  0.00    , 'qFM ',  'chmr08',       
     d  139, 20, 9, 2, 1,   50,  0.00    , 'qFM ',  'chmr09',       
     e  140, 20,10, 2, 1,   10,  0.00    , 'qFM ',  'chmr10',       
     f  141, 21, 2, 2, 1,   10,  0.00    , 'qFM ',  'putt01',       
     g  142, 22, 1, 2, 1,   10,  0.00    , 'qFM ',  'bang01',       
     h  143, 17, 2, 2, 1,   10, 0.100E-03, 'qFM ',  'ashr02',       
     i  144, 17, 3, 2, 1,   10, 0.100E-03, 'qFM ',  'ashr03',       
     j  145, 17, 4, 2, 1,   10, 0.100E-03, 'qFM ',  'ashr04',       
     k  146,  1, 3, 2, 1,   10, 0.100E-03, 'qFM ',  'qlin03',       
     l  147,  1, 4, 2, 1,   10,  0.00    , 'qFM ',  'qlin04',       
     m  148, 23, 1, 2, 1,   10,  10.0    , 'qFM ',  'heat01',       
     n  149, 24, 1, 2, 1,   10,  0.00    , 'qFM ',  'gsoc01',       
     o  150, 23, 2, 2, 1,   50,  50.0    , 'qFM ',  'heat02'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk = 151,165) /
     a  151, 25, 1, 2, 1,   25,  0.00    , 'qFM ',  'arao01',       
     b  152, 25, 2, 2, 1,   26, -1.00    , 'bFM ',  'arao02',       
     c  153, 17, 5, 2, 1,   40, 0.100E-05, 'qFM ',  'ashr05',       
     d  154, 26, 1, 0, 1,    2,  2.00    , 'bFM ',  'imp201',       
     e  155, 26,-1, 0, 1,    2,  2.00    , 'bFM ',  'imp202',       
     f  156, 26, 2, 0, 1,    2,  2.00    , 'qFM ',  'imp203',       
     g  157, 26,-2, 0, 1,    2,  2.00    , 'qFM ',  'imp204',       
     h  158, 27, 1, 2, 1,   10,  2.00    , 'qFM ',  'jmp201',       
     i  159, 27, 2, 2, 1,   10,  2.00    , 'qFM ',  'jmp202',       
     j  160, 28, 1, 2, 1,   10,  1.00    , 'qFM ',  'brn201',       
     k  161, 28, 1, 2, 1,   10,  2.00    , 'qFM ',  'brn202',       
     l  162, 28, 2, 2, 1,   10,  1.00    , 'qFM ',  'brn203',       
     m  163, 28, 2, 2, 1,   10,  2.00    , 'bFM ',  'brn204',       
     n  164, 29, 1, 3, 1,   21,  0.00    , 'qFM ',  'alpr01',       
     o  165, 30, 1, 2, 1,   25,  0.00    , 'qFM ',  'mirv01'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk = 166,180) /
     a  166, 31, 1, 2, 1,   10,  0.00    , 'qFM ',  'robo01',       
     b  167, 31, 2, 2, 1,   10,  0.00    , 'qFM ',  'robo02',       
     c  168, 31, 3, 2, 1,   10,  0.00    , 'qFM ',  'robo03',       
     d  169, 31, 4, 2, 1,    5,  0.00    , 'qFM ',  'robo04',       
     e  170, 32, 1, 2, 2,   20,  0.00    , 'qFM ',  'skwz01',       
     f  171, 32, 2, 2, 2,   20,  0.00    , 'bFM ',  'skwz02',       
     g  172, 32, 3, 2, 2,   20,  0.00    , 'bFM ',  'skwz03',       
     h  173, 32, 4, 2, 2,   20,  0.00    , 'bFM ',  'skwz04',       
     i  174, 33, 1, 2, 1,   10,  0.00    , 'qFM ',  'dlay01',       
     j  175, 34, 1, 2, 1,    2,  0.00    , 'qFM ',  'nzym01',       
     k  176, 38, 2, 2, 1,   50,  0.00    , 'qFM ',  'pndl02',       
     l  177, 36, 1, 2, 1,   10,  0.00    , 'qFM ',  'lwbr01',       
     m  178, 37, 1, 2, 1,   10,  0.00    , 'qFM ',  'aqua01',       
     n  179, 38, 1, 2, 1,   10,  0.00    , 'qFM ',  'pndl01',       
     o  180, 39, 1, 2, 1,   10,  0.00    , 'bFM ',  'tran01'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk = 181,195) /
     a  181, 40, 1, 2, 1,   10,  0.00    , 'qFM ',  'zrml01',       
     b  182, 41, 1, 2, 1,   10,  20.0    , 'qFM ',  'hdae01',       
     c  183, 42, 1, 2, 1,   10, 0.100E-01, 'qFM ',  'cran01',       
     d  184, 43, 1, 2, 1,   10,  0.00    , 'qFM ',  'vpol01',       
     e  185, 43, 2, 2, 1,   10,  0.00    , 'qFM ',  'vpol02',       
     f  186, 43, 3, 2, 1,   10,  0.00    , 'qFM ',  'vpol03',       
     g  187, 43, 4, 2, 1,   10,  0.00    , 'qFM ',  'vpol04',       
     h  188, 43, 5, 2, 1,   10,  0.00    , 'qFM ',  'vpol05',       
     i  189, 43, 6, 2, 1,   10,  0.00    , 'qFM ',  'vpol06',       
     j  190, 44, 1, 2, 1,   10,  0.00    , 'bFM ',  'ffrb01',       
     k  191, 45, 1, 2, 1,   50,  0.00    , 'bFM ',  'jshi01',       
     l  192, 45, 2, 2, 1,   50,  0.00    , 'qFM ',  'jshi02',       
     m  193, 46, 1, 2, 1,   50,  0.00    , 'qFM ',  'lnht01',       
     n  194, 46, 2, 2, 1,   50,  0.00    , 'qFM ',  'lnht02',       
     o  195, 47, 1, 2, 2,   25,  0.00    , 'bFM ',  'chan01'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk = 196,210) /
     a  196, 30, 2, 5, 1,   25,  0.00    , 'qFM ',  'mirv02',       
     b  197, 30, 3, 2, 1,   25,  0.00    , 'qFM ',  'mirv03',       
     c  198, 48, 3, 2, 1,   20,  0.00    , 'qFM ',  'ssmd01',       
     d  199, 49, 1, 2, 1,   10,  0.00    , 'qFM ',  'ltsp01',       
     e  200, 49, 2, 2, 1,   10,  0.00    , 'qFM ',  'ltsp02',       
     f  201, 28,-1, 2, 1,   50,  1.00    , 'bFM ',  'brn205',       
     g  202, 28,-2, 2, 1,   50,  2.00    , 'qFM ',  'brn206',       
     h  203, 50, 1, 2, 1,   10,  0.00    , 'qFM ',  'dlt301',       
     i  204, 51, 2, 3, 1,   10,  0.00    , 'qFM ',  'rcsp01',       
     j  205, 51, 4, 3, 1,   10,  0.00    , 'qFM ',  'rcsp02',       
     k  206, 51, 6, 3, 1,   10,  0.00    , 'qFM ',  'rcsp03',       
     l  207, 51, 8, 3, 1,   10,  0.00    , 'qFM ',  'rcsp04',       
     m  208,  5, 3, 2, 1,   10, 0.100E-02, 'qFM ',  'brgr02',       
     n  209, 52, 1, 2, 1,   10,  0.00    , 'qFM ',  'aomp01',       
     o  210, 52, 2, 2, 1,   10,  0.00    , 'qFM ',  'aomp02'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk = 211,225) /
     a  211, 53, 1, 2, 1,   10,  0.00    , 'bFM ',  'tb2s01',       
     b  212, 54, 1, 3, 1,  100,  0.00    , 'qFM ',  'fhoc01',       
     c  213, 55, 1, 2, 1,    3,  0.00    , 'qFM ',  'mrck01',       
     d  214, 56, 1, 2, 1,   10,  0.00    , 'bFM ',  'asyq01',       
     e  215, 56, 2, 2, 1,   10,  0.00    , 'qFM ',  'asyq02',       
     f  216, 57, 1, 2, 1,   10,  0.00    , 'qFM ',  'rayl01',       
     g  217, 57, 2, 2, 1,   10,  0.00    , 'qFM ',  'rayl02',       
     h  218, 57, 3, 2, 1,   10,  20.0    , 'qFM ',  'rayl03',       
     i  219, 57, 4, 2, 1,   10, 0.100E-01, 'qFM ',  'rayl04',       
     j  220, 57, 5, 2, 1,   10,  0.00    , 'qFM ',  'rayl05',       
     k  221, 58,-1, 2, 1,   10, 0.100    , 'qFM ',  'medi01',       
     l  222, 58, 1, 2, 1,   10, 0.100    , 'qFM ',  'medi02',       
     m  223, 58,-2, 2, 1,   10, 0.200    , 'qFM ',  'medi03',       
     n  224, 58, 2, 2, 1,   10, 0.200    , 'qFM ',  'medi04',       
     o  225, 58,-3, 2, 1,   10, 0.500    , 'qFM ',  'medi05'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk = 226,240) /
     a  226, 58, 3, 2, 1,   10, 0.500    , 'qFM ',  'medi06',       
     b  227, 59, 1, 2, 1,   10,  0.00    , 'bFM ',  'soar01',       
     c  228, 60, 1, 2, 1,   10,  0.00    , 'qFM ',  'pnav01',       
     d  229, 60, 2, 2, 1,   10,  0.00    , 'qFM ',  'pnav02',       
     e  230, 61, 3, 2, 1,   10, 0.500E-01, 'qFM ',  'foic01',       
     f  231, 62, 1, 2, 2,   10,  0.00    , 'qFM ',  'ccgo01',       
     g  232, 62, 2, 2, 2,   10,  0.00    , 'qFM ',  'ccgo02',       
     h  233, 62,-2, 3, 1,   10,  0.00    , 'bFM ',  'ccgo03',       
     i  234, 62, 3, 2, 2,   10,  0.00    , 'qFM ',  'ccgo04',       
     j  235, 62, 4, 2, 2,   10,  0.00    , 'bFM ',  'ccgo05',       
     k  236, 63, 1, 3, 1,   10, 0.100    , 'qFM ',  'fonr01',       
     l  237, 64, 1, 2, 2,   10,  1.00    , 'bFM ',  'apin01',       
     m  238, 65, 1, 3, 1,   50,  1.00    , 'qFM ',  'lrnz01',       
     n  239, 65, 2, 3, 1,   50,  1.00    , 'qM  ',  'lrnz02',       
     o  240, 66, 1, 3, 1,   50,  0.00    , 'qM  ',  'rslr01'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk = 241,255) /
     a  241, 66, 2, 3, 1,   50,  0.00    , 'bFM ',  'rslr02',       
     b  242, 67, 1, 3, 1,   10, 0.100E-01, 'qFM ',  'rntr01',       
     c  243, 68, 1, 2, 1,   10, 0.500E-01, 'qFM ',  'blrs01',       
     d  244, 68, 2, 2, 1,   10, 0.500E-01, 'qFM ',  'blrs02',       
     e  245, 68, 3, 2, 1,   10, 0.500E-01, 'qFM ',  'blrs03',       
     f  246, 68, 4, 2, 1,   10, 0.500E-01, 'qFM ',  'blrs04',       
     g  247, 68, 5, 2, 1,   10, 0.500E-01, 'qFM ',  'blrs05',       
     h  248, 68, 6, 2, 1,   10, 0.500E-01, 'qFM ',  'blrs06',       
     i  249, 69, 1, 2, 2,   12, 0.100E-01, 'qFM ',  'zwag01',       
     j  250, 70, 1, 2, 1,  100,  0.00    , 'qFM ',  'hiv101',       
     k  251, 71, 3, 2, 1,   10,  0.00    , 'qFM ',  'cgro01',       
     l  252, 72, 1, 2, 1,   10,  0.00    , 'bFM ',  'stab01',       
     m  253, 72, 2, 2, 1,   10,  0.00    , 'bFM ',  'stab02',       
     n  254, 72, 3, 2, 1,   10,  0.00    , 'qFM ',  'stab03',       
     o  255, 73, 1, 2, 1,   10,  0.00    , 'qFM ',  'stwt01'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk = 256,270) /
     a  256, 74, 1, 2, 1,   10,  0.00    , 'qFM ',  'stsp01',       
     b  257, 74, 2, 2, 1,   10,  0.00    , 'qFM ',  'stsp02',       
     c  258, 75, 1, 2, 1,   10,  0.00    , 'qFM ',  'knbr01',       
     d  259, 76, 1, 2, 1,   10,  0.00    , 'qFM ',  'knba01',       
     e  260,  4, 4, 2, 1,   50,  0.00    , 'qFM ',  'gdrd08',       
     f  261,  4, 4, 2, 1,   50, -1.00    , 'qFM ',  'gdrd09',       
     g  262, 77, 1, 2, 1,  150,  2.00    , 'qFM ',  'lthr01',       
     h  263, 15, 2, 2, 1,   25,  600.    , 'qFM ',  'aotv02',       
     i  264, 78, 1, 2, 1,   20,  0.00    , 'qFM ',  'lbrp01',       
     j  265, 78, 2, 2, 1,   20,  0.00    , 'qFM ',  'lbrp02',       
     k  266, 78,-1, 2, 1,   20,  0.00    , 'qFM ',  'lbrp03',       
     l  267, 78,-2, 2, 1,   20,  0.00    , 'qFM ',  'lbrp04',       
     m  268, 79, 1, 2, 1,   20,  0.00    , 'qFM ',  'dock01',       
     n  269, 79, 2, 2, 1,   20,  0.00    , 'qFM ',  'dock02',       
     o  270, 79, 3, 2, 1,   20,  0.00    , 'qFM ',  'dock03'/       
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk = 271,285) /
     a  271,  2, 4, 5, 2,    2,  0.00    , 'qFM ',  'lnts13',       
     b  272,  4, 5, 5, 2,    2,  0.00    , 'qFM ',  'gdrd10',       
     c  273, 14, 3, 2, 1,  100,  0.00    , 'qFM ',  'hang03',       
     d  274, 17, 6, 2, 1,   20,  0.00    , 'bFM ',  'ashr06',       
     e  275, 21, 2, 5, 2,    2,  0.00    , 'qFM ',  'putt02',       
     f  276, 35, 1, 2, 1,   10,  0.00    , 'qFM ',  'rivr01',       
     g  277, 35, 2, 2, 1,   10,  0.00    , 'qFM ',  'rivr02',       
     h  278, 43, 7, 5, 2,    2,  0.00    , 'qFM ',  'vpol07',       
     i  279, 47, 2, 2, 2,   25,  0.00    , 'bFM ',  'chan02',       
     j  280, 47, 3, 2, 2,   25,  0.00    , 'bFM ',  'chan03',       
     k  281, 80, 1, 2, 1,    2,  0.00    , 'qFM ',  'goll01',
     l  282, 80, 2, 2, 1,    2,  0.00    , 'qFM ',  'goll02',
     m  283, 80, 3, 2, 1,    2,  0.00    , 'qFM ',  'goll03',
     n  284, 81, 1, 2, 1,    2,  0.00    , 'qFM ',  'cstr01',
     o  285, 81, 3, 2, 1,    5,  0.00    , 'qFM ',  'cstr02'/
C
      data (itprob(kk),itpode(kk),itpcls(kk),itpmet(kk),itpstg(kk),
     $  itpgrd(kk),tpseed(kk),itpalg(kk),prname(kk), kk = 286,300) /
     a  286, 81, 4, 2, 1,    2,  0.00    , 'qFM ',  'cstr03',
     b  287, 82, 4, 2, 1,    2,  0.00    , 'qFM ',  'fish01',
     c  288, 83, 4, 2, 1,    2,  0.00    , 'qFM ',  'cst201',
     d  289, 84, 1, 2, 1,    3,  0.00    , 'qFM ',  'stgl01',
     e  290, 85, 1, 2, 1,    2,  0.00    , 'qFM ',  'mncx01',
     f  291, 85, 2, 2, 1,    2,  0.00    , 'qFM ',  'mncx02',
     g  292, 85, 4, 2, 1,    2,  0.00    , 'qFM ',  'mncx03',
     h  293, 86, 1, 2, 1,    5,  0.00    , 'qFM ',  'pdly01',
     i  294, 87, 1, 2, 1,   10,  0.00    , 'qFM ',  'rbrm01',
     j  295, 88, 1, 2, 1,   10,  0.00    , 'qFM ',  'tumr01',
     k  296, 88, 2, 2, 1,   10,  0.00    , 'qFM ',  'tumr02',
     l  297, 88, 3, 2, 1,   10,  0.00    , 'qFM ',  'tumr03',
     m  298, 89, 3, 2, 1,   20,  0.00    , 'qFM ',  'lbri01',
     n  299, 89, 4, 2, 1,   20,  0.00    , 'qFM ',  'lbri02',
     o  300,  0, 0, 0, 0,    0,  0.00    , 'FM  ',  'zzzbad'/
c
C
      END
