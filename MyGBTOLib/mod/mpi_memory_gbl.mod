  S  Ð   k820309    ,          2021.5.0    fïc                                                                                                          
       /home/nico/Workspace/Progs/Ints_From_gbtolib/MyGBTOLib/source/mpi_memory_mod.F90 MPI_MEMORY_GBL                      @                              
       WP EP WP_BYTES EP1_BYTES LONGINT                      @                              
       MPIINT MPIADDR LOCAL_RANK SHARED_ENABLED SHARED_COMMUNICATOR MPI_XERMSG MPI_MOD_RANK                      @                              
       LEVEL3          @ 0        @                              
       C_PTR C_F_POINTER C_NULL_PTR            0                                        
       INT32 INT64                                                        u #MPI_MEMORY_ALLOCATE_INTEGER_ALC    #MPI_MEMORY_ALLOCATE_INTEGER_PTR                                                           u #MPI_MEMORY_DEALLOCATE_INTEGER_ALC    #MPI_MEMORY_DEALLOCATE_INTEGER_PTR 	                                                          u #MPI_MEMORY_ALLOCATE_INT32_2DIM_ALC 
   #MPI_MEMORY_ALLOCATE_INT32_2DIM_PTR    #MPI_MEMORY_ALLOCATE_INT64_2DIM_ALC    #MPI_MEMORY_ALLOCATE_INT64_2DIM_PTR                                                           u #MPI_MEMORY_DEALLOCATE_INT32_2DIM_ALC    #MPI_MEMORY_DEALLOCATE_INT32_2DIM_PTR    #MPI_MEMORY_DEALLOCATE_INT64_2DIM_ALC    #MPI_MEMORY_DEALLOCATE_INT64_2DIM_PTR                                                           u #MPI_MEMORY_ALLOCATE_REAL_WP_ALC    #MPI_MEMORY_ALLOCATE_REAL_EP_ALC    #MPI_MEMORY_ALLOCATE_REAL_WP_PTR    #MPI_MEMORY_ALLOCATE_REAL_EP_PTR                                                           u #MPI_MEMORY_DEALLOCATE_REAL_WP_ALC    #MPI_MEMORY_DEALLOCATE_REAL_EP_ALC    #MPI_MEMORY_DEALLOCATE_REAL_WP_PTR    #MPI_MEMORY_DEALLOCATE_REAL_EP_PTR                                                           u #MPI_MEMORY_ALLOCATE_REAL_2DIM_WP_ALC    #MPI_MEMORY_ALLOCATE_REAL_2DIM_EP_ALC    #MPI_MEMORY_ALLOCATE_REAL_2DIM_WP_PTR    #MPI_MEMORY_ALLOCATE_REAL_2DIM_EP_PTR                                                           u #MPI_MEMORY_DEALLOCATE_REAL_2DIM_WP_ALC    #MPI_MEMORY_DEALLOCATE_REAL_2DIM_EP_ALC    #MPI_MEMORY_DEALLOCATE_REAL_2DIM_WP_PTR     #MPI_MEMORY_DEALLOCATE_REAL_2DIM_EP_PTR !                                                        u #C_F_POINTER_SCALAR_NEW "   #C_F_POINTER_ARRAY1 #   #C_F_POINTER_ARRAY2 $   #C_F_POINTER_ARRAY4 %   #C_F_POINTER_ARRAY8 &                     @                          '     '                    #PTR (                 D                             (                                    @                        )     '                    #PTR *                 D                            *                                                                +                                                                                                      ,                                                                                                      -                                                      8                                             .                                                      16                                             /                                                                                                      0                                                                                                      1                                                                                                   2                                                       3                                                      4            #         @                                  5                    #MOD_NAME 6   #ROUTINE_NAME 7   #ERR_MSG 8   #ERR 9   #LEVEL :             
                                6                    1           
                                7                    1           
                                8                    1           
                                  9                     
                                  :           #         @                                  ;                    #RANK <   #COMM =                                             <                      
                                =                                                       >                                                        ?                                                          X#C_PTR '      n                                      0  #C_PTR '                                                                      @                                                      4                                            A                                                      8#         @     @          @                "     	               #CPTR B   #FPTR C   n                             Å              Cc_f_pointer_set_scalar                          
                                 B                   #C_PTR '                    @                       C            #         @     @          @                #     	               #CPTR D   #FPTR E   #SHAPE F   n                          É              Cc_f_pointer_set_desc1                             
                                 D                   #C_PTR '                   @                      E                                  &                                                     
                                F                                 &                                           #         @     @          @                $     	               #CPTR G   #FPTR H   #SHAPE I   n                          Ì              Cc_f_pointer_set_desc2                             
                                 G                   #C_PTR '                   @                      H                                  &                                                     
                                I                                 &                                           #         @     @         @                %     	               #CPTR J   #FPTR K   #SHAPE L   n                          Ï              Cc_f_pointer_set_desc4                             
                                 J                   #C_PTR '                   @                      K                                  &                                                     
                                L                                 &                                           #         @     @          @                &     	               #CPTR M   #FPTR N   #SHAPE O   n                          Ò              Cc_f_pointer_set_desc8                             
                                 M                   #C_PTR '                   @                      N                    	              &                                                     
                                O                    
             &                                                                                       P                                                       0%         @    X                                                       #ARRAY Q   #NELEM R   #COMM S           
D                                 Q                                  &                                                     
                                  R                     
                                S           %         @    X                                                       #ARRAY T   #NELEM U   #COMM V            
D P                               T                                  &                                                     
                                  U                     
 @                              V           #         @      X                                                 #ARRAY W   #NELEM X   #WINDOW Y   #COMM Z           
D                                 W                                  &                                                     
                                  X                     
                                  Y                     
                                Z           #         @      X                             	                    #ARRAY [   #NELEM \   #WINDOW ]   #COMM ^            
D                                [                                  &                                                     
                                  \                     
  @                               ]                     
 @                              ^           %         @    X                            
                           #ARRAY _   #NELEM1 `   #NELEM2 a   #COMM b           
D                                _                                  &                   &                                                     
                                  `                     
                                  a                     
                                b           %         @    X                                                       #ARRAY c   #NELEM1 d   #NELEM2 e   #COMM f            
D P                              c                                  &                   &                                                     
                                  d                     
                                  e                     
 @                              f           %         @    X                                                       #ARRAY g   #NELEM1 h   #NELEM2 i   #COMM j           
D                                g                                  &                   &                                                     
                                  h                     
                                  i                     
                                j           %         @    X                                                       #ARRAY k   #NELEM1 l   #NELEM2 m   #COMM n            
D P                              k                                  &                   &                                                     
                                  l                     
                                  m                     
 @                              n           #         @      X                                                 #ARRAY o   #NELEM p   #WINDOW q   #COMM r           
D                                o                    
              &                   &                                                     
                                  p                     
                                  q                     
                                r           #         @      X                                                 #ARRAY s   #NELEM t   #WINDOW u   #COMM v            
D                               s                                  &                   &                                                     
                                  t                     
  @                               u                     
 @                              v           #         @      X                                                 #ARRAY w   #NELEM x   #WINDOW y   #COMM z           
D                                w                                  &                   &                                                     
                                  x                     
                                  y                     
                                z           #         @      X                                                 #ARRAY {   #NELEM |   #WINDOW }   #COMM ~            
D                               {                                  &                   &                                                     
                                  |                     
  @                               }                     
 @                              ~           %         @    X                                                       #ARRAY    #NELEM    #COMM            
D                                                   
               &                                                     
                                                       
                                           %         @    X                                                       #ARRAY    #NELEM    #COMM            
D                                                                  &                                                     
                                                       
                                           %         @    X                                                       #ARRAY    #NELEM    #COMM             
D P                                                 
               &                                                     
                                                       
 @                                         %         @    X                                                       #ARRAY    #NELEM    #COMM             
D P                                                                &                                                     
                                                       
 @                                         #         @      X                                                 #ARRAY    #NELEM    #WINDOW    #COMM            
D                                                   
               &                                                     
                                                       
                                                       
                                           #         @      X                                                 #ARRAY    #NELEM    #WINDOW    #COMM            
D                                                    	              &                                                     
                                                       
                                                       
                                           #         @      X                                                 #ARRAY    #NELEM    #WINDOW    #COMM             
D                                                  
               &                                                     
                                                       
  @                                                    
 @                                         #         @      X                                                 #ARRAY    #NELEM    #WINDOW    #COMM             
D                                                                 &                                                     
                                                       
  @                                                    
                                           %         @    X                                                       #ARRAY    #NELEM1    #NELEM2    #COMM            
D                                                   
               &                   &                                                     
                                                       
                                                       
                                           %         @    X                                                       #ARRAY    #NELEM1     #NELEM2 ¡   #COMM ¢           
D                                                                  &                   &                                                     
                                                        
                                  ¡                     
                                ¢           %         @    X                                                       #ARRAY £   #NELEM1 ¤   #NELEM2 ¥   #COMM ¦            
D P                              £                   
               &                   &                                                     
                                  ¤                     
                                  ¥                     
 @                              ¦           %         @    X                                                       #ARRAY §   #NELEM1 ¨   #NELEM2 ©   #COMM ª            
D P                              §                                  &                   &                                                     
                                  ¨                     
                                  ©                     
 @                              ª           #         @      X                                                 #ARRAY «   #NELEM ¬   #WINDOW ­   #COMM ®           
D                                «                   
               &                   &                                                     
                                  ¬                     
                                  ­                     
                                ®           #         @      X                                                 #ARRAY ¯   #NELEM °   #WINDOW ±   #COMM ²           
D                                ¯                                  &                   &                                                     
                                  °                     
                                  ±                     
                                ²           #         @      X                                                  #ARRAY ³   #NELEM ´   #WINDOW µ   #COMM ¶            
D                               ³                   
               &                   &                                                     
                                  ´                     
  @                               µ                     
 @                              ¶           #         @      X                             !                    #ARRAY ·   #NELEM ¸   #WINDOW ¹   #COMM º            
D                               ·                                  &                   &                                                     
                                  ¸                     
  @                               ¹                     
 @                              º           #         @                                   »                     #         @                                  ¼                    #ALLOC_SIZE ½   #GROUPCOMM ¾   #BASEPTR ¿   #WIN À             
                                 ½                     
                                 ¾                     D                                ¿                    #C_PTR '             D                                À            #         @                                  Á                    #WINDOW Â   #COMM Ã             
                                  Â                     
                                Ã           #         @                                  Ä                    #WINDOW Å   #IERROR Æ             
                                  Å                     
                                Æ                   h      fn#fn      a   J  PRECISN_GBL    i     J  MPI_GBL    þ  G   J  CONST_GBL    E  ]   J  ISO_C_BINDING     ¢  L   J  ISO_FORTRAN_ENV 0   î         gen@MPI_MEMORY_ALLOCATE_INTEGER 2   x         gen@MPI_MEMORY_DEALLOCATE_INTEGER 5     à       gen@MPI_MEMORY_ALLOCATE_INTEGER_2DIM 7   æ  è       gen@MPI_MEMORY_DEALLOCATE_INTEGER_2DIM -   Î  Ô       gen@MPI_MEMORY_ALLOCATE_REAL /   ¢  Ü       gen@MPI_MEMORY_DEALLOCATE_REAL 2   ~  è       gen@MPI_MEMORY_ALLOCATE_REAL_2DIM 4   f  ð       gen@MPI_MEMORY_DEALLOCATE_REAL_2DIM .   V	  ¼       gen@C_F_POINTER+ISO_C_BINDING $   
  Y       C_PTR+ISO_C_BINDING ,   k
  H   %   C_PTR%PTR+ISO_C_BINDING=PTR $   ³
  Y       C_PTR+ISO_C_BINDING ,     H   %   C_PTR%PTR+ISO_C_BINDING=PTR    T  p       WP+PRECISN_GBL    Ä  p       EP+PRECISN_GBL %   4  q       WP_BYTES+PRECISN_GBL &   ¥  r       EP1_BYTES+PRECISN_GBL $     p       LONGINT+PRECISN_GBL      p       MPIINT+MPI_GBL     ÷  p       MPIADDR+MPI_GBL #   g  @       LOCAL_RANK+MPI_GBL '   §  @       SHARED_ENABLED+MPI_GBL ,   ç  @       SHARED_COMMUNICATOR+MPI_GBL #   '         MPI_XERMSG+MPI_GBL ,   °  L   a   MPI_XERMSG%MOD_NAME+MPI_GBL 0   ü  L   a   MPI_XERMSG%ROUTINE_NAME+MPI_GBL +   H  L   a   MPI_XERMSG%ERR_MSG+MPI_GBL '     @   a   MPI_XERMSG%ERR+MPI_GBL )   Ô  @   a   MPI_XERMSG%LEVEL+MPI_GBL %     \       MPI_MOD_RANK+MPI_GBL *   p  @   a   MPI_MOD_RANK%RANK+MPI_GBL *   °  @   a   MPI_MOD_RANK%COMM+MPI_GBL !   ð  @       LEVEL3+CONST_GBL )   0  ×       C_NULL_PTR+ISO_C_BINDING &     q       INT32+ISO_FORTRAN_ENV &   x  q       INT64+ISO_FORTRAN_ENV 5   é  ·      C_F_POINTER_SCALAR_NEW+ISO_C_BINDING :      S   a   C_F_POINTER_SCALAR_NEW%CPTR+ISO_C_BINDING :   ó  @   a   C_F_POINTER_SCALAR_NEW%FPTR+ISO_C_BINDING 1   3  Á      C_F_POINTER_ARRAY1+ISO_C_BINDING 6   ô  S   a   C_F_POINTER_ARRAY1%CPTR+ISO_C_BINDING 6   G     a   C_F_POINTER_ARRAY1%FPTR+ISO_C_BINDING 7   Ó     a   C_F_POINTER_ARRAY1%SHAPE+ISO_C_BINDING 1   _  Á      C_F_POINTER_ARRAY2+ISO_C_BINDING 6      S   a   C_F_POINTER_ARRAY2%CPTR+ISO_C_BINDING 6   s     a   C_F_POINTER_ARRAY2%FPTR+ISO_C_BINDING 7   ÿ     a   C_F_POINTER_ARRAY2%SHAPE+ISO_C_BINDING 1     Á      C_F_POINTER_ARRAY4+ISO_C_BINDING 6   L  S   a   C_F_POINTER_ARRAY4%CPTR+ISO_C_BINDING 6        a   C_F_POINTER_ARRAY4%FPTR+ISO_C_BINDING 7   +     a   C_F_POINTER_ARRAY4%SHAPE+ISO_C_BINDING 1   ·  Á      C_F_POINTER_ARRAY8+ISO_C_BINDING 6   x  S   a   C_F_POINTER_ARRAY8%CPTR+ISO_C_BINDING 6   Ë     a   C_F_POINTER_ARRAY8%FPTR+ISO_C_BINDING 7   W     a   C_F_POINTER_ARRAY8%SHAPE+ISO_C_BINDING    ã  q       LOCAL_MASTER 0   T  p       MPI_MEMORY_ALLOCATE_INTEGER_ALC 6   Ä     a   MPI_MEMORY_ALLOCATE_INTEGER_ALC%ARRAY 6   P  @   a   MPI_MEMORY_ALLOCATE_INTEGER_ALC%NELEM 5     @   a   MPI_MEMORY_ALLOCATE_INTEGER_ALC%COMM 0   Ð  p       MPI_MEMORY_ALLOCATE_INTEGER_PTR 6   @      a   MPI_MEMORY_ALLOCATE_INTEGER_PTR%ARRAY 6   Ì   @   a   MPI_MEMORY_ALLOCATE_INTEGER_PTR%NELEM 5   !  @   a   MPI_MEMORY_ALLOCATE_INTEGER_PTR%COMM 2   L!  t       MPI_MEMORY_DEALLOCATE_INTEGER_ALC 8   À!     a   MPI_MEMORY_DEALLOCATE_INTEGER_ALC%ARRAY 8   L"  @   a   MPI_MEMORY_DEALLOCATE_INTEGER_ALC%NELEM 9   "  @   a   MPI_MEMORY_DEALLOCATE_INTEGER_ALC%WINDOW 7   Ì"  @   a   MPI_MEMORY_DEALLOCATE_INTEGER_ALC%COMM 2   #  t       MPI_MEMORY_DEALLOCATE_INTEGER_PTR 8   #     a   MPI_MEMORY_DEALLOCATE_INTEGER_PTR%ARRAY 8   $  @   a   MPI_MEMORY_DEALLOCATE_INTEGER_PTR%NELEM 9   L$  @   a   MPI_MEMORY_DEALLOCATE_INTEGER_PTR%WINDOW 7   $  @   a   MPI_MEMORY_DEALLOCATE_INTEGER_PTR%COMM 3   Ì$  }       MPI_MEMORY_ALLOCATE_INT32_2DIM_ALC 9   I%  ¤   a   MPI_MEMORY_ALLOCATE_INT32_2DIM_ALC%ARRAY :   í%  @   a   MPI_MEMORY_ALLOCATE_INT32_2DIM_ALC%NELEM1 :   -&  @   a   MPI_MEMORY_ALLOCATE_INT32_2DIM_ALC%NELEM2 8   m&  @   a   MPI_MEMORY_ALLOCATE_INT32_2DIM_ALC%COMM 3   ­&  }       MPI_MEMORY_ALLOCATE_INT32_2DIM_PTR 9   *'  ¤   a   MPI_MEMORY_ALLOCATE_INT32_2DIM_PTR%ARRAY :   Î'  @   a   MPI_MEMORY_ALLOCATE_INT32_2DIM_PTR%NELEM1 :   (  @   a   MPI_MEMORY_ALLOCATE_INT32_2DIM_PTR%NELEM2 8   N(  @   a   MPI_MEMORY_ALLOCATE_INT32_2DIM_PTR%COMM 3   (  }       MPI_MEMORY_ALLOCATE_INT64_2DIM_ALC 9   )  ¤   a   MPI_MEMORY_ALLOCATE_INT64_2DIM_ALC%ARRAY :   ¯)  @   a   MPI_MEMORY_ALLOCATE_INT64_2DIM_ALC%NELEM1 :   ï)  @   a   MPI_MEMORY_ALLOCATE_INT64_2DIM_ALC%NELEM2 8   /*  @   a   MPI_MEMORY_ALLOCATE_INT64_2DIM_ALC%COMM 3   o*  }       MPI_MEMORY_ALLOCATE_INT64_2DIM_PTR 9   ì*  ¤   a   MPI_MEMORY_ALLOCATE_INT64_2DIM_PTR%ARRAY :   +  @   a   MPI_MEMORY_ALLOCATE_INT64_2DIM_PTR%NELEM1 :   Ð+  @   a   MPI_MEMORY_ALLOCATE_INT64_2DIM_PTR%NELEM2 8   ,  @   a   MPI_MEMORY_ALLOCATE_INT64_2DIM_PTR%COMM 5   P,  t       MPI_MEMORY_DEALLOCATE_INT32_2DIM_ALC ;   Ä,  ¤   a   MPI_MEMORY_DEALLOCATE_INT32_2DIM_ALC%ARRAY ;   h-  @   a   MPI_MEMORY_DEALLOCATE_INT32_2DIM_ALC%NELEM <   ¨-  @   a   MPI_MEMORY_DEALLOCATE_INT32_2DIM_ALC%WINDOW :   è-  @   a   MPI_MEMORY_DEALLOCATE_INT32_2DIM_ALC%COMM 5   (.  t       MPI_MEMORY_DEALLOCATE_INT32_2DIM_PTR ;   .  ¤   a   MPI_MEMORY_DEALLOCATE_INT32_2DIM_PTR%ARRAY ;   @/  @   a   MPI_MEMORY_DEALLOCATE_INT32_2DIM_PTR%NELEM <   /  @   a   MPI_MEMORY_DEALLOCATE_INT32_2DIM_PTR%WINDOW :   À/  @   a   MPI_MEMORY_DEALLOCATE_INT32_2DIM_PTR%COMM 5    0  t       MPI_MEMORY_DEALLOCATE_INT64_2DIM_ALC ;   t0  ¤   a   MPI_MEMORY_DEALLOCATE_INT64_2DIM_ALC%ARRAY ;   1  @   a   MPI_MEMORY_DEALLOCATE_INT64_2DIM_ALC%NELEM <   X1  @   a   MPI_MEMORY_DEALLOCATE_INT64_2DIM_ALC%WINDOW :   1  @   a   MPI_MEMORY_DEALLOCATE_INT64_2DIM_ALC%COMM 5   Ø1  t       MPI_MEMORY_DEALLOCATE_INT64_2DIM_PTR ;   L2  ¤   a   MPI_MEMORY_DEALLOCATE_INT64_2DIM_PTR%ARRAY ;   ð2  @   a   MPI_MEMORY_DEALLOCATE_INT64_2DIM_PTR%NELEM <   03  @   a   MPI_MEMORY_DEALLOCATE_INT64_2DIM_PTR%WINDOW :   p3  @   a   MPI_MEMORY_DEALLOCATE_INT64_2DIM_PTR%COMM 0   °3  p       MPI_MEMORY_ALLOCATE_REAL_WP_ALC 6    4     a   MPI_MEMORY_ALLOCATE_REAL_WP_ALC%ARRAY 6   ¬4  @   a   MPI_MEMORY_ALLOCATE_REAL_WP_ALC%NELEM 5   ì4  @   a   MPI_MEMORY_ALLOCATE_REAL_WP_ALC%COMM 0   ,5  p       MPI_MEMORY_ALLOCATE_REAL_EP_ALC 6   5     a   MPI_MEMORY_ALLOCATE_REAL_EP_ALC%ARRAY 6   (6  @   a   MPI_MEMORY_ALLOCATE_REAL_EP_ALC%NELEM 5   h6  @   a   MPI_MEMORY_ALLOCATE_REAL_EP_ALC%COMM 0   ¨6  p       MPI_MEMORY_ALLOCATE_REAL_WP_PTR 6   7     a   MPI_MEMORY_ALLOCATE_REAL_WP_PTR%ARRAY 6   ¤7  @   a   MPI_MEMORY_ALLOCATE_REAL_WP_PTR%NELEM 5   ä7  @   a   MPI_MEMORY_ALLOCATE_REAL_WP_PTR%COMM 0   $8  p       MPI_MEMORY_ALLOCATE_REAL_EP_PTR 6   8     a   MPI_MEMORY_ALLOCATE_REAL_EP_PTR%ARRAY 6    9  @   a   MPI_MEMORY_ALLOCATE_REAL_EP_PTR%NELEM 5   `9  @   a   MPI_MEMORY_ALLOCATE_REAL_EP_PTR%COMM 2    9  t       MPI_MEMORY_DEALLOCATE_REAL_WP_ALC 8   :     a   MPI_MEMORY_DEALLOCATE_REAL_WP_ALC%ARRAY 8    :  @   a   MPI_MEMORY_DEALLOCATE_REAL_WP_ALC%NELEM 9   à:  @   a   MPI_MEMORY_DEALLOCATE_REAL_WP_ALC%WINDOW 7    ;  @   a   MPI_MEMORY_DEALLOCATE_REAL_WP_ALC%COMM 2   `;  t       MPI_MEMORY_DEALLOCATE_REAL_EP_ALC 8   Ô;     a   MPI_MEMORY_DEALLOCATE_REAL_EP_ALC%ARRAY 8   `<  @   a   MPI_MEMORY_DEALLOCATE_REAL_EP_ALC%NELEM 9    <  @   a   MPI_MEMORY_DEALLOCATE_REAL_EP_ALC%WINDOW 7   à<  @   a   MPI_MEMORY_DEALLOCATE_REAL_EP_ALC%COMM 2    =  t       MPI_MEMORY_DEALLOCATE_REAL_WP_PTR 8   =     a   MPI_MEMORY_DEALLOCATE_REAL_WP_PTR%ARRAY 8    >  @   a   MPI_MEMORY_DEALLOCATE_REAL_WP_PTR%NELEM 9   `>  @   a   MPI_MEMORY_DEALLOCATE_REAL_WP_PTR%WINDOW 7    >  @   a   MPI_MEMORY_DEALLOCATE_REAL_WP_PTR%COMM 2   à>  t       MPI_MEMORY_DEALLOCATE_REAL_EP_PTR 8   T?     a   MPI_MEMORY_DEALLOCATE_REAL_EP_PTR%ARRAY 8   à?  @   a   MPI_MEMORY_DEALLOCATE_REAL_EP_PTR%NELEM 9    @  @   a   MPI_MEMORY_DEALLOCATE_REAL_EP_PTR%WINDOW 7   `@  @   a   MPI_MEMORY_DEALLOCATE_REAL_EP_PTR%COMM 5    @  }       MPI_MEMORY_ALLOCATE_REAL_2DIM_WP_ALC ;   A  ¤   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_WP_ALC%ARRAY <   ÁA  @   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_WP_ALC%NELEM1 <   B  @   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_WP_ALC%NELEM2 :   AB  @   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_WP_ALC%COMM 5   B  }       MPI_MEMORY_ALLOCATE_REAL_2DIM_EP_ALC ;   þB  ¤   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_EP_ALC%ARRAY <   ¢C  @   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_EP_ALC%NELEM1 <   âC  @   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_EP_ALC%NELEM2 :   "D  @   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_EP_ALC%COMM 5   bD  }       MPI_MEMORY_ALLOCATE_REAL_2DIM_WP_PTR ;   ßD  ¤   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_WP_PTR%ARRAY <   E  @   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_WP_PTR%NELEM1 <   ÃE  @   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_WP_PTR%NELEM2 :   F  @   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_WP_PTR%COMM 5   CF  }       MPI_MEMORY_ALLOCATE_REAL_2DIM_EP_PTR ;   ÀF  ¤   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_EP_PTR%ARRAY <   dG  @   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_EP_PTR%NELEM1 <   ¤G  @   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_EP_PTR%NELEM2 :   äG  @   a   MPI_MEMORY_ALLOCATE_REAL_2DIM_EP_PTR%COMM 7   $H  t       MPI_MEMORY_DEALLOCATE_REAL_2DIM_WP_ALC =   H  ¤   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_WP_ALC%ARRAY =   <I  @   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_WP_ALC%NELEM >   |I  @   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_WP_ALC%WINDOW <   ¼I  @   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_WP_ALC%COMM 7   üI  t       MPI_MEMORY_DEALLOCATE_REAL_2DIM_EP_ALC =   pJ  ¤   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_EP_ALC%ARRAY =   K  @   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_EP_ALC%NELEM >   TK  @   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_EP_ALC%WINDOW <   K  @   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_EP_ALC%COMM 7   ÔK  t       MPI_MEMORY_DEALLOCATE_REAL_2DIM_WP_PTR =   HL  ¤   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_WP_PTR%ARRAY =   ìL  @   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_WP_PTR%NELEM >   ,M  @   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_WP_PTR%WINDOW <   lM  @   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_WP_PTR%COMM 7   ¬M  t       MPI_MEMORY_DEALLOCATE_REAL_2DIM_EP_PTR =    N  ¤   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_EP_PTR%ARRAY =   ÄN  @   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_EP_PTR%NELEM >   O  @   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_EP_PTR%WINDOW <   DO  @   a   MPI_MEMORY_DEALLOCATE_REAL_2DIM_EP_PTR%COMM !   O  H       MPI_MEMORY_SETUP 1   ÌO  }       MPI_MEMORY_ALLOCATE_SHARED_BYTES <   IP  @   a   MPI_MEMORY_ALLOCATE_SHARED_BYTES%ALLOC_SIZE ;   P  @   a   MPI_MEMORY_ALLOCATE_SHARED_BYTES%GROUPCOMM 9   ÉP  S   a   MPI_MEMORY_ALLOCATE_SHARED_BYTES%BASEPTR 5   Q  @   a   MPI_MEMORY_ALLOCATE_SHARED_BYTES%WIN '   \Q  ^       MPI_MEMORY_SYNCHRONIZE .   ºQ  @   a   MPI_MEMORY_SYNCHRONIZE%WINDOW ,   úQ  @   a   MPI_MEMORY_SYNCHRONIZE%COMM $   :R  `       MPI_MEMORY_WIN_FREE +   R  @   a   MPI_MEMORY_WIN_FREE%WINDOW +   ÚR  @   a   MPI_MEMORY_WIN_FREE%IERROR 