NAME min_flip_model
ROWS
 N  OBJ
 L  R0      
 L  R1      
 L  R2      
 L  R3      
 L  R4      
 L  R5      
 L  R6      
 L  R7      
 L  R8      
 L  R9      
 L  R10     
 L  R11     
 L  R12     
 L  R13     
 L  R14     
 G  R15     
COLUMNS
    MARKER    'MARKER'                 'INTORG'
    X[0,0]    OBJ       1
    X[0,0]    R0        -1
    X[0,0]    R3        1
    X[0,0]    R6        1
    X[0,0]    R11       0.5
    X[0,1]    R0        1
    X[0,1]    R3        -1
    X[0,1]    R6        1
    X[0,1]    R10       -1
    X[0,1]    R13       0.5
    X[1,0]    R1        -1
    X[1,0]    R4        1
    X[1,0]    R7        1
    X[1,0]    R10       -1
    X[1,0]    R12       0.5
    X[1,1]    OBJ       1
    X[1,1]    R1        1
    X[1,1]    R4        -1
    X[1,1]    R7        1
    X[1,1]    R14       0.5
    X[2,0]    R2        -1
    X[2,0]    R5        1
    X[2,0]    R8        1
    X[2,0]    R10       -1
    X[2,0]    R11       -0.5
    X[2,0]    R12       -0.5
    X[2,1]    R2        1
    X[2,1]    R5        -1
    X[2,1]    R8        1
    X[2,1]    R10       -1
    X[2,1]    R13       -0.5
    X[2,1]    R14       -0.5
    B01[0,0]  OBJ       0
    B01[0,1]  R0        -1
    B01[0,1]  R1        -1
    B01[0,1]  R2        -1
    B01[0,1]  R9        1
    B01[1,0]  OBJ       0
    B01[1,1]  OBJ       0
    B10[0,0]  OBJ       0
    B10[0,1]  R3        -1
    B10[0,1]  R4        -1
    B10[0,1]  R5        -1
    B10[0,1]  R9        1
    B10[1,0]  OBJ       0
    B10[1,1]  OBJ       0
    B11[0,0]  OBJ       0
    B11[0,1]  R6        -1
    B11[0,1]  R7        -1
    B11[0,1]  R8        -1
    B11[0,1]  R9        1
    B11[1,0]  OBJ       0
    B11[1,1]  OBJ       0
    z[0,0]    R11       1
    z[0,0]    R15       1
    z[0,1]    R12       1
    z[0,1]    R15       1
    z[1,0]    R13       1
    z[1,0]    R15       1
    z[1,1]    R14       1
    z[1,1]    R15       1
    MARKER    'MARKER'                 'INTEND'
RHS
    RHS1      R9        2
    RHS1      R10       -2
    RHS1      R11       0.5
    RHS1      R12       0.5
    RHS1      R13       0.5
    RHS1      R14       0.5
    RHS1      R15       1
BOUNDS
 BV BND1      X[0,0]  
 BV BND1      X[0,1]  
 BV BND1      X[1,0]  
 BV BND1      X[1,1]  
 BV BND1      X[2,0]  
 BV BND1      X[2,1]  
 BV BND1      B01[0,0]
 BV BND1      B01[0,1]
 BV BND1      B01[1,0]
 BV BND1      B01[1,1]
 BV BND1      B10[0,0]
 BV BND1      B10[0,1]
 BV BND1      B10[1,0]
 BV BND1      B10[1,1]
 BV BND1      B11[0,0]
 BV BND1      B11[0,1]
 BV BND1      B11[1,0]
 BV BND1      B11[1,1]
 BV BND1      z[0,0]  
 BV BND1      z[0,1]  
 BV BND1      z[1,0]  
 BV BND1      z[1,1]  
ENDATA
