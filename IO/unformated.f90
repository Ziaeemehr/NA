program WriteUnformated
DOUBLE PRECISION :: D
INTEGER :: I,J

    I = 1024*1024
    J = -1
    D = 10.0D0
    OPEN(8,FILE='xample.out',STATUS='NEW',FORM='UNFORMATTED')
    WRITE(8)I,J,D
    CLOSE(8,STATUS='KEEP')
    STOP 'End of program'
END program

