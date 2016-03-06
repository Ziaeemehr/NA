program ReadUnformated
DOUBLE PRECISION :: D
INTEGER :: I,J

    OPEN(8,FILE='xample.out',STATUS='OLD',FORM='UNFORMATTED')
    READ(8)I,J,D
    CLOSE(8,STATUS='KEEP')
    print *, I,J,D
END program

