PROGRAM READ_CG

  IMPLICIT NONE
  INTEGER :: it
  REAL(8) :: x(16384), y(16384), u(16384)
  CHARACTER(15) :: header

  OPEN (100, FILE='CG_result.plt', STATUS='OLD',FORM='UNFORMATTED',access='stream')

  READ(100) header
  print*,header
  DO it= 1,16384
   READ(100,END=999) x(it),y(it),u(it)
   print*,x(it),y(it),u(it)
  ENDDO

999 CONTINUE

END PROGRAM READ_CG
