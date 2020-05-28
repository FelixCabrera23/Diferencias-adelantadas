! 2020 - 19 - 5
! parabolicas.f90
! Félix Cabrera (walberto.cabrera@gmail.com)

! Solución de Ecuaciónes parciales Parabolicas
! SOLUCION DE LA ECUACIÓN DE CALOR
! POR EL METODO DE DIFERENCIAS REGRECIVAS

! Codificación del texto UTF-8
! Compiladores probados: GNU Fortran (Ubuntu 9.2.1-9ubuntu2) 9.2.1 2019008

! Algoritmo tomado del libro de analisis numerico de Burden
! Algoritmo 12.2 pgs 557-558

! Requiere:

! instrucciónes de compilación:
! gfortran -Wall -pedantic -std=f95 -c Fx.f90
! gfortran -Wall -pedantic -std=f95 -c funciones.f90
! gfortran -Wall -pedantic -std=f95 -c parabolicas.f90
! gfortran -Wall -pedantic -std=f95 -o calor parabolicas.o funciones.o Fx.o
! ./calor

PROGRAM parabolica
  USE funciones
  IMPLICIT NONE
  
  ! Definimos variables principales
  INTEGER(8):: m,N
  REAL(8):: l_s,T ! extremo l y tiempo máx T
  REAL(8), ALLOCATABLE:: x(:),w(:)
  ! Variables auxiliares
  REAL(8):: h,k,alpha,lambda
  INTEGER(8):: i,j
  REAL(8), ALLOCATABLE:: l(:),u(:),z(:)
  ! Variables de control
  INTEGER(4):: err

  ! Abrimos el archivo de configuración (12)
  OPEN (12, FILE='parametros.config', STATUS='old', IOSTAT=err)
  IF (err .ne. 0) STOP 'parametros.config is missing'

  READ(12,*) l_s
  READ(12,*) T
  READ(12,*) alpha
  READ(12,*) m
  READ(12,*) N
  CLOSE(12)
  
  ALLOCATE(w(m),x(m),l(m),u(m),z(m))
  
  ! Paso 1
  h = l_s/REAL(m,8)
  k = T/REAL(N,8)
  lambda =((alpha**2.)*k)/(h**2.)

  ! Paso 2  ! Definimos los valores iniciales
  DO i=1, m-1
    w(i) = Fx(i*h)
  END DO
  
  ! Paso 3
  l(1)=1.+(2.*lambda)
  u(1)=(-lambda)/l(1)

  ! Paso 4
  DO i=2, M-2
    l(i)=1.+(2.*lambda)+(lambda*u(i-1))
    u(i)=(-lambda)/l(i)
  END DO
    
  ! Paso 5
  l(m-1)=1.+(2.*lambda)+(lambda*u(m-2))
  
  ! Paso 6
  DO j=1, N ! para los pasos 7 -11
    ! Paso 7
    t=j*k
    z(1)=w(1)/l(1)
    ! Paso 8
    DO i=2, m-1
      z(i)=(w(i)+(lambda*z(i-1)))/l(i)
    END DO
    ! Paso 9
    w(m-1)=z(m-1)
    ! Paso 10
    i=m-2
    DO WHILE (i>=1)
      w(i)=z(i)-(u(i)*w(i+1))
      i=i-1
    END DO
    ! Paso 11
    IF (j==N) THEN
      DO i=1, m-1
        x=i*h
        PRINT *, x(i),t,w(i)
      END DO
    END IF
  END DO
  ! Paso 12 Fin
  
  DEALLOCATE(w,x,l,u,z)
    
END PROGRAM parabolica















