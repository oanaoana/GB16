!*****************************************************************************!
!* Copyright (c) 2015-2016, Intel Corporation                                *!
!* All rights reserved.                                                      *!
!*                                                                           *!
!* Redistribution and use in source and binary forms, with or without        *!
!* modification, are permitted provided that the following conditions        *!
!* are met:                                                                  *!
!* 1. Redistributions of source code must retain the above copyright         *!
!*    notice, this list of conditions and the following disclaimer.          *!
!* 2. Redistributions in binary form must reproduce the above copyright      *!
!*    notice, this list of conditions and the following disclaimer in the    *!
!*    documentation and/or other materials provided with the distribution.   *!
!* 3. Neither the name of the copyright holder nor the names of its          *!
!*    contributors may be used to endorse or promote products derived        *!
!*    from this software without specific prior written permission.          *!
!*                                                                           *!
!* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       *!
!* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         *!
!* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     *!
!* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      *!
!* HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    *!
!* SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  *!
!* TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    *!
!* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *!
!* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *!
!* NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *!
!* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *!
!*****************************************************************************!
!* Hans Pabst (Intel Corp.), Alexander Heinecke (Intel Corp.)                *!
!*****************************************************************************!

PROGRAM nekbox_libxsmm_bench
  USE :: LIBXSMM
  USE :: MPI
  IMPLICIT NONE

  INTEGER, PARAMETER :: T = KIND(0D0)

  REAL(T), ALLOCATABLE, TARGET :: a(:,:,:), b(:,:), c(:,:)
  REAL(T), ALLOCATABLE, TARGET :: bufa(:), bufb(:)
  !DIR$ ATTRIBUTES ALIGN:LIBXSMM_ALIGNMENT :: a, b, c
  TYPE(LIBXSMM_DMMFUNCTION) :: xmm
  DOUBLE PRECISION :: duration, gflops, gbytes, sysgflops, scale
  DOUBLE PRECISION :: mingflops, maxgflops, avggflops
  DOUBLE PRECISION :: mingbytes, maxgbytes, avggbytes
  DOUBLE PRECISION :: minar, maxar, avgar
  INTEGER(8) :: i, r, reps, start
  INTEGER :: m, n, k, s, mpierror, rank, procs, arreps, arwarm
  INTEGER :: s_low, s_high, s_step, s_loop

  CALL MPI_Init( mpierror )
  CALL MPI_Comm_size ( MPI_COMM_WORLD, procs, mpierror )
  CALL MPI_Comm_rank ( MPI_COMM_WORLD, rank,  mpierror )

  m = 32
  n = 32
  k = 32

  s_low = 0
  s_high = 768
  s_step = 8

  arreps = 10000
  arwarm = 100

  CALL libxsmm_init()
  CALL libxsmm_dispatch(xmm, m, n, k)

  
  ALLOCATE(bufa(1024))
  ALLOCATE(bufb(1024))
  bufa(:) = 1.0
  bufb(:) = 1.0


  do s_loop = s_low, s_high, s_step
    s = max(1, s_loop)
     
    reps = 6000000/s
    ALLOCATE(a(m,k,s))
    ALLOCATE(b(k,n))
    ALLOCATE(c(m,n))
    
    scale = (1D0 / s)
    DO r = 1, s
      CALL init(42, a(:,:,r), scale, r - 1)  
    END DO
    CALL init(24, b(:,:), scale, 257)
    CALL init(55, c(:,:), scale, 258)
 
    ! Measure GEMM
 
    CALL MPI_Barrier( MPI_COMM_WORLD, mpierror )
    start = libxsmm_timer_tick()
    DO i = 1, reps
      DO r = 1, s
        CALL libxsmm_call(xmm, a(:,:,r), b(:,:), c(:,:))
      END DO
    END DO
    duration = libxsmm_timer_duration(start, libxsmm_timer_tick())
    CALL MPI_Barrier( MPI_COMM_WORLD, mpierror )
 
    IF (0.LT.duration) THEN
      gflops = (2D0 * reps * m * n * k * s * 1D-9 / duration)
      gbytes = (8D0 * reps * m * k * s * 1D-9 / duration)
 
      CALL MPI_Reduce(gflops, mingflops, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD, mpierror)
      CALL MPI_Reduce(gflops, maxgflops, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, mpierror)
      CALL MPI_Reduce(gflops, avggflops, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, mpierror)
 
      CALL MPI_Reduce(gbytes, mingbytes, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD, mpierror)
      CALL MPI_Reduce(gbytes, maxgbytes, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, mpierror)
      CALL MPI_Reduce(gbytes, avggbytes, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, mpierror)
 
      avggflops = avggflops/procs
      avggbytes = avggbytes/procs
 
      CALL MPI_Allreduce( gflops, sysgflops, 1, MPI_DOUBLE, &
        MPI_SUM, MPI_COMM_WORLD, mpierror ) 
 
      IF (0.eq.rank) THEN
        write(*,"(A,I5,A,I5,A,I5,A,I5)") "m=", m, "    n=", n, "    k=", k, "    s=", s
        IF (s.gt.200) THEN
          WRITE(*, "(A,F4.1,A,A,F4.1,A)") &
            "xsmm min performance: ", mingflops, " GFLOPS", &
            "xsmm min bandwidth  : ", mingbytes, " GB/s"
          WRITE(*, "(A,F4.1,A,A,F4.1,A)") &
            "xsmm max performance: ", maxgflops, " GFLOPS", &
            "xsmm max bandwidth  : ", maxgbytes, " GB/s"
          WRITE(*, "(A,F4.1,A,A,F4.1,A)") &
            "xsmm avg performance: ", avggflops, " GFLOPS", &
            "xsmm avg bandwidth  : ", avggbytes, " GB/s"
        ELSE
          WRITE(*, "(A,F4.1,A)") &
            "xsmm min performance: ", mingflops, " GFLOPS"
          WRITE(*, "(A,F4.1,A)") &
            "xsmm max performance: ", maxgflops, " GFLOPS"
          WRITE(*, "(A,F4.1,A)") &
            "xsmm avg performance: ", avggflops, " GFLOPS"
        END IF
        WRITE(*, "(A,F12.1,A,I5)") &
           "system xsmm performance: ", sysgflops, " GFLOPS for s=", s
      ENDIF
    ENDIF 
 
    DEALLOCATE(a)
    DEALLOCATE(b)
    DEALLOCATE(c)

  enddo

  CALL MPI_Barrier( MPI_COMM_WORLD, mpierror )

  ! Measure allreduce
  s = 1
  DO r = 1, 11  
    IF (0.eq.rank) THEN
      WRITE(*, "(A,I10,A)") &
        "allreduce test for: ", s, " double"
    ENDIF

    CALL MPI_Barrier( MPI_COMM_WORLD, mpierror )
    DO i = 1, arwarm
      CALL MPI_Allreduce( bufa, bufb, s, MPI_DOUBLE, &
        MPI_SUM, MPI_COMM_WORLD, mpierror )
      CALL MPI_Barrier( MPI_COMM_WORLD, mpierror )
    END DO

    start = libxsmm_timer_tick()
    DO i = 1, arreps
      CALL MPI_Allreduce( bufa, bufb, s, MPI_DOUBLE, &
        MPI_SUM, MPI_COMM_WORLD, mpierror )
      CALL MPI_Barrier( MPI_COMM_WORLD, mpierror )
    END DO
    duration = libxsmm_timer_duration(start, libxsmm_timer_tick())
    duration = duration / arreps
    
    CALL MPI_Reduce(duration, minar, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD, mpierror)
    CALL MPI_Reduce(duration, maxar, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD, mpierror)
    CALL MPI_Reduce(duration, avgar, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD, mpierror)
    avgar = avgar / procs

    IF (0.eq.rank) THEN
      WRITE(*, "(A,F20.10,A)") &
        "allreduce min time: ", minar*1000000.0, " [us]"
      WRITE(*, "(A,F20.10,A)") &
        "allreduce max time: ", maxar*1000000.0, " [us]"
      WRITE(*, "(A,F20.10,A)") &
        "allreduce avg time: ", avgar*1000000.0, " [us]"
    ENDIF

    CALL MPI_Barrier( MPI_COMM_WORLD, mpierror )
    s = s*2
  END DO

!  do s = 0, procs-1
!    if (rank == s) then
!      write(*,*) rank, rank/32, gflops
!    endif
!    CALL MPI_Barrier( MPI_COMM_WORLD, mpierror )
!  enddo
    

  CALL libxsmm_finalize()
  CALL MPI_Finalize( mpierror )
CONTAINS
  PURE SUBROUTINE init(seed, matrix, scale, n)
    INTEGER, INTENT(IN) :: seed
    REAL(T), INTENT(OUT) :: matrix(:,:)
    REAL(8), INTENT(IN) :: scale
    INTEGER(8), INTENT(IN), OPTIONAL :: n
    INTEGER(8) :: minval, addval, maxval
    INTEGER :: ld, i, j
    REAL(8) :: value, norm
    ld = UBOUND(matrix, 1) - LBOUND(matrix, 1) + 1
    minval = MERGE(n, 0_8, PRESENT(n)) + seed
    addval = (UBOUND(matrix, 1) - LBOUND(matrix, 1)) * ld + (UBOUND(matrix, 2) - LBOUND(matrix, 2))
    maxval = MAX(ABS(minval), addval)
    norm = MERGE(scale / maxval, scale, 0.NE.maxval)
    DO j = LBOUND(matrix, 2), UBOUND(matrix, 2)
      DO i = LBOUND(matrix, 1), LBOUND(matrix, 1) + UBOUND(matrix, 1) - 1
        value = (i - LBOUND(matrix, 1)) * ld + (j - LBOUND(matrix, 2)) + minval
        matrix(i,j) = norm * (value - 0.5D0 * addval)
      END DO
    END DO
  END SUBROUTINE
END PROGRAM
